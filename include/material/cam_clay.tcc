//! Constructor with id and material properties
template <unsigned Tdim>
mpm::CamClay<Tdim>::CamClay(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // General parameters
    // Density
    density_ = material_properties["density"].template get<double>();
    // Young's modulus
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    // Initial porosity
    e0_ = material_properties["e0"].template get<double>();
    // Properties
    properties_ = material_properties;
    // Cam Clay parameters
    // M
    m_value_ = material_properties["m_value"].template get<double>();
    // Initial preconsolidation pressure
    p0_ = material_properties["p0"].template get<double>();
    // Lambda
    lambda_ = material_properties["lambda"].template get<double>();
    // Kappa
    kappa_ = material_properties["kappa"].template get<double>();
  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::CamClay<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// Stress invariants
                               // J2
                               {"j2", 0.},
                               // J3
                               {"j3", 0.},
                               // Volumetric stress
                               {"p", 0.},
                               // Deviatoric stress
                               {"q", 0.},
                               // Lode's angle
                               {"theta", 0.},
                               // Preconsolidation pressure
                               {"pc", p0_},
                               // Porosity
                               {"porosity", e0_},
                               // Consistency multiplier
                               {"multiplier", 0.}};

  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_elastic_tensor(
    const mpm::dense_map* state_vars) {
  // Bulk modulus
  if ((*state_vars).at("p") < 0)
    bulk_modulus_ = fabs((1 + (*state_vars).at("porosity")) / kappa_ *
                         (*state_vars).at("p"));
  else
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  // Shear modulus
  shear_modulus_ =
      3 * bulk_modulus_ * (1 - 2 * poisson_ratio_) / (2 * (1 + poisson_ratio_));
  // Components in stiffness matrix
  const double G = shear_modulus_;
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  // compute elastic stiffness matrix
  // clang-format off
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on
  return true;
}

//! Compute stress invariants
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_stress_invariants(const Vector6d& stress,
                                                   mpm::dense_map* state_vars) {
  // Compute volumetic stress
  (*state_vars).at("p") = (stress(0) + stress(1) + stress(2)) / 3.;
  // Compute the deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - (*state_vars).at("p");
  dev_stress(1) = stress(1) - (*state_vars).at("p");
  dev_stress(2) = stress(2) - (*state_vars).at("p");
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }
  // Compute J2
  (*state_vars)["j2"] =
      (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
       pow((stress(0) - stress(2)), 2)) /
          6.0 +
      pow(stress(3), 2);
  if (Tdim == 3) (*state_vars)["j2"] += pow(stress(4), 2) + pow(stress(5), 2);
  // Compute J3
  (*state_vars)["j3"] = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
                        (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    (*state_vars)["j3"] +=
        ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
         (dev_stress(0) * pow(dev_stress(4), 2)) -
         (dev_stress(1) * pow(dev_stress(5), 2)));
  // Compute theta value (Lode angle)
  double theta_val = 0.;
  if (fabs((*state_vars).at("j2")) > 0.0)
    theta_val = (3. * sqrt(3.) / 2.) *
                ((*state_vars).at("j3") / pow((*state_vars).at("j2"), 1.5));
  // Check theta value
  if (theta_val > 1.0) theta_val = 1.0;
  if (theta_val < -1.0) theta_val = -1.0;
  // Compute theta
  (*state_vars)["theta"] = (1. / 3.) * acos(theta_val);
  // Check theta
  if ((*state_vars).at("theta") > M_PI / 3.) (*state_vars)["theta"] = M_PI / 3.;
  if ((*state_vars).at("theta") < 0.0) (*state_vars)["theta"] = 0.;
  // Compute q
  (*state_vars)["q"] = sqrt(3 * ((*state_vars).at("j2")));

  return true;
}

//! Compute yield function and yield state
template <unsigned Tdim>
typename mpm::CamClay<Tdim>::FailureState
    mpm::CamClay<Tdim>::compute_yield_state(double* yield_function,
                                            const mpm::dense_map* state_vars) {
  // Tolerance for yield function
  const double Tolerance = 1.E-5;
  // Get stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  // Plastic volumetic strain
  const double pc = (*state_vars).at("pc");
  // Initialise yield status (0: elastic, 1: yield)
  auto yield_type = FailureState::Elastic;
  // Compute yield functions
  (*yield_function) = q * q / (m_value_ * m_value_) + p * (p - pc);
  // Tension failure
  if ((*yield_function) > Tolerance) yield_type = FailureState::Yield;
  return yield_type;
}

//! Compute dF/dSigma and dY/dSigma
template <unsigned Tdim>
void mpm::CamClay<Tdim>::compute_df_dy(const mpm::dense_map* state_vars,
                                       const Vector6d& stress,
                                       Vector6d* df_dsigma,
                                       Vector6d* dy_dsigma) {
  // Stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  // Preconsolidation pressure
  const double pc = (*state_vars).at("pc");
  // Compute deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - p;
  dev_stress(1) = stress(1) - p;
  dev_stress(2) = stress(2) - p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }
  // Compute dF / dp
  double df_dp = 2 * p - pc;
  // Compute dF / dq
  double df_dq = 2 * q / (m_value_ * m_value_);
  // Compute dF / dpc
  double df_dpc = -p;
  // Compute dp / dSigma
  Vector6d dp_dsigma = Vector6d::Zero();
  dp_dsigma(0) = 1. / 3.;
  dp_dsigma(1) = 1. / 3.;
  dp_dsigma(2) = 1. / 3.;
  // Compute dq / dSigma
  Vector6d dq_dsigma = Vector6d::Zero();
  double multiplier = 1.;
  if (fabs(q) > 0.) multiplier = 1.5 / q;
  dq_dsigma = multiplier * dev_stress;
  if (Tdim == 2) {
    dq_dsigma(4) = 0.;
    dq_dsigma(5) = 0.;
  }
  // Compute dF / dSigma
  (*df_dsigma) = (df_dp * dp_dsigma) + (df_dq * dq_dsigma);
  if (Tdim == 2) {
    (*df_dsigma)(4) = 0.;
    (*df_dsigma)(5) = 0.;
  }
  // Compute dY / dsigma
  (*dy_dsigma) = (*df_dsigma);
}

//! Compute dF/dmul
template <unsigned Tdim>
void mpm::CamClay<Tdim>::compute_df_dmul(const mpm::dense_map* state_vars,
                                         double* df_dmul) {
  // Stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  // Preconsolidation pressure
  const double pc = (*state_vars).at("pc");
  // Compute dF / dp
  double df_dp = 2 * p - pc;
  // Compute dF / dq
  double df_dq = 2 * q / (m_value_ * m_value_);
  // Compute dF / dpc
  double df_dpc = -p;
  // Upsilon
  double upsilon = (1 + (*state_vars).at("porosity")) / (lambda_ - kappa_);
  // A_den
  double a_den =
      1 + (2 * bulk_modulus_ + upsilon * pc) * (*state_vars).at("multiplier");
  // Compute dp / dmultiplier
  double dp_dmul = -bulk_modulus_;
  // Compute dpc / dmultiplier
  double dpc_dmul = upsilon * pc;
  if (fabs(a_den) > 1.E-5) {
    dp_dmul *= ((2 * p - pc) / a_den);
    dpc_dmul *= ((2 * p - pc) / a_den);
  }
  // B_den
  double b_den = (*state_vars).at("multiplier") +
                 m_value_ * m_value_ / (6 * shear_modulus_);
  // Compute dq / dmultiplier
  double dq_dmul = -q;
  if (fabs(b_den) > 1.E-5) dq_dmul / b_den;
  // Compute dF / dmultiplier
  (*df_dmul) = (df_dp * dp_dmul) + (df_dq * dq_dmul) + (df_dpc * dpc_dmul);
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::CamClay<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  // Tolerance for yield function
  const double Ftolerance = 1.E-5;
  // Tolerance for preconsolidation function
  const double Gtolerance = 1.E-5;
  // Maximum iteration step number
  const int itrstep = 100;
  // Maximum subiteration step number
  const int substep = 100;
  // Update volumetric stress
  (*state_vars).at("p") = (stress(0) + stress(1) + stress(2)) / 3.;
  // Set elastic tensor
  this->compute_elastic_tensor(state_vars);
  //-------------------------------------------------------------------------
  // Elastic step
  // Compute trial stress
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  // Compute trial stress invariants
  this->compute_stress_invariants(trial_stress, state_vars);
  // Initialise value for yield function
  double f_function;
  auto yield_type = this->compute_yield_state(&f_function, state_vars);
  // Return the updated stress in elastic state
  if (yield_type == FailureState::Elastic) return trial_stress;
  //-------------------------------------------------------------------------
  // Plastic step
  // Counters for interations
  int counter_f = 0;
  int counter_g = 0;
  // Initialise consistency multiplier
  (*state_vars).at("multiplier") = 0.;
  // Volumetric trial stress
  const double p_trial = (*state_vars).at("p");
  // Deviatoric trial stress
  const double q_trial = (*state_vars).at("q");
  // Upsilon
  const double upsilon =
      (1 + (*state_vars).at("porosity")) / (lambda_ - kappa_);
  // Preconsolidation pressure of last step
  const double pc_n = (*state_vars).at("pc");
  // Compute dF / dmul
  double df_dmul = 0.;
  this->compute_df_dmul(state_vars, &df_dmul);
  // Iteration for consistency multiplier
  while (fabs(f_function) > Ftolerance && counter_f < itrstep) {
    // Update consistency multiplier
    if (fabs(df_dmul) > 1.E-5)
      (*state_vars).at("multiplier") -= f_function / df_dmul;
    // Update volumetric stress
    (*state_vars).at("p") =
        p_trial - bulk_modulus_ * (*state_vars).at("multiplier") *
                      (2 * (*state_vars).at("p") - (*state_vars).at("pc"));
    // Update deviatoric stress
    (*state_vars).at("q") =
        q_trial / (1 + 6 * shear_modulus_ * (*state_vars).at("multiplier") /
                           (m_value_ * m_value_));

    // Exponential index
    double e_index = upsilon * (*state_vars).at("multiplier") *
                     (2 * p_trial - (*state_vars).at("pc")) /
                     (1 + 2 * (*state_vars).at("multiplier") * bulk_modulus_);
    // Compute consistency multiplier function
    double g_function = pc_n * exp(e_index) - (*state_vars).at("pc");
    // Subiteraction for preconsolidation pressure
    while (fabs(g_function) > Gtolerance && counter_g < substep) {
      // Update preconsolidation pressure
      (*state_vars).at("pc") =
          pc_n * exp(upsilon * (*state_vars).at("multiplier") *
                     (2 * (*state_vars).at("p") - (*state_vars).at("pc")));
      // Update exponential index
      e_index = upsilon * (*state_vars).at("multiplier") *
                (2 * p_trial - (*state_vars).at("pc")) /
                (1 + 2 * (*state_vars).at("multiplier") * bulk_modulus_);
      // Update consistency multiplier function
      g_function = pc_n * exp(e_index) - (*state_vars).at("pc");
      // Counter subiteration step
      ++counter_g;
    }
    // Update yield function
    yield_type = this->compute_yield_state(&f_function, state_vars);
    // Compute dF / dmul
    this->compute_df_dmul(state_vars, &df_dmul);
    // Counter iteration step
    ++counter_f;
  }
  // Compute trial stress invariants
  this->compute_stress_invariants(trial_stress, state_vars);
  //! Compute dF/dSigma and dY/dSigma
  Vector6d df_dsigma = Vector6d::Zero();
  Vector6d dy_dsigma = Vector6d::Zero();
  this->compute_df_dy(state_vars, stress, &df_dsigma, &dy_dsigma);
  // Incremental Volumetric strain
  double dvstrain = dstrain(0) + dstrain(1) + dstrain(2);
  // Update porosity
  (*state_vars).at("porosity") += dvstrain * (1 + e0_);
  // Compte current stress
  Vector6d updated_stress =
      trial_stress - ((*state_vars).at("multiplier") * this->de_ * dy_dsigma);

  return updated_stress;
}
