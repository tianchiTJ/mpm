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
    // Properties
    properties_ = material_properties;
    // Cam Clay parameters
    // Reference mean pressure
    p_ref_ = material_properties["p_ref"].template get<double>();
    // Reference void ratio
    e_ref_ = material_properties["e_ref"].template get<double>();
    // Initial preconsolidation pressure
    pc0_ = material_properties["pc0"].template get<double>();
    // OCR
    ocr_ = material_properties["ocr"].template get<double>();
    // M
    m_ = material_properties["m"].template get<double>();
    // Lambda
    lambda_ = material_properties["lambda"].template get<double>();
    // Kappa
    kappa_ = material_properties["kappa"].template get<double>();
    // Ellipticity
    ellipticity_ = material_properties["ellipticity"].template get<double>();
    // Initial void ratio
    e0_ = e_ref_ - lambda_ * log(pc0_ / p_ref_) - kappa_ * log(ocr_);
  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::CamClay<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {
      // Elastic modulus
      // Bulk modulus
      {"bulk_modulus", youngs_modulus_},
      // Shear modulus
      {"shear_modulus", 3 * youngs_modulus_ * (1 - 2 * poisson_ratio_) /
                            (2 * (1 + poisson_ratio_))},
      // Stress invariants
      // Volumetric stress
      {"p", 0.},
      // Deviatoric stress
      {"q", 0.},
      // Lode's angle
      {"theta", 0.},
      // Lode's angle multiplier
      {"lode_multiplier", 1.},
      // Cam clay parameters
      // Preconsolidation pressure
      {"pc", pc0_},
      // void_ratio
      {"void_ratio", e0_},
      // Consistency multiplier
      {"multiplier", 0.},
      {"f_function", 0.}};

  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_elastic_tensor(mpm::dense_map* state_vars) {
  // Bulk modulus
  if ((*state_vars).at("p") > 0) {
    (*state_vars).at("bulk_modulus") =
        (1 + (*state_vars).at("void_ratio")) / kappa_ * (*state_vars).at("p");
    // Shear modulus
    (*state_vars).at("shear_modulus") = 3 * (*state_vars).at("bulk_modulus") *
                                        (1 - 2 * poisson_ratio_) /
                                        (2 * (1 + poisson_ratio_));
  }
  // Components in stiffness matrix
  const double G = (*state_vars).at("shear_modulus");
  const double a1 = (*state_vars).at("bulk_modulus") + (4.0 / 3.0) * G;
  const double a2 = (*state_vars).at("bulk_modulus") - (2.0 / 3.0) * G;
  // Compute elastic stiffness matrix
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

//! Compute plastic tensor
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_plastic_tensor(const Vector6d& stress,
                                                mpm::dense_map* state_vars) {
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  // Preconsolidation pressure
  const double pc = (*state_vars).at("pc");
  // Compute dF / dp
  const double df_dp = 2 * p - pc;
  // Compute dF / dq
  const double df_dq = 2 * q / (m_ * m_);
  // Upsilon
  const double upsilon =
      (1 + (*state_vars).at("void_ratio")) / (lambda_ - kappa_);
  // chi
  const double chi = (*state_vars).at("bulk_modulus") * (df_dp * df_dp) +
                     3 * (*state_vars).at("shear_modulus") * (df_dq * df_dq) +
                     upsilon * p * pc * df_dp;
  // Coefficients in plastic stiffness matrix
  const double a1 = pow(((*state_vars).at("bulk_modulus") * df_dp), 2);
  const double a2 = -sqrt(6) * (*state_vars).at("bulk_modulus") * df_dp *
                    (*state_vars).at("shear_modulus") * df_dq;
  const double a3 = 6 * pow(((*state_vars).at("shear_modulus") * df_dq), 2);
  // Compute the deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) + p;
  dev_stress(1) = stress(1) + p;
  dev_stress(2) = stress(2) + p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }
  Matrix6x6 n_l = Matrix6x6::Zero();
  Matrix6x6 l_n = Matrix6x6::Zero();
  Matrix6x6 l_l = Matrix6x6::Zero();
  Matrix6x6 n_n = Matrix6x6::Zero();
  double xi = q / sqrt(1.5);
  // lxn
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) l_n(i, j) = dev_stress(j) / xi;
  }
  // nxl
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) n_l(i, j) = dev_stress(i) / xi;
  }
  // nxn
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      n_n(i, j) = dev_stress(i) * dev_stress(j) / (xi * xi);
      if (i > 2 && j > 2) n_n(i, j) /= 2;
    }
  }
  // lxl
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) l_l(i, j) = 1;
  }
  // Compute plastic tensor
  this->dp_ = 1. / chi * (a1 * l_l + a2 * (n_l + l_n) + a3 * (n_n));

  return true;
}

//! Compute stress invariants
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_stress_invariants(const Vector6d& stress,
                                                   mpm::dense_map* state_vars) {
  // Compute volumetic stress
  (*state_vars).at("p") = -(stress(0) + stress(1) + stress(2)) / 3.;
  // Compute the deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) + (*state_vars).at("p");
  dev_stress(1) = stress(1) + (*state_vars).at("p");
  dev_stress(2) = stress(2) + (*state_vars).at("p");
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }
  // Compute J2
  double j2 =
      (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
       pow((stress(0) - stress(2)), 2)) /
          6.0 +
      pow(stress(3), 2);
  if (Tdim == 3) j2 += pow(stress(4), 2) + pow(stress(5), 2);
  // Compute q
  (*state_vars)["q"] = sqrt(3 * j2);
  // Compute J3
  double j3 = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
              (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    j3 += ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
           (dev_stress(0) * pow(dev_stress(4), 2)) -
           (dev_stress(1) * pow(dev_stress(5), 2)));
  // Compute theta value (Lode angle)
  double theta_val = 0.;
  if (fabs(j2) > 0.0) theta_val = (3. * sqrt(3.) / 2.) * (j3 / pow(j2, 1.5));
  // Check theta value
  if (theta_val > 1.0) theta_val = 1.0;
  if (theta_val < -1.0) theta_val = -1.0;
  // Compute theta
  (*state_vars).at("theta") = (1. / 3.) * acos(theta_val);
  // Check theta
  if ((*state_vars).at("theta") > M_PI / 3.)
    (*state_vars).at("theta") = M_PI / 3.;
  if ((*state_vars).at("theta") < 0.0) (*state_vars).at("theta") = 0.;

  return true;
}

//! Compute lode angle effect
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_lode_multiplier(mpm::dense_map* state_vars) {
  // Get stress invariants
  const double theta = (*state_vars).at("theta");
  (*state_vars).at("lode_multiplier") = 1.;
  double num = 4 * (1 - pow(ellipticity_, 2)) * pow(cos(theta), 2);
  double root = num + 5 * pow(ellipticity_, 2) - 4 * ellipticity_;
  if (root > 0) {
    double den = 2 * (1 - pow(ellipticity_, 2)) * cos(theta) +
                 (2 * ellipticity_ - 1) * sqrt(root);
    if (fabs(den) > 1.E-10) {
      num += pow((2 * ellipticity_ - 1), 2);
      (*state_vars).at("lode_multiplier") = num / den;
    }
  }
  return true;
}

//! Compute yield function and yield state
template <unsigned Tdim>
typename mpm::CamClay<Tdim>::FailureState
    mpm::CamClay<Tdim>::compute_yield_state(double* yield_function,
                                            const mpm::dense_map* state_vars) {
  // Get stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  const double theta = (*state_vars).at("theta");
  const double lode_multiplier = (*state_vars).at("lode_multiplier");
  // Plastic volumetic strain
  const double pc = (*state_vars).at("pc");
  // Initialise yield status (0: elastic, 1: yield)
  auto yield_type = FailureState::Elastic;
  // Compute yield functions
  (*yield_function) =
      q * q * lode_multiplier * lode_multiplier / (m_ * m_) + p * (p - pc);
  // Tension failure
  if ((*yield_function) > 1.E-22) yield_type = FailureState::Yield;

  return yield_type;
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
  // Lode angle multiplier
  const double lode_multiplier = (*state_vars).at("lode_multiplier");
  // Compute dF / dp
  double df_dp = 2 * p - pc;
  // Compute dF / dq
  double df_dq = 2 * q / (m_ * m_) * lode_multiplier * lode_multiplier;
  // Compute dF / dpc
  double df_dpc = -p;
  // Upsilon
  double upsilon = (1 + (*state_vars).at("void_ratio")) / (lambda_ - kappa_);
  // A_den
  double a_den = 1 + (2 * (*state_vars).at("bulk_modulus") + upsilon * pc) *
                         (*state_vars).at("multiplier");
  // Compute dp / dmultiplier
  double dp_dmul = -(*state_vars).at("bulk_modulus") * (2 * p - pc) / a_den;
  // Compute dpc / dmultiplier
  double dpc_dmul = upsilon * pc * (2 * p - pc) / a_den;
  // Compute dq / dmultiplier
  double dq_dmul = -q / ((*state_vars).at("multiplier") +
                         m_ * m_ /
                             (6 * (*state_vars).at("shear_modulus") *
                              lode_multiplier * lode_multiplier));
  // Compute dF / dmultiplier
  (*df_dmul) = (df_dp * dp_dmul) + (df_dq * dq_dmul) + (df_dpc * dpc_dmul);
}

//! Compute dg/dpc
template <unsigned Tdim>
void mpm::CamClay<Tdim>::compute_dg_dpc(const mpm::dense_map* state_vars,
                                        const double pc_n, const double p_trial,
                                        double* g_function, double* dg_dpc) {
  // Upsilon
  const double upsilon =
      (1 + (*state_vars).at("void_ratio")) / (lambda_ - kappa_);
  // Exponential index
  double e_index = upsilon * (*state_vars).at("multiplier") *
                   (2 * p_trial - (*state_vars).at("pc")) /
                   (1 + 2 * (*state_vars).at("multiplier") *
                            (*state_vars).at("bulk_modulus"));
  // Compute consistency multiplier function
  (*g_function) = pc_n * exp(e_index) - (*state_vars).at("pc");
  // Compute dG / dpc
  (*dg_dpc) = pc_n * exp(e_index) *
                  (-upsilon * (*state_vars).at("multiplier") /
                   (1 + 2 * (*state_vars).at("multiplier") *
                            (*state_vars).at("bulk_modulus"))) -
              1;
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
  (*state_vars).at("p") = -(stress(0) + stress(1) + stress(2)) / 3.;
  // Set elastic tensor
  this->compute_elastic_tensor(state_vars);
  //-------------------------------------------------------------------------
  // Elastic step
  // Compute trial stress
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  // Compute trial stress invariants
  this->compute_stress_invariants(trial_stress, state_vars);
  // Compute lode multiplier
  if (ellipticity_ != 1) this->compute_lode_multiplier(state_vars);
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
  // Preconsolidation pressure of last step
  const double pc_n = (*state_vars).at("pc");
  // Initialise dF / dmul
  double df_dmul = 0;
  // Iteration for consistency multiplier
  while (fabs(f_function) > Ftolerance && counter_f < itrstep) {
    // Compute dF / dmul
    this->compute_df_dmul(state_vars, &df_dmul);
    // Update consistency multiplier
    // if (fabs(df_dmul) > 1.E-5)
    (*state_vars).at("multiplier") -= f_function / df_dmul;
    // Initialise G and dG / dpc
    double g_function = 0;
    double dg_dpc = 0;
    // Compute G and dG / dpc
    this->compute_dg_dpc(state_vars, pc_n, p_trial, &g_function, &dg_dpc);
    // Subiteraction for preconsolidation pressure
    while (fabs(g_function) > Gtolerance && counter_g < substep) {
      // Update preconsolidation pressure
      (*state_vars).at("pc") -= g_function / dg_dpc;
      // Update G and dG / dpc
      this->compute_dg_dpc(state_vars, pc_n, p_trial, &g_function, &dg_dpc);
      // Counter subiteration step
      ++counter_g;
    }
    // Update volumetric stress
    (*state_vars).at("p") = (p_trial + (*state_vars).at("bulk_modulus") *
                                           (*state_vars).at("multiplier") *
                                           (*state_vars).at("pc")) /
                            (1 + 2 * (*state_vars).at("bulk_modulus") *
                                     (*state_vars).at("multiplier"));
    // Update deviatoric stress
    (*state_vars).at("q") =
        q_trial / (1 + 6 * (*state_vars).at("shear_modulus") *
                           (*state_vars).at("multiplier") *
                           (*state_vars).at("lode_multiplier") *
                           (*state_vars).at("lode_multiplier") / (m_ * m_));

    // Update yield function
    yield_type = this->compute_yield_state(&f_function, state_vars);
    // Counter iteration step
    ++counter_f;
    // Update lode angle
    if (ellipticity_ != 1) {
      // Set plastic tensor
      this->compute_plastic_tensor(stress, state_vars);
      // Modified dstrain
      Vector6d dstrain_m = dstrain;
      dstrain_m(3) *= 0.5;
      dstrain_m(4) *= 0.5;
      dstrain_m(5) *= 0.5;
      // Compte current stress
      Vector6d updated_stress = trial_stress - this->dp_ * dstrain_m;
      // Compute stress invariants
      this->compute_stress_invariants(updated_stress, state_vars);
      // Compute lode multiplier
      this->compute_lode_multiplier(state_vars);
    }
  }

  // Set plastic tensor
  this->compute_plastic_tensor(stress, state_vars);
  // Check yield function
  (*state_vars).at("f_function") = f_function;
  // Modified dstrain
  Vector6d dstrain_m = dstrain;
  dstrain_m(3) *= 0.5;
  dstrain_m(4) *= 0.5;
  dstrain_m(5) *= 0.5;
  // Compte current stress
  Vector6d updated_stress = trial_stress - this->dp_ * dstrain_m;
  // Incremental Volumetric strain
  double dvstrain = dstrain(0) + dstrain(1) + dstrain(2);
  // Update void_ratio
  (*state_vars).at("void_ratio") += dvstrain * (1 + e0_);

  return updated_stress;
}
