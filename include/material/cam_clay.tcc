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
      // {"theta", 0.},
      // Cam clay parameters
      // Preconsolidation pressure
      {"pc", p0_},
      // Porosity
      {"porosity", e0_},
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
        (1 + (*state_vars).at("porosity")) / kappa_ * (*state_vars).at("p");
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
  const double df_dq = 2 * q / (m_value_ * m_value_);
  // Upsilon
  const double upsilon =
      (1 + (*state_vars).at("porosity")) / (lambda_ - kappa_);
  // chi
  const double chi = (*state_vars).at("bulk_modulus") * (df_dp * df_dp) +
                     3 * (*state_vars).at("shear_modulus") * (df_dq * df_dq) +
                     upsilon * p * pc * df_dp;
  // Coefficients in plastic stiffness matrix
  const double a1 = pow(((*state_vars).at("bulk_modulus") * df_dp), 2);
  const double a2 = -sqrt(1.5) * (*state_vars).at("bulk_modulus") * df_dp *
                    (*state_vars).at("shear_modulus") * df_dq;
  const double a3 = 1.5 * pow(((*state_vars).at("shear_modulus") * df_dq), 2);
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
  // Plastic volumetic strain
  const double pc = (*state_vars).at("pc");
  // Initialise yield status (0: elastic, 1: yield)
  auto yield_type = FailureState::Elastic;
  // Compute yield functions
  (*yield_function) = q * q / (m_value_ * m_value_) + p * (p - pc);
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
  // Compute dF / dp
  double df_dp = 2 * p - pc;
  // Compute dF / dq
  double df_dq = 2 * q / (m_value_ * m_value_);
  // Compute dF / dpc
  double df_dpc = -p;
  // Upsilon
  double upsilon = (1 + (*state_vars).at("porosity")) / (lambda_ - kappa_);
  // A_den
  double a_den = 1 + (2 * (*state_vars).at("bulk_modulus") + upsilon * pc) *
                         (*state_vars).at("multiplier");
  // Compute dp / dmultiplier
  double dp_dmul =
      -(*state_vars).at("bulk_modulus") * (2 * p - pc) / (a_den + 1.E-10);
  // Compute dpc / dmultiplier
  double dpc_dmul = upsilon * pc * (2 * p - pc) / (a_den + 1.E-10);
  // Compute dq / dmultiplier
  double dq_dmul =
      -q /
      ((*state_vars).at("multiplier") +
       m_value_ * m_value_ / (6 * (*state_vars).at("shear_modulus")) + 1.E-10);
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
      (1 + (*state_vars).at("porosity")) / (lambda_ - kappa_);
  // Exponential index
  double e_index =
      upsilon * (*state_vars).at("multiplier") *
      (2 * p_trial - (*state_vars).at("pc")) /
      (1 +
       2 * (*state_vars).at("multiplier") * (*state_vars).at("bulk_modulus") +
       1.E-10);
  // Compute consistency multiplier function
  (*g_function) = pc_n * exp(e_index) - (*state_vars).at("pc");
  // Compute dG / dpc
  (*dg_dpc) = pc_n * exp(e_index) *
                  (-upsilon * (*state_vars).at("multiplier") /
                   (1 +
                    2 * (*state_vars).at("multiplier") *
                        (*state_vars).at("bulk_modulus") +
                    1.E-10)) -
              1;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::CamClay<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  // Tolerance for yield function
  const double Ftolerance = 1.E-1;
  // Tolerance for preconsolidation function
  const double Gtolerance = 1.E-1;
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
      // if (fabs(dg_dpc) > 1.E-5)
      (*state_vars).at("pc") -= g_function / dg_dpc;
      // Update G and dG / dpc
      this->compute_dg_dpc(state_vars, pc_n, p_trial, &g_function, &dg_dpc);
      // Counter subiteration step
      ++counter_g;
    }
    // Update volumetric stress
    (*state_vars).at("p") =
        (p_trial + (*state_vars).at("bulk_modulus") *
                       (*state_vars).at("multiplier") *
                       (*state_vars).at("pc")) /
        (1 +
         2 * (*state_vars).at("bulk_modulus") * (*state_vars).at("multiplier") +
         1.E-10);
    // Update deviatoric stress
    (*state_vars).at("q") =
        q_trial / (1 +
                   6 * (*state_vars).at("shear_modulus") *
                       (*state_vars).at("multiplier") / (m_value_ * m_value_) +
                   1.E-10);

    // Update yield function
    yield_type = this->compute_yield_state(&f_function, state_vars);
    // Counter iteration step
    ++counter_f;
  }

  // Set plastic tensor
  this->compute_plastic_tensor(stress, state_vars);
  // Check yield function
  (*state_vars).at("f_function") = f_function;
  // Compte current stress
  Vector6d updated_stress = trial_stress - this->dp_ * dstrain;
  // Incremental Volumetric strain
  double dvstrain = dstrain(0) + dstrain(1) + dstrain(2);
  // Update porosity
  (*state_vars).at("porosity") += dvstrain * (1 + e0_);

  return updated_stress;
}
