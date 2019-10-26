//! Constructor with id and material properties
template <unsigned Tdim>
mpm::CamClayFiniteStrain<Tdim>::CamClayFiniteStrain(
    unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // General parameters
    // Density
    density_ = material_properties["density"].template get<double>();
    // Young's modulus
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    // Constant part of shear modulus
    mu0_ = material_properties["mu0"].template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    // Initial status
    // Initial porosity
    e0_ = material_properties["e0"].template get<double>();
    // Initial volumetic stress
    p0_ = material_properties["p0"].template get<double>();
    // Initial elastic volumetic strain
    evstrain0_ = material_properties["evstrain0"].template get<double>();
    // Initial plastic volumetic strain
    pvstrain0_ = material_properties["pvstrain0"].template get<double>();
    // Properties
    properties_ = material_properties;
    // Cam Clay parameters
    // M
    m_value_ = material_properties["m_value"].template get<double>();
    // Initial preconsolidation pressure
    pc0_ = material_properties["pc0"].template get<double>();
    // Alpha
    alpha_ = material_properties["alpha"].template get<double>();
    // Lambda
    lambda_ = material_properties["lambda"].template get<double>();
    // Kappa
    kappa_ = material_properties["kappa"].template get<double>();
    // Lambda in finite deformation
    lambda_fs_ = lambda_ / (1. - lambda_);
    // Kappa in finite deformation
    kappa_fs_ = kappa_ / (1. - kappa_);

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::CamClayFiniteStrain<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// Elastic modulus
                               // Bulk modulus
                               {"bulk_modulus", youngs_modulus_},
                               // Shear modulus
                               {"shear_modulus", mu0_},
                               // Stress invariants
                               // Volumetric stress
                               {"p", 0.},
                               // Deviatoric stress
                               {"q", 0.},
                               // Lode's angle
                               // {"theta", 0.},
                               // Elastic strain
                               // Elastic strain components
                               {"estrain0", 0.},
                               {"estrain1", 0.},
                               {"estrain2", 0.},
                               {"estrain3", 0.},
                               {"estrain4", 0.},
                               {"estrain5", 0.},
                               // Elastic volumetric strain
                               {"evstrain", evstrain0_},
                               // Elastic deviatoric strain
                               {"esstrain", 0.},
                               // Cam clay parameters
                               // Omega
                               {"omega", 0.},
                               // consistency parameter
                               {"delta_phi", 0.},
                               // Preconsolidation pressure
                               {"pc", pc0_},
                               // Porosity
                               {"porosity", e0_}};

  return state_vars;
}

//! Compute elastic modulus
template <unsigned Tdim>
bool mpm::CamClayFiniteStrain<Tdim>::compute_elastic_modulus(
    mpm::dense_map* state_vars) {
  bool status = false;
  // Compute slastic modulus
  if ((*state_vars).at("p") > 0) {
    // Compute bulk modulus
    (*state_vars).at("bulk_modulus") =
        (1 + (*state_vars).at("porosity")) / kappa_fs_ * (*state_vars).at("p");
    // Compute shear modulus
    (*state_vars).at("shear_modulus") =
        mu0_ - alpha_ * p0_ * exp((*state_vars).at("omega"));
    // Computation status
    status = true;
  }
  return status;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::CamClayFiniteStrain<Tdim>::compute_elastic_tensor(
    mpm::dense_map* state_vars) {
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
bool mpm::CamClayFiniteStrain<Tdim>::compute_plastic_tensor(
    const Vector6d& stress, mpm::dense_map* state_vars) {
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

//! Compute strain invariants
template <unsigned Tdim>
bool mpm::CamClayFiniteStrain<Tdim>::compute_strain_invariants(
    const Vector6d& dstrain, mpm::dense_map* state_vars) {
  // Compute elastic volumetic strain
  (*state_vars).at("evstrain") -= (dstrain(0) + dstrain(1) + dstrain(2));
  // Update elastic strain
  (*state_vars).at("estrain0") += dstrain(0);
  (*state_vars).at("estrain1") += dstrain(1);
  (*state_vars).at("estrain2") += dstrain(2);
  (*state_vars).at("estrain3") += dstrain(3) / 2.;
  if (Tdim == 3) {
    (*state_vars).at("estrain4") += strain(4) / 2.;
    (*state_vars).at("estrain5") += strain(5) / 2.;
  }
  // Compute the elastic deviatoric strain
  Vector6d edev_strain = Vector6d::Zero();
  edev_strain(0) = (*state_vars).at("estrain0") + (*state_vars).at("evstrain");
  edev_strain(1) = (*state_vars).at("estrain1") + (*state_vars).at("evstrain");
  edev_strain(2) = (*state_vars).at("estrain2") + (*state_vars).at("evstrain");
  edev_strain(3) = (*state_vars).at("estrain3");
  if (Tdim == 3) {
    edev_strain(4) = (*state_vars).at("estrain4");
    edev_strain(5) = (*state_vars).at("estrain5");
  }
  // Compute I2
  double i2 = (pow((edev_strain(0) - edev_strain(1)), 2) +
               pow((edev_strain(1) - edev_strain(2)), 2) +
               pow((edev_strain(0) - edev_strain(2)), 2)) /
                  6.0 +
              pow(edev_strain(3), 2);
  if (Tdim == 3) i2 += pow(edev_strain(4), 2) + pow(edev_strain(5), 2);
  // Compute deviatoric strain
  (*state_vars).at("esstrain") = sqrt(4. / 3. * i2);

  return true;
}

//! Compute stress invariants
template <unsigned Tdim>
bool mpm::CamClayFiniteStrain<Tdim>::compute_stress_invariants(
    const Vector6d& stress, mpm::dense_map* state_vars) {
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
  (*state_vars).at("q") = sqrt(3 * j2);

  return true;
}

//! Compute yield function and yield state
template <unsigned Tdim>
typename mpm::CamClayFiniteStrain<Tdim>::FailureState
    mpm::CamClayFiniteStrain<Tdim>::compute_yield_state(
        double* yield_function, const mpm::dense_map* state_vars) {
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
void mpm::CamClayFiniteStrain<Tdim>::compute_df_ddphi(
    const mpm::dense_map* state_vars, double* df_dp, double* df_dq,
    double* df_dpc, double* df_ddphi) {
  // Stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  // Preconsolidation pressure
  const double pc = (*state_vars).at("pc");
  // Compute dF / dp
  (*df_dp) = 2 * p - pc;
  // Compute dF / dq
  (*df_dq) = 2 * q / (m_value_ * m_value_);
  // Compute dF / dpc
  (*df_dpc) = -p;
  // Upsilon
  double upsilon = (1 + (*state_vars).at("porosity")) / (lambda_ - kappa_);
  // A_den
  double a_den = 1 + (2 * (*state_vars).at("bulk_modulus") + upsilon * pc) *
                         (*state_vars).at("delta_phi");
  // Compute dp / dmultiplier
  double dp_dmul =
      -(*state_vars).at("bulk_modulus") * (2 * p - pc) / (a_den + 1.E-10);
  // Compute dpc / dmultiplier
  double dpc_dmul = upsilon * pc * (2 * p - pc) / (a_den + 1.E-10);
  // Compute dq / dmultiplier
  double dq_dmul =
      -q /
      ((*state_vars).at("delta_phi") +
       m_value_ * m_value_ / (6 * (*state_vars).at("shear_modulus")) + 1.E-10);
  // Compute dF / ddphitiplier
  (*df_ddphi) =
      ((*df_dp) * dp_dmul) + ((*df_dq) * dq_dmul) + ((*df_dpc) * dpc_dmul);
}

//! Compute hessian Matrix
template <unsigned Tdim>
void mpm::CamClayFiniteStrain<Tdim>::compute_hessian_matrix(
    const mpm::dense_map* state_vars, const double* df_dp, const double* df_dq,
    const double* df_dpc, Matrix3x3* operator_a) {
  // Stress invariants
  const double p = (*state_vars).at("p");
  // Compute plastic hardening modulus
  const double k_p = (*state_vars).at("pc") / (lambda_fs_ - kappa_fs_);
  // Initialise Hessian matrix
  Matrix3x3 hessian_psi, hessian_psi_f;
  // Compute Hessian matrix of Psi
  hessian_psi(0, 0) = -p / kappa_fs_;
  hessian_psi(1, 1) = 3 * (*state_vars).at("shear_modulus");
  hessian_psi(0, 1) = hessian_psi(1, 0) =
      3 * p0_ * alpha_ * (*state_vars).at("esstrain") / kappa_fs_ *
      exp((*state_vars).at("omega"));
  /// Compute Hessian matrix of yield function
  hessian_psi_f(0, 0) = 2 * hessian_psi(0, 0);
  hessian_psi_f(1, 1) = 2 * hessian_psi(1, 1) / (m_value_ * m_value_);
  hessian_psi_f(0, 1) = 2 * hessian_psi(0, 1);
  hessian_psi_f(1, 0) = 2 * hessian_psi(1, 0) / (m_value_ * m_value_);
  // Compute A matrix
  (*operator_a)(0, 0) =
      1 + (*state_vars).at("delta_phi") * (hessian_psi_f(0, 0) - k_p);
  (*operator_a)(0, 1) = (*state_vars).at("delta_phi") * hessian_psi_f(0, 1);
  (*operator_a)(0, 2) = (*df_dp);
  (*operator_a)(1, 0) = (*state_vars).at("delta_phi") * hessian_psi_f(1, 0);
  (*operator_a)(1, 1) = 1 + (*state_vars).at("delta_phi") * hessian_psi_f(1, 1);
  (*operator_a)(1, 2) = (*df_dq);
  (*operator_a)(2, 0) = hessian_psi(0, 0) * (*df_dp) +
                        hessian_psi(1, 0) * (*df_dq) + k_p * (*df_dpc);
  (*operator_a)(2, 1) =
      hessian_psi(0, 1) * (*df_dp) + hessian_psi(1, 1) * (*df_dq);
  (*operator_a)(2, 2) = 0.;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::CamClayFiniteStrain<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  // Phase
  const unsigned phase = 0;
  // Maximum iteration step number
  const int itrstep = 100;
  // Maximum subiteration step number
  const int substep = 100;
  // Tolerance for yield function
  const double Ftolerance = 1.E-1;
  // Get current strain vector
  const Vector6d strain = (*ptr)->strain(phase);
  // Compute Jacobian
  const double jacobian = exp(-strain(0) - strain(1) - strain(2));
  // Update volumetric stress
  (*state_vars).at("p") = -jacobian * (stress(0) + stress(1) + stress(2)) / 3.;
  // Compute elastic modulus
  this->compute_elastic_modulus(state_vars);
  // Set elastic tensor
  this->compute_elastic_tensor(state_vars);
  //-------------------------------------------------------------------------
  // Elastic step
  // Compute incremental stress
  const Vector6d dstress = this->de_ * dstrain / jacobian;
  // Compute trial stress
  const Vector6d trial_stress = stress + dstress;
  // Compute trial stress invariants
  this->compute_stress_invariants(trial_stress * jacobian, state_vars);
  // Compute trial strain invariants
  this->compute_strain_invariants(dstrain, state_vars);
  // Initialise value for yield function
  double f_function;
  // Check yield
  auto yield_type = this->compute_yield_state(&f_function, state_vars);
  // Return the updated stress in elastic state
  if (yield_type == FailureState::Elastic) return trial_stress;
  //-------------------------------------------------------------------------
  // Plastic step
  // Counters for interations
  int counter_f = 0;
  int counter_g = 0;
  // Initialise consistency parameter
  (*state_vars).at("delta_phi") = 0.;
  // Volumetric trial strain
  const double evstrain_trial = (*state_vars).at("evstrain");
  // Deviatoric trial strain
  const double esstrain_trial = (*state_vars).at("esstrain");
  // Volumetric trial stress
  const double p_trial = (*state_vars).at("p");
  // Deviatoric trial stress
  const double q_trial = (*state_vars).at("q");
  // Preconsolidation pressure of last step
  const double pc_n = (*state_vars).at("pc");
  // Initialise dF / dp, dF / dq, dF / dpc, dF / ddphi
  double df_dp, df_dq, df_dpc, df_ddphi;
  // Initialise residual vector
  // 0 for evstrain, 1 for esstrain, 2 for yield function
  Vector3d r_vector = Vector3d::Zero();
  // Iteration for consistency parameter
  while (fabs(f_function) > Ftolerance && counter_f < itrstep) {
    // Compute dF / ddphi
    this->compute_df_ddphi(state_vars, &df_dp, &df_dq, &df_dpc, &df_ddphi);
    // Compute residual vector
    r_vector(0) = (*state_vars).at("evstrain") - evstrain_trial +
                  (*state_vars).at("delta_phi") * df_dp;
    r_vector(1) = (*state_vars).at("esstrain") - esstrain_trial +
                  (*state_vars).at("delta_phi") * df_dq;
    r_vector(2) = f_function;
    // Initialise A Matrix
    Matrix3x3 operator_a;
    // Compute Hessian matrix and A Matrix
    this->compute_hessian_matrix(state_vars, &df_dp, &df_dq, &df_dpc,
                                 &operator_a);
    // Compute unkown vector
    Vector3d delta_x = -operator_a.inverse() * r_vector;
    // Update elastic volumetric strain
    (*state_vars).at("evstrain") += delta_x(0);
    // Update elastic deviatoric strain
    (*state_vars).at("esstrain") += delta_x(1);
    // Update consistency parameter
    (*state_vars).at("delta_phi") += delta_x(2);
    // Update pc
    (*state_vars).at("pc") =
        pc_n * exp(-1. / (lambda_fs_ - kappa_fs_) *
                   (evstrain_trial - (*state_vars).at("evstrain")));
    // Compute omega
    (*state_vars).at("omega") =
        -((*state_vars).at("evstrain") - evstrain0_) / kappa_fs_;
    // Update elastic modulus
    this->compute_elastic_modulus(state_vars);
    // Update volumetric stress
    (*state_vars).at("p") = p0_ * exp((*state_vars).at("omega")) *
                            (1 + 3 * alpha_ / (2 * kappa_fs_) *
                                     pow((*state_vars).at("esstrain"), 2));
    // Update deviatoric stress
    (*state_vars).at("q") =
        3 * (*state_vars).at("shear_modulus") * (*state_vars).at("esstrain");
    // Update yield function
    yield_type = this->compute_yield_state(&f_function, state_vars);
    // Counter iteration step
    ++counter_f;
  }
  // Set plastic tensor
  this->compute_plastic_tensor(stress * jacobian, state_vars);
  // Compute current stress
  Vector6d updated_stress = trial_stress - this->dp_ * dstrain / jacobian;
  // Update porosity
  (*state_vars).at("porosity") =
      (1 + e0_) * pow((pc0_ / (*state_vars).at("pc")), lambda_) - 1.;

  return updated_stress;
}
