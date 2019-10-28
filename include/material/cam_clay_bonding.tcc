//! Critical state ratio (M(thera)) is calculated based on Jefferies and Shuttle
//! (2011) Constructor with id and material properties
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
    // M (or Mtc in "Three invariants type")
    m_ = material_properties["m"].template get<double>();
    // Lambda
    lambda_ = material_properties["lambda"].template get<double>();
    // Kappa
    kappa_ = material_properties["kappa"].template get<double>();
    // Initial void ratio
    e0_ = e_ref_ - lambda_ * log(pc0_ / p_ref_) - kappa_ * log(ocr_);
    // Three invariants
    three_invariants_ =
        material_properties["three_invariants"].template get<bool>();
    // Bonding properties
    // Bonding status
    bonding_ = material_properties["bonding"].template get<bool>();
    // Bonding material constants
    if (bonding_) {
      // Hydrate saturation
      s_h_ = material_properties["s_h"].template get<bool>();
      // Material constants a
      mc_a_ = material_properties["mc_a"].template get<bool>();
      // Material constants a
      mc_b_ = material_properties["mc_b"].template get<bool>();
      // Material constants a
      mc_c_ = material_properties["mc_c"].template get<bool>();
      // Material constants a
      mc_d_ = material_properties["mc_d"].template get<bool>();
      // Degradation
      degradation_ = material_properties["degradation"].template get<bool>();
    }
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
      // Cam clay parameters
      // Preconsolidation pressure
      {"pc", pc0_},
      // void_ratio
      {"void_ratio", e0_},
      // Consistency multiplier
      {"multiplier", 0.},
      // M_theta
      {"m_theta", m_},
      // Yield function
      {"f_function", 0.},
      // Incremental plastic strain
      // Incremental plastic volumetic strain
      {"dpvstrain", 0.},
      // Incremental plastic deviatoric strain
      {"dpdstrain", 0.},
      // Bonding parameters
      // Chi
      {"chi", 1.},
      // Pcd
      {"pcd", 0.},
      // Pcc
      {"pcc", 0.}};
  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::CamClay<Tdim>::compute_elastic_tensor(mpm::dense_map* state_vars) {
  // Compute bulk modulus
  if ((*state_vars).at("p") > 0)
    (*state_vars).at("bulk_modulus") =
        (1 + (*state_vars).at("void_ratio")) / kappa_ * (*state_vars).at("p");
  // Compute Shear modulus
  (*state_vars).at("shear_modulus") = 3 * (*state_vars).at("bulk_modulus") *
                                      (1 - 2 * poisson_ratio_) /
                                      (2 * (1 + poisson_ratio_));
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
  // Bonding parameters
  const double pcc = (*state_vars).at("pcc");
  const double pcd = (*state_vars).at("pcd");
  // Compute dF / dp
  const double df_dp = 2 * p - pc - pcd;
  // Compute dF / dq
  const double df_dq = 2 * q / (m_ * m_);
  // Upsilon
  const double upsilon =
      (1 + (*state_vars).at("void_ratio")) / (lambda_ - kappa_);
  // Coefficients in plastic stiffness matrix
  const double a1 = pow(((*state_vars).at("bulk_modulus") * df_dp), 2);
  const double a2 = -sqrt(6) * (*state_vars).at("bulk_modulus") * df_dp *
                    (*state_vars).at("shear_modulus") * df_dq;
  const double a3 = 6 * pow(((*state_vars).at("shear_modulus") * df_dq), 2);
  // hardening parameter
  double hardening = (*state_vars).at("bulk_modulus") * (df_dp * df_dp) +
                     3 * (*state_vars).at("shear_modulus") * (df_dq * df_dq) +
                     upsilon * p * pc * df_dp;
  // Compute bonding
  if (bonding_) {
    // Compute ds_h_mec / dpdstrain
    const double ds_dpds = -s_h_ * degradation_ * (*state_vars).at("chi");
    // Compute pcd hardening parameter
    const double hardening_pcd =
        mc_a_ * mc_b_ * pow((*state_vars).at("chi") * s_h_, mc_b_ - 1) *
        ds_dpds;
    // Compute pcc hardening parameter
    const double hardening_pcc =
        mc_c_ * mc_d_ * pow((*state_vars).at("chi") * s_h_, mc_d_ - 1) *
        ds_dpds;
    // Update hardening parameter
    hardening += (hardening_pcd + hardening_pcc);
  }
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
  this->dp_ = 1. / hardening * (a1 * l_l + a2 * (n_l + l_n) + a3 * (n_n));

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
  (*state_vars).at("j2") =
      (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
       pow((stress(0) - stress(2)), 2)) /
          6.0 +
      pow(stress(3), 2);
  if (Tdim == 3)
    (*state_vars).at("j2") += pow(stress(4), 2) + pow(stress(5), 2);
  // Compute q
  (*state_vars)["q"] = sqrt(3 * (*state_vars).at("j2"));
  // Compute J3
  (*state_vars).at("j3") = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
                           (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    (*state_vars).at("j3") +=
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
  (*state_vars).at("theta") = (1. / 3.) * acos(theta_val);
  // Check theta
  if ((*state_vars).at("theta") > M_PI / 3.)
    (*state_vars).at("theta") = M_PI / 3.;
  if ((*state_vars).at("theta") < 0.0) (*state_vars).at("theta") = 0.;

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
  const double m_theta = (*state_vars).at("m_theta");
  // Plastic volumetic strain
  const double pc = (*state_vars).at("pc");
  // Get bonding parameters
  const double pcd = (*state_vars).at("pcd");
  const double pcc = (*state_vars).at("pcc");
  // Initialise yield status (0: elastic, 1: yield)
  auto yield_type = FailureState::Elastic;
  // Compute yield functions
  (*yield_function) =
      q * q / (m_theta * m_theta) + (p + pcc) * (p - pc - pcd - pcc);
  // Tension failure
  if ((*yield_function) > 1.E-22) yield_type = FailureState::Yield;

  return yield_type;
}

//! Compute bonding parameters
template <unsigned Tdim>
void mpm::CamClay<Tdim>::compute_bonding_parameters(
    const double* chi, mpm::dense_map* state_vars) {
  // Compute chi
  (*state_vars).at("chi") =
      (*chi) + (-degradation_ * (*state_vars).at("chi") *
                (*state_vars).at("chi") * (*state_vars).at("dpdstrain"));
  // Compute mechanical hydrate saturation
  const double s_h_mec = (*state_vars).at("chi") * s_h_;
  // Compute pcd
  (*state_vars).at("pcd") = mc_a_ * pow(s_h_mec, mc_b_);
  // Compute pcc
  (*state_vars).at("pcc") = mc_c_ * pow(s_h_mec, mc_d_);
}

//! Compute dF/dmul
template <unsigned Tdim>
void mpm::CamClay<Tdim>::compute_df_dmul(const mpm::dense_map* state_vars,
                                         double* df_dmul) {
  // Stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  const double m_theta = (*state_vars).at("m_theta");
  // Preconsolidation pressure
  const double pc = (*state_vars).at("pc");
  // Get bonding parameters
  const double pcd = (*state_vars).at("pcd");
  const double pcc = (*state_vars).at("pcc");
  // Get elastic modulus
  const double e_b = (*state_vars).at("bulk_modulus");
  const double e_s = (*state_vars).at("shear_modulus");
  // Get consistency parameter
  const double mul = (*state_vars).at("multiplier");
  // Compute dF / dp
  double df_dp = 2 * p - pc - pcd;
  // Compute dF / dq
  double df_dq = 2 * q / (m_theta * m_theta);
  // Compute dF / dpc
  double df_dpc = -(p + pcc);
  // Upsilon
  double upsilon = (1 + (*state_vars).at("void_ratio")) / (lambda_ - kappa_);
  // A_den
  double a_den = 1 + (2 * e_b + upsilon * pc) * mul;
  // Compute dp / dmul
  double dp_dmul = -e_b * (2 * p - pc - pcd) / a_den;
  // Compute dpc / dmul
  double dpc_dmul = upsilon * pc * (2 * p - pc - pcd) / a_den;
  // Compute dq / dmul
  double dq_dmul = -q / (mul + m_theta * m_theta / (6 * e_s));
  // Compute dF / dmul
  if (!bonding_)
    (*df_dmul) = (df_dp * dp_dmul) + (df_dq * dq_dmul) + (df_dpc * dpc_dmul);
  // Compute bonding parameters
  if (bonding_) {
    // Compute dF / dpcd
    double df_dpcd = -p - pcc;
    // Compute dF / dpcc
    double df_dpcc = -2 * pcc - pc - pcd;
    // Compute dpcd / dmul
    double dpcd_dmul =
        -2 * mc_b_ * degradation_ * q / (6 * e_s * mul + pow(m_, 2)) * pcd;
    // Compute dpcc / dmul
    double dpcc_dmul =
        -2 * mc_d_ * degradation_ * q / (6 * e_s * mul + pow(m_, 2)) * pcc;
    // Compute dp / dmul
    dp_dmul +=
        (e_b * mul *
         (pc * upsilon * mul * (1 + pc * upsilon * mul - 2 * e_b * mul) + 1) /
         ((1 + 2 * e_b * mul) * (1 + pc * upsilon * mul) *
          (1 + 2 * e_b * mul + pc * upsilon * mul)) *
         dpcd_dmul);
    // Compute dpc / dmul
    dpc_dmul +=
        (pc * upsilon * mul *
         (2 * e_b * mul * (pc * upsilon * mul - 2 * e_b * mul - 1) - 1) /
         (pow((1 + 2 * e_b * mul), 2) * pow((1 + pc * upsilon * mul), 2)) *
         dpcd_dmul);
    // Compute dF / dmul
    (*df_dmul) = (df_dp * dp_dmul) + (df_dq * dq_dmul) + (df_dpc * dpc_dmul) +
                 (df_dpcd * dpcd_dmul) + (df_dpcc * dpcc_dmul);
  }
}

//! Compute dF/dSigma
template <unsigned Tdim>
void mpm::CamClay<Tdim>::compute_df_dsigma(const mpm::dense_map* state_vars,
                                           const Vector6d& stress,
                                           Vector6d* dtheta_dsigma,
                                           Vector6d* df_dsigma) {
  // Stress invariants
  const double p = (*state_vars).at("p");
  const double q = (*state_vars).at("q");
  const double j2 = (*state_vars).at("j2");
  const double j3 = (*state_vars).at("j3");
  const double m_theta = (*state_vars).at("m_theta");
  // Preconsolidation pressure
  const double pc = (*state_vars).at("pc");
  // Bonding parameters
  const double pcd = (*state_vars).at("pcd");
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

  // Compute dF / dp
  double df_dp = 2 * p - pc - pcd;
  // Compute dF / dq
  double df_dq = 2 * q / (m_theta * m_theta);
  // Initialise dF / dMtheta
  double df_dm = -2 * q * q / pow(m_theta, 3);
  // Initialise dMtheta / dtheta
  double dm_dtheta = -1.5 * pow((*state_vars).at("m_theta"), 3) *
                     sin(1.5 * (*state_vars).at("theta")) /
                     pow(3 + (*state_vars).at("m_theta") *
                                 cos(1.5 * (*state_vars).at("theta")),
                         2);
  // Compute dp / dSigma
  Vector6d dp_dsigma = Vector6d::Zero();
  dp_dsigma(0) = -1. / 3.;
  dp_dsigma(1) = -1. / 3.;
  dp_dsigma(2) = -1. / 3.;
  // Compute dq / dSigma
  Vector6d dq_dsigma = Vector6d::Zero();
  double multiplier = 1.;
  if (fabs(q) > 0.) multiplier = 1.5 / q;
  dq_dsigma = multiplier * dev_stress;
  if (Tdim == 2) {
    dq_dsigma(4) = 0.;
    dq_dsigma(5) = 0.;
  }
  // Compute r
  double r_val = 0.;
  if (fabs(j2) > 1.E-22) r_val = (3. * sqrt(3.) / 2.) * (j3 / pow(j2, 1.5));
  // Compute dTheta / dr
  double divider = 1 - (r_val * r_val);
  if (divider <= 0.) divider = 1.E-3;
  double dtheta_dr = -1 / (3. * sqrt(divider));
  // Compute dr / dJ2
  double dr_dj2 = (-9 * sqrt(3.) / 4.) * j3;
  if (fabs(j2) > 1.E-22) dr_dj2 = dr_dj2 / pow(j2, 2.5);
  // Compute dr / dJ3
  double dr_dj3 = 1.5 * sqrt(3.);
  if (fabs(j2) > 1.E-22) dr_dj3 = dr_dj3 / pow(j2, 1.5);
  // Compute dJ2 / dSigma
  Vector6d dj2_dsigma = dev_stress;
  // Compute dJ3 / dSigma
  Eigen::Matrix<double, 3, 1> dev1;
  dev1(0) = dev_stress(0);
  dev1(1) = dev_stress(3);
  dev1(2) = dev_stress(5);
  Eigen::Matrix<double, 3, 1> dev2;
  dev2(0) = dev_stress(3);
  dev2(1) = dev_stress(1);
  dev2(2) = dev_stress(4);
  Eigen::Matrix<double, 3, 1> dev3;
  dev3(0) = dev_stress(5);
  dev3(1) = dev_stress(4);
  dev3(2) = dev_stress(2);
  Vector6d dj3_dsigma = Vector6d::Zero();
  dj3_dsigma(0) = dev1.dot(dev1) - (2. / 3.) * j2;
  dj3_dsigma(1) = dev2.dot(dev2) - (2. / 3.) * j2;
  dj3_dsigma(2) = dev3.dot(dev3) - (2. / 3.) * j2;
  dj3_dsigma(3) = dev1.dot(dev2);
  if (Tdim == 3) {
    dj3_dsigma(4) = dev2.dot(dev3);
    dj3_dsigma(5) = dev1.dot(dev3);
  }
  // Compute dtheta / dsigma
  (*dtheta_dsigma) =
      dtheta_dr * ((dr_dj2 * dj2_dsigma) + (dr_dj3 * dj3_dsigma));
  if (Tdim == 2) {
    (*dtheta_dsigma)(4) = 0.;
    (*dtheta_dsigma)(5) = 0.;
  }
  (*df_dsigma) = (df_dp * dp_dsigma) + (df_dq * dq_dsigma) +
                 (df_dm * dm_dtheta * (*dtheta_dsigma));
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
  double e_index =
      upsilon * (*state_vars).at("multiplier") *
      (2 * p_trial - (*state_vars).at("pc") - (*state_vars).at("pcd")) /
      (1 +
       2 * (*state_vars).at("multiplier") * (*state_vars).at("bulk_modulus"));
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
  // Compute current mean pressure
  (*state_vars).at("p") = -(stress(0) + stress(1) + stress(2)) / 3.;
  // Set elastic tensor
  this->compute_elastic_tensor(state_vars);
  //-------------------------------------------------------------------------
  // Elastic step
  // Compute trial stress
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  // Compute trial stress invariants
  this->compute_stress_invariants(trial_stress, state_vars);
  // Compute bonding parameters
  const double chi_n = (*state_vars).at("chi");
  if (bonding_) this->compute_bonding_parameters(&chi_n, state_vars);
  // Update Mtheta
  if (three_invariants_)
    (*state_vars).at("m_theta") =
        m_ - pow(m_, 2) / (3 + m_ * cos(1.5 * (*state_vars).at("theta")));
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
  // Lode angle of trial stress
  const double theta_trial = (*state_vars).at("theta");
  // Preconsolidation pressure of last step
  const double pc_n = (*state_vars).at("pc");
  // Initialise dF / dmul
  double df_dmul = 0;
  // Initialise dtheta / dsigma
  Vector6d dtheta_dsigma = Vector6d::Zero();
  // Initialise dF / dsigma
  Vector6d df_dsigma = Vector6d::Zero();
  // Initialise dtheta
  double dtheta = 0.;
  if (three_invariants_) {
    // Compute dF / dsigma for three invariants
    this->compute_df_dsigma(state_vars, trial_stress, &dtheta_dsigma,
                            &df_dsigma);
    // Compute dtheta
    dtheta = dtheta_dsigma.transpose() * df_dsigma;
  }
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
                           (*state_vars).at("multiplier") / (m_ * m_));
    // Compute plastic volumetic strain
    (*state_vars).at("dpvstrain") =
        (*state_vars).at("multiplier") *
        (2 * (*state_vars).at("p") - (*state_vars).at("pc") -
         (*state_vars).at("pcd"));
    // Compute plastic deviatoric strain
    (*state_vars).at("dpdstrain") =
        (*state_vars).at("multiplier") *
        (2 * (*state_vars).at("q") / pow((*state_vars).at("m_theta"), 2));
    // Update bonding parameters
    if (bonding_) this->compute_bonding_parameters(&chi_n, state_vars);
    // Compute three invariants parameters
    if (three_invariants_) {
      // Update lode angle
      (*state_vars).at("theta") =
          theta_trial - (*state_vars).at("multiplier") * dtheta;
      // Update Mtheta
      (*state_vars).at("m_theta") =
          m_ - pow(m_, 2) / (3 + m_ * cos(1.5 * (*state_vars).at("theta")));
    }
    // Update yield function
    yield_type = this->compute_yield_state(&f_function, state_vars);
    // Counter iteration step
    ++counter_f;
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
