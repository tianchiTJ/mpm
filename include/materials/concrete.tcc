//! Read material properties
template <unsigned Tdim>
mpm::Concrete<Tdim>::Concrete(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    tension_strength_ =
        material_properties.at("tension_strength").template get<double>();
    compression_strength_ =
        material_properties.at("compression_strength").template get<double>();
    // Calculate bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

    properties_ = material_properties;
    // Set elastic tensor
    this->compute_elastic_tensor();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::Concrete<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// Broken property
                               {"broken", 0.}};
  return state_vars;
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::Concrete<Tdim>::state_variables() const {
  const std::vector<std::string> state_vars = {"broken"};
  return state_vars;
}

//! Return elastic tensor
template <unsigned Tdim>
bool mpm::Concrete<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  // clang-format off
  // compute elasticityTensor
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on
  return true;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Concrete<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  // Initialize update stress
  Eigen::Matrix<double, 6, 1> updated_stress;
  updated_stress.setZero();

  if ((*state_vars).at("broken") != std::numeric_limits<double>::max()) {
    // Compute update stress
    updated_stress = stress + this->de_ * dstrain;
    // Check tension failure
    if (updated_stress.leftCols(Tdim).maxCoeff() > tension_strength_ ||
        updated_stress.leftCols(Tdim).maxCoeff() < compression_strength_) {
      // Set stress zero
      updated_stress.setZero();
      // Broken property
      (*state_vars).at("broken") = std::numeric_limits<double>::max();
    }
  }

  return updated_stress;
}
