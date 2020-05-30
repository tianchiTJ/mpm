#ifndef MPM_MATERIAL_CONCRETE_H_
#define MPM_MATERIAL_CONCRETE_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! Concrete class
//! \brief Concrete material model
//! \details Concrete class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Concrete : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  Concrete(unsigned id, const Json& material_properties);

  //! Destructor
  ~Concrete() override{};

  //! Delete copy constructor
  Concrete(const Concrete&) = delete;

  //! Delete assignement operator
  Concrete& operator=(const Concrete&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! State variables
  std::vector<std::string> state_variables() const override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Compute elastic tensor
  bool compute_elastic_tensor();

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Tension strength
  double tension_strength_{std::numeric_limits<double>::max()};
  //! Compression strength
  double compression_strength_{std::numeric_limits<double>::max()};
};  // Concrete class
}  // namespace mpm

#include "concrete.tcc"

#endif  // MPM_MATERIAL_CONCRETE_H_
