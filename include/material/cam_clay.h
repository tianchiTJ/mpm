#ifndef MPM_MATERIAL_CAM_CLAY_H_
#define MPM_MATERIAL_CAM_CLAY_H_

#include <limits>

#include <cmath>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! CamClay class
//! \brief Cam Clay material model
//! \tparam Tdim Dimension
template <unsigned Tdim>
class CamClay : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Failure state
  enum FailureState { Elastic = 0, Yield = 1 };

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  CamClay(unsigned id, const Json& material_properties);

  //! Destructor
  ~CamClay() override{};

  //! Delete copy constructor
  CamClay(const CamClay&) = delete;

  //! Delete assignement operator
  CamClay& operator=(const CamClay&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! Thermodynamic pressure
  //! \param[in] volumetric_strain dVolumetric_strain
  //! \retval pressure Pressure for volumetric strain
  double thermodynamic_pressure(double volumetric_strain) const override {
    return 0;
  };

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

  //! Compute stress invariants (j2, j3, rho, theta, and epsilon)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of stress invariants
  bool compute_stress_invariants(const Vector6d& stress,
                                 mpm::dense_map* state_vars);

  //! Compute yield function and yield state
  //! \param[in] state_vars History-dependent state variables
  //! \retval yield_type Yield type (elastic, shear or tensile)
  FailureState compute_yield_state(double* yield_function,
                                   const mpm::dense_map* state_vars);

  //! Compute dF/dSigma and dY/dSigma
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] stress Stress
  //! \param[in] df_dsigma dF/dSigma
  //! \param[in] dp_dsigma dY/dSigma
  void compute_df_dy(const mpm::dense_map* state_vars, const Vector6d& stress,
                     Vector6d* df_dsigma, Vector6d* dy_dsigma);

  //! Compute dF/dmul
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] df_dmul dF / dmultiplier
  void compute_df_dmul(const mpm::dense_map* state_vars, double* df_dmul);

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Compute elastic tensor
  bool compute_elastic_tensor(mpm::dense_map* state_vars);
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  // General parameters
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Initial porosity
  double e0_{std::numeric_limits<double>::max()};
  // Cam Clay parameters
  //! M
  double m_value_{std::numeric_limits<double>::max()};
  //! Initial preconsolidation pressure
  double p0_{std::numeric_limits<double>::max()};
  //! Lambda
  double lambda_{std::numeric_limits<double>::max()};
  //! Kappa
  double kappa_{std::numeric_limits<double>::max()};

};  // CamClay class
}  // namespace mpm

#include "cam_clay.tcc"

#endif  // MPM_MATERIAL_CAM_CLAY_H_
