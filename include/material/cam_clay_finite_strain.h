#ifndef MPM_MATERIAL_CAM_CLAY_FINITE_STRAIN_H_
#define MPM_MATERIAL_CAM_CLAY_FINITE_STRAIN_H_

#include <limits>

#include <cmath>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! CamClay class (With finite strain)
//! Base on Borja, Ronaldo I., and Claudio Tamagnini. "Cam-Clay plasticity Part
//! III: Extension of the infinitesimal model to include finite strains."
//! Computer Methods in Applied Mechanics and Engineering 155.1-2 (1998): 73-95.
//! Authors: Tianchi Zhao
//! \tparam Tdim Dimension
template <unsigned Tdim>
class CamClayFiniteStrain : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector3d = Eigen::Matrix<double, 3, 1>;
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 2 x 2
  using Matrix2x2 = Eigen::Matrix<double, 2, 2>;
  //! Define a Matrix of 3 x 3
  using Matrix3x3 = Eigen::Matrix<double, 3, 3>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Failure state
  enum FailureState { Elastic = 0, Yield = 1 };

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  CamClayFiniteStrain(unsigned id, const Json& material_properties);

  //! Destructor
  ~CamClayFiniteStrain() override{};

  //! Delete copy constructor
  CamClayFiniteStrain(const CamClayFiniteStrain&) = delete;

  //! Delete assignement operator
  CamClayFiniteStrain& operator=(const CamClayFiniteStrain&) = delete;

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

  //! Compute elastic modulus
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of elastic modulus
  bool compute_elastic_modulus(mpm::dense_map* state_vars);

  //! Compute elastic tensor
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of elastic tensor
  bool compute_elastic_tensor(mpm::dense_map* state_vars);

  //! Compute plastic tensor
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of plastic tensor
  bool compute_plastic_tensor(const Vector6d& stress,
                              mpm::dense_map* state_vars);

  //! Compute strain invariants
  //! \param[in] dstress Stress
  //! \param[in] dstrain Strain
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of strain invariants
  bool compute_strain_invariants(const Vector6d& dstrain,
                                 mpm::dense_map* state_vars);

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

  //! Compute dF/dmul
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] df_dp dF / dp
  //! \param[in] df_dq dF / dq
  //! \param[in] df_dpc dF / dpc
  //! \param[in] df_dmul dF / ddphi
  void compute_df_ddphi(const mpm::dense_map* state_vars, double* df_dp,
                        double* df_dq, double* df_dpc, double* df_ddphi);

  //! Compute hessian Matrix
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] df_dp dF / dp
  //! \param[in] df_dq dF / dq
  //! \param[in] df_dpc dF / dpc
  //! \param[in] operator_a A matrix
  void compute_hessian_matrix(const mpm::dense_map* state_vars,
                              const double* df_dp, const double* df_dq,
                              const double* df_dpc, Matrix3x3* operator_a);

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Plastic stiffness matrix
  Matrix6x6 dp_;
  // General parameters
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  // Constant part of shear modulus
  double mu0_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Porosity
  double e0_{std::numeric_limits<double>::max()};
  //! Volumetic stress
  double p0_{std::numeric_limits<double>::max()};
  //! Elastic volumetic strain
  double evstrain0_{std::numeric_limits<double>::max()};
  //! Plastic volumetic strain
  double pvstrain0_{std::numeric_limits<double>::max()};
  // Cam Clay parameters
  //! M
  double m_value_{std::numeric_limits<double>::max()};
  //! Initial preconsolidation pressure
  double pc0_{std::numeric_limits<double>::max()};
  //! Alpha
  double alpha_{std::numeric_limits<double>::max()};
  //! Lambda
  double lambda_{std::numeric_limits<double>::max()};
  //! Kappa
  double kappa_{std::numeric_limits<double>::max()};
  //! Lambda in finite deformation
  double lambda_fs_{std::numeric_limits<double>::max()};
  //! Kappa in finite deformation
  double kappa_fs_{std::numeric_limits<double>::max()};

};  // CamClayFiniteStrain class
}  // namespace mpm

#include "cam_clay_finite_strain.tcc"

#endif  // MPM_MATERIAL_CAM_CLAY_Finite_Srain_H_
