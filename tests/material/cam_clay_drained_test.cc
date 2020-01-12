#include <fstream>
#include <limits>
#include <vector>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/cam_clay.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

using Json = nlohmann::json;

//! Check drained condition with bonding for hardening in 3D
TEST_CASE("drained condition hardening is checked in 3D",
          "[material][cam_clay][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;
  // Dimension
  const unsigned Dim = 3;
  /*
  // Incremental axial stress
  const double dstrain_axial = -1.E-5;
  // Number of steps
  const unsigned long long nsteps = 200000;
  // Interval of output step;
  const unsigned output_steps = 100;
  // Confining pressure
  const double confining_p = 8000 * 6894.75729;
  // Initial void ratio
  const double e0 = 0.38;
  // Lambda
  const double lambda = 0.2;
  // Kappa
  const double kappa = 0.008;
  // OCR
  const double ocr = 1;
  // pc0
  const double p_ref = 1000000;
  // Reference void ratio
  const double e_ref =
      e0 + lambda * log(confining_p / p_ref) + kappa * log(ocr);
  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);
  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 2000.;
  jmaterial["youngs_modulus"] = 1.0E+10;
  jmaterial["poisson_ratio"] = 0.17;
  jmaterial["p_ref"] = p_ref;
  jmaterial["pc0"] = confining_p;
  jmaterial["ocr"] = ocr;
  jmaterial["m"] = 1.2;
  jmaterial["lambda"] = lambda;
  jmaterial["kappa"] = kappa;
  jmaterial["three_invariants"] = false;
  jmaterial["bonding"] = true;
  jmaterial["e_ref"] = e_ref;
  // Bonding parameters
  jmaterial["s_h"] = 1.;
  jmaterial["mc_a"] = 3.E6;
  jmaterial["mc_b"] = 1.;
  jmaterial["mc_c"] = 3.E6;
  jmaterial["mc_d"] = 1.;
  jmaterial["m_degradation"] = 1.;
  */
  // Input file
  std::ifstream input_file;
  input_file.open("./Exxon_data/BMCC-parameters.json");
  // JSON object
  Json json_ = Json::parse(input_file);
  // Incremental axial stress
  const double dstrain_axial = json_["dstrain_axial"];
  // Number of steps
  const unsigned long long nsteps = json_["nsteps"];
  // Interval of output step;
  const unsigned output_steps = json_["output_steps"];
  // Confining pressure
  double confining_p0 = json_["confining_p"];
  confining_p0 *= 6894.75729;
  // Confining pressure multipliers
  const std::vector<double> multipliers = json_["multipliers"];
  // Initial void ratio
  const double e0 = json_["e0"];
  // Lambda
  const double lambda = json_["lambda"];
  // Kappa
  const double kappa = json_["kappa"];
  // pc0
  const double pc0 = json_["pc0"];
  // Iterate over multipliers
  for (int j = 0; j < multipliers.size(); ++j) {
    // Confining pressure
    const double confining_p = confining_p0 * multipliers[j];
    // Add particle
    mpm::Index pid = 0;
    Eigen::Matrix<double, Dim, 1> coords;
    coords.setZero();
    auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = json_["density"];
    jmaterial["youngs_modulus"] = json_["youngs_modulus"];
    jmaterial["poisson_ratio"] = json_["poisson_ratio"];
    jmaterial["pc0"] = pc0;
    jmaterial["m"] = json_["m"];
    jmaterial["lambda"] = lambda;
    jmaterial["kappa"] = kappa;
    jmaterial["three_invariants"] = json_["three_invariants"];
    if (confining_p < (4000 * 6894.75729))
      jmaterial["bonding"] = false;
    else
      jmaterial["bonding"] = json_["bonding"];
    jmaterial["e0"] = e0;
    // Bonding parameters
    jmaterial["s_h"] = json_["s_h"];
    jmaterial["mc_a"] = json_["mc_a"];
    jmaterial["mc_b"] = json_["mc_b"];
    jmaterial["mc_c"] = json_["mc_c"];
    jmaterial["mc_d"] = json_["mc_d"];
    jmaterial["m_degradation"] = json_["m_degradation"];
    jmaterial["m_shear"] = json_["m_shear"];
    // Material id
    unsigned id = 0;
    // Create material
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "CamClay3D", std::move(id), jmaterial);
    // Pointer of the MCC
    auto cam_clay = std::make_shared<mpm::CamClay<Dim>>(id, jmaterial);
    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    // Initialise dpstrain
    mpm::Material<Dim>::Vector6d dpstrain;
    dpstrain.setZero();
    // Initialise dstress
    mpm::Material<Dim>::Vector6d dstress;
    dstress.setZero();
    // Initialise strain
    mpm::Material<Dim>::Vector6d strain;
    mpm::Material<Dim>::Vector6d pstrain;
    strain.setZero();
    pstrain.setZero();
    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    // Initial stress status
    stress.setZero();
    // Define the output files
    std::ofstream fPc, fStrain, fStress, fVoid;
    fPc.open("MCC_bonded_drained_pc_pcc_pcd_" + std::to_string(j + 1) + ".txt");
    fStrain.open("MCC_bonded_drained_strain_" + std::to_string(j + 1) + ".txt");
    fStress.open("MCC_bonded_drained_p_q_" + std::to_string(j + 1) + ".txt");
    fVoid.open("MCC_bonded_drained_void_ratio_" + std::to_string(j + 1) +
               ".txt");
    stress(0) = stress(1) = stress(2) = -confining_p;
    // Compute stress invariants
    mpm::dense_map state_vars = material->initialise_state_variables();
    // Initialise chi
    double chi = 1.;
    // Compute initial bonding parameters
    cam_clay->compute_bonding_parameters(chi, &state_vars);
    // Initialise vector n
    mpm::Material<Dim>::Vector6d n;
    n.setZero();
    // Computation loop
    for (int i = 0; i < nsteps; ++i) {
      // Compute elastic tensor
      cam_clay->compute_elastic_tensor(&state_vars);
      // Compute stress invariants
      cam_clay->compute_stress_invariants(stress, n, &state_vars);
      // Check yield
      typename mpm::CamClay<Dim>::FailureState yield =
          cam_clay->compute_yield_state(&state_vars);
      // Elastic status
      if (yield == mpm::CamClay<Dim>::FailureState::Elastic) {
        // Compute dstrain
        dstrain(0) = dstrain(1) = -cam_clay->de()(0, 2) * dstrain_axial /
                                  (cam_clay->de()(0, 0) + cam_clay->de()(0, 1));
        dstrain(2) = dstrain_axial;
        // Compute dstress_axial
        dstress(2) =
            (cam_clay->de()(2, 0) + cam_clay->de()(2, 1)) * dstrain(0) +
            cam_clay->de()(2, 2) * dstrain_axial;

      } else {
        // Compute plastic stiffness
        cam_clay->compute_plastic_tensor(stress, &state_vars);
        /*
        // Compute dstrain
        dstrain = (cam_clay->de() - cam_clay->dp()).inverse() * dstress;
        // Compute dpstrain
        dpstrain = dstrain - cam_clay->de().inverse() * dstress;
        // Compute pstrain
        pstrain += dpstrain;
        */
        // Compute dstrain
        dstrain(0) = dstrain(1) = -(cam_clay->de() - cam_clay->dp())(0, 2) *
                                  dstrain_axial /
                                  ((cam_clay->de() - cam_clay->dp())(0, 0) +
                                   (cam_clay->de() - cam_clay->dp())(0, 1));
        dstrain(2) = dstrain_axial;
        // Compute dstress_axial
        dstress(2) = ((cam_clay->de() - cam_clay->dp())(2, 0) +
                      (cam_clay->de() - cam_clay->dp())(2, 1)) *
                         dstrain(0) +
                     (cam_clay->de() - cam_clay->dp())(2, 2) * dstrain_axial;
        // Compute dpstrain
        dpstrain = dstrain - cam_clay->de().inverse() * dstress;
        // Compute pstrain
        pstrain += dpstrain;
        // Compute plastic volumetrix strain
        state_vars.at("dpvstrain") = (dpstrain(0) + dpstrain(1) + dpstrain(2));
        // Compute plastic deviatoric strain
        state_vars.at("dpdstrain") = 2. / 3. *
                                     sqrt(pow((dpstrain(0) - dpstrain(1)), 2) +
                                          pow((dpstrain(0) - dpstrain(2)), 2) +
                                          pow((dpstrain(2) - dpstrain(1)), 2));
        // Update pc
        state_vars.at("pc") *= exp(
            -state_vars.at("dpvstrain") * (1 + state_vars.at("void_ratio")) /
            (material->property("lambda") - material->property("kappa")));
        // Record chi at last step
        chi = state_vars.at("chi");
        // Update pcc and pcd
        cam_clay->compute_bonding_parameters(chi, &state_vars);
      }
      // Compute strain
      strain += dstrain;
      // Update void ratio
      state_vars.at("void_ratio") +=
          ((dstrain(0) + dstrain(1) + dstrain(2)) * (1 + e0));
      // Compute stress
      stress += dstress;
      // Write output
      if (i % output_steps == 0) {
        fPc << state_vars.at("pc") << "\t" << state_vars.at("pcc") << "\t"
            << state_vars.at("pcd") << "\n";
        fStrain << strain(0) << "\t" << strain(1) << "\t" << strain(2) << "\t"
                << state_vars.at("dpvstrain") << "\t"
                << state_vars.at("dpdstrain") << "\n";
        fStress << state_vars.at("p") << "\t" << state_vars.at("q") << "\n";
        fVoid << state_vars.at("void_ratio") << "\n";
      }
    }
    fPc.close();
    fStrain.close();
    fStress.close();
    fVoid.close();
  }
}