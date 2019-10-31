#include <fstream>
#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/cam_clay.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

//! Check Undrained condition in 2D
TEST_CASE("Undrained condition is checked in 3D", "[material][cam_clay][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.275;
  jmaterial["p_ref"] = 100000;
  jmaterial["e_ref"] = 1.12;
  jmaterial["pc0"] = 200000;
  jmaterial["ocr"] = 1.;
  jmaterial["m"] = 1.2;
  jmaterial["lambda"] = 0.4;
  jmaterial["kappa"] = 0.05;
  jmaterial["three_invariants"] = false;
  jmaterial["bonding"] = false;

  SECTION("CamClay check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "CamClay3D", std::move(id), jmaterial);

    auto cam_clay = std::make_shared<mpm::CamClay<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000;
    stress(1) = -200000;
    stress(2) = -200000;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute stress invariants
    mpm::dense_map state_vars = material->initialise_state_variables();
    cam_clay->compute_stress_invariants(stress, &state_vars);

    // Compute elastic modulus
    cam_clay->compute_elastic_tensor(&state_vars);
    REQUIRE(state_vars.at("bulk_modulus") ==
            Approx(7370964.511104).epsilon(Tolerance));
    REQUIRE(state_vars.at("shear_modulus") ==
            Approx(3902275.329408).epsilon(Tolerance));

    // Initialise strain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = 0.00050000;
    dstrain(1) = 0.00050000;
    dstrain(2) = -0.00100000;
    dstrain(3) = 0.0000000;
    dstrain(4) = 0.0000000;
    dstrain(5) = 0.0000000;

    // Define the steps
    const int iter = 100;
    std::ofstream fStress, fVoid;
    fStress.open("mcc_undrained_p_q.txt");
    fVoid.open("mcc_undrained_void_ratio.txt");
    // Compute updated stress
    for (int i = 0; i < iter; i++) {
      // Compute stress
      stress = cam_clay->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      //! Initialise writing of inputfiles
      fStress << state_vars.at("p") << "\t" << state_vars.at("q") << "\n";
      fVoid << state_vars.at("void_ratio") << "\n";
    }
    fStress.close();
    fVoid.close();
    // Check stressees
    // REQUIRE(stress(0) == Approx(-200000).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(-200000).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(-200000).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));
  }
}

//! Check Drained condition in 2D
TEST_CASE("Drained condition is checked in 3D", "[material][cam_clay][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.275;
  jmaterial["p_ref"] = 100000;
  jmaterial["e_ref"] = 1.12;
  jmaterial["pc0"] = 200000;
  jmaterial["ocr"] = 1.;
  jmaterial["m"] = 1.2;
  jmaterial["lambda"] = 0.4;
  jmaterial["kappa"] = 0.05;
  jmaterial["three_invariants"] = false;
  jmaterial["bonding"] = false;

  SECTION("CamClay check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "CamClay3D", std::move(id), jmaterial);

    auto cam_clay = std::make_shared<mpm::CamClay<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000;
    stress(1) = -200000;
    stress(2) = -200000;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute stress invariants
    mpm::dense_map state_vars = material->initialise_state_variables();
    cam_clay->compute_stress_invariants(stress, &state_vars);

    // Compute elastic modulus
    cam_clay->compute_elastic_tensor(&state_vars);
    REQUIRE(state_vars.at("bulk_modulus") ==
            Approx(7370964.511104).epsilon(Tolerance));
    REQUIRE(state_vars.at("shear_modulus") ==
            Approx(3902275.329408).epsilon(Tolerance));

    // Initialise strain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.0001;

    // Define the steps
    const int iter = 100;
    std::ofstream fStress, fVoid;
    fStress.open("mcc_drained_p_q.txt");
    fVoid.open("mcc_drained_void_ratio.txt");
    // Compute updated stress
    for (int i = 0; i < iter; i++) {
      double G = state_vars.at("shear_modulus");
      double a1 = state_vars.at("bulk_modulus") + (4.0 / 3.0) * G;
      double a2 = state_vars.at("bulk_modulus") - (2.0 / 3.0) * G;
      dstrain(0) = a2 * dstrain(2) / (-a1 - a2);
      dstrain(1) = a2 * dstrain(2) / (-a1 - a2);
      // Compute stress
      stress = cam_clay->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);
      //! Initialise writing of inputfiles
      fStress << state_vars.at("p") << "\t" << state_vars.at("q") << "\n";
      fVoid << state_vars.at("void_ratio") << "\n";
    }
    fStress.close();
    fVoid.close();
  }
}
