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

//! Check undrained condition for hardening in 3D
TEST_CASE("Undrained condition hardening is checked in 3D",
          "[material][cam_clay][3D]") {
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
    mpm::Material<Dim>::Vector6d n;
    cam_clay->compute_stress_invariants(stress, n, &state_vars);

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
    fStress.open("mcc_undrained_hardening_p_q.txt");
    fVoid.open("mcc_undrained_hardening_void_ratio.txt");
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

//! Check undrained condition for softening in 3D
TEST_CASE("Undrained condition softening is checked in 3D",
          "[material][cam_clay][3D]") {
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
  jmaterial["ocr"] = 4.;
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
    stress(0) = -50000;
    stress(1) = -50000;
    stress(2) = -50000;
    REQUIRE(stress(0) == Approx(-50000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-50000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-50000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute stress invariants
    mpm::dense_map state_vars = material->initialise_state_variables();
    mpm::Material<Dim>::Vector6d n;
    cam_clay->compute_stress_invariants(stress, n, &state_vars);

    // Compute elastic modulus
    // cam_clay->compute_elastic_tensor(&state_vars);
    // REQUIRE(state_vars.at("bulk_modulus") ==
    //        Approx(7370964.511104).epsilon(Tolerance));
    // REQUIRE(state_vars.at("shear_modulus") ==
    //        Approx(3902275.329408).epsilon(Tolerance));

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
    const int iter = 200;
    std::ofstream fStress, fVoid;
    fStress.open("mcc_undrained_softening_p_q.txt");
    fVoid.open("mcc_undrained_softening_void_ratio.txt");
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

//! Check drained condition for hardening in 3D
TEST_CASE("drained condition hardening is checked in 3D",
          "[material][cam_clay][3D]") {
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
  jmaterial["pc0"] = 300000;
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
    stress(0) = -300000;
    stress(1) = -300000;
    stress(2) = -300000;
    REQUIRE(stress(0) == Approx(-300000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-300000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-300000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute stress invariants
    mpm::dense_map state_vars = material->initialise_state_variables();
    mpm::Material<Dim>::Vector6d n;
    cam_clay->compute_stress_invariants(stress, n, &state_vars);

    // Compute elastic modulus
    cam_clay->compute_elastic_tensor(&state_vars);
    REQUIRE(state_vars.at("bulk_modulus") ==
            Approx(7370964.511104).epsilon(Tolerance));
    REQUIRE(state_vars.at("shear_modulus") ==
            Approx(3902275.329408).epsilon(Tolerance));

    // Initialise dstrain
    Eigen::Matrix<double, 262, 1> dstrainv;
    dstrainv << -0.025, 0.029, 0.032, 0.033, 0.03, 0.034, 0.03, 0.032, 0.031,
        0.033, 0.03, 0.03, 0.032, 0.031, 0.032, 0.033, 0.03, 0.031, 0.032,
        0.032, 0.031, 0.029, 0.033, 0.031, 0.031, 0.032, 0.03, 0.032, 0.031,
        0.032, 0.032, 0.031, 0.032, 0.03, 0.031, 0.033, 0.031, 0.032, 0.03,
        0.033, 0.032, 0.026, 0.029, 0.032, 0.032, 0.032, 0.032, 0.032, 0.031,
        0.031, 0.033, 0.03, 0.03, 0.031, 0.031, 0.032, 0.032, 0.03, 0.031,
        0.031, 0.032, 0.031, 0.03, 0.03, 0.032, 0.031, 0.059, 0.057, 0.058,
        0.057, 0.057, 0.058, 0.058, 0.055, 0.056, 0.055, 0.058, 0.056, 0.054,
        0.057, 0.058, 0.058, 0.056, 0.055, 0.057, 0.057, 0.057, 0.056, 0.058,
        0.056, 0.058, 0.059, 0.058, 0.057, 0.058, 0.058, 0.058, 0.058, 0.059,
        0.058, 0.059, 0.06, 0.057, 0.06, 0.058000000000001, 0.058, 0.059, 0.058,
        0.058, 0.058, 0.057, 0.058, 0.057, 0.058, 0.056, 0.058, 0.055, 0.057,
        0.058, 0.109, 0.1, 0.102, 0.106, 0.102, 0.104, 0.101999999999999,
        0.103000000000001, 0.1, 0.104, 0.103, 0.105, 0.101999999999999,
        0.108000000000001, 0.104, 0.101, 0.103, 0.102, 0.105, 0.103,
        0.103000000000001, 0.103, 0.104, 0.108, 0.103, 0.104, 0.105, 0.105,
        0.102000000000001, 0.103999999999999, 0.105, 0.106999999999999, 0.106,
        0.1, 0.107000000000001, 0.103999999999999, 0.105, 0.105, 0.102,
        0.100999999999999, 0.103, 0.102, 0.108000000000001, 0.103999999999999,
        0.101000000000001, 0.106, 0.105, 0.107999999999999, 0.104000000000001,
        0.102, 0.103999999999999, 0.105, 0.103999999999999, 0.107000000000001,
        0.1, 0.103, 0.103999999999999, 0.103, 0.105, 0.108000000000001, 0.103,
        0.106999999999999, 0.103000000000002, 0.106999999999999, 0.109,
        0.103999999999999, 0.107000000000001, 0.105, 0.107999999999999, 0.102,
        0.104000000000001, 0.107999999999999, 0.104000000000001,
        0.106999999999999, 0.108000000000001, 0.103, 0.103999999999999, 0.102,
        0.103, 0.109, 0.104000000000001, 0.102, 0.106, 0.103999999999999, 0.106,
        0.101000000000001, 0.106, 0.103999999999999, 0.108000000000001, 0.106,
        0.102, 0.105, 0.104999999999999, 0.104000000000001, 0.105,
        0.109999999999999, 0.100999999999999, 0.103, 0.104000000000001,
        0.103999999999999, 0.106, 0.105, 0.104000000000001, 0.106, 0.106, 0.102,
        0.105999999999998, 0.106000000000002, 0.106999999999999,
        0.106999999999999, 0.101000000000003, 0.103999999999999, 0.105,
        0.102999999999998, 0.105, 0.108000000000001, 0.102999999999998, 0.102,
        0.102, 0.104000000000003, 0.099999999999998, 0.102, 0.105,
        0.103000000000002, 0.107999999999997, 0.105, 0.108000000000001,
        0.106999999999999, 0.105, 0.106000000000002, 0.108000000000001,
        0.101999999999997, 0.104000000000003, 0.102, 0.103999999999999,
        0.103999999999999, 0.106000000000002, 0.100999999999999,
        0.108000000000001, 0.102999999999998, 0.100000000000001,
        0.099999999999998, 0.102;

    Eigen::Matrix<double, 262, 1> dstrainh;
    dstrainh << 0.008333333333333, -0.009666666666667, -0.011,
        -0.011333333333333, -0.011, -0.015, -0.013666666666667,
        -0.014666666666667, -0.014666666666667, -0.016, -0.014333333333333,
        -0.014666666666667, -0.014666666666667, -0.014666666666667,
        -0.014666666666667, -0.015, -0.013666666666667, -0.013666666666667,
        -0.014333333333333, -0.014, -0.013333333333333, -0.012333333333333,
        -0.013666666666667, -0.013, -0.012666666666667, -0.012666666666667,
        -0.012, -0.012333333333333, -0.012, -0.012333333333333, -0.012,
        -0.011333333333333, -0.011666666666667, -0.011, -0.011,
        -0.011666666666667, -0.010666666666667, -0.011, -0.01, -0.011,
        -0.010666666666667, -0.008666666666667, -0.009666666666667,
        -0.010333333333333, -0.01, -0.010333333333333, -0.009666666666667,
        -0.01, -0.009333333333333, -0.009333333333333, -0.01,
        -0.008666666666667, -0.008666666666667, -0.009333333333333, -0.009,
        -0.009, -0.009333333333333, -0.008333333333333, -0.008666666666667,
        -0.008333333333333, -0.009, -0.008333333333333, -0.008333333333333,
        -0.007666666666667, -0.008666666666667, -0.008333333333333,
        -0.015666666666667, -0.014666666666667, -0.015, -0.014333333333333,
        -0.014, -0.014333333333333, -0.014333333333333, -0.013,
        -0.013333333333333, -0.013, -0.013666666666667, -0.013,
        -0.012333333333333, -0.013, -0.013333333333334, -0.013,
        -0.012666666666667, -0.012, -0.012666666666667, -0.012333333333333,
        -0.012666666666667, -0.012, -0.012666666666667, -0.012,
        -0.012333333333333, -0.012666666666667, -0.012333333333334, -0.012,
        -0.012, -0.012333333333334, -0.012, -0.012, -0.012333333333333, -0.012,
        -0.012, -0.012666666666667, -0.011666666666667, -0.012666666666666,
        -0.011666666666667, -0.011666666666667, -0.011666666666667,
        -0.011666666666667, -0.011666666666667, -0.011666666666667, -0.011,
        -0.011333333333333, -0.011, -0.011333333333334, -0.011, -0.011,
        -0.010333333333333, -0.011, -0.011, -0.020666666666666,
        -0.018666666666667, -0.019, -0.02, -0.019, -0.019666666666667, -0.019,
        -0.019, -0.019, -0.019666666666666, -0.019, -0.02, -0.019,
        -0.020666666666667, -0.019666666666667, -0.018666666666667,
        -0.019333333333333, -0.019333333333333, -0.019666666666667, -0.019,
        -0.019333333333334, -0.019333333333333, -0.02, -0.020333333333333,
        -0.019333333333334, -0.019666666666666, -0.020333333333334,
        -0.020333333333333, -0.019666666666667, -0.02, -0.02,
        -0.020666666666667, -0.020666666666666, -0.019, -0.021, -0.02,
        -0.020333333333333, -0.020333333333334, -0.019666666666667,
        -0.019333333333333, -0.02, -0.02, -0.021333333333334,
        -0.020333333333333, -0.019666666666667, -0.021333333333333,
        -0.020666666666667, -0.021333333333333, -0.021, -0.020333333333333,
        -0.020666666666666, -0.021, -0.020666666666667, -0.021666666666667,
        -0.019666666666666, -0.020666666666667, -0.021, -0.020666666666667,
        -0.021666666666667, -0.022, -0.021333333333334, -0.022333333333333,
        -0.021000000000001, -0.022666666666666, -0.023, -0.022333333333334,
        -0.023, -0.022, -0.023333333333333, -0.022, -0.022333333333334,
        -0.023666666666666, -0.022, -0.023666666666667, -0.023333333333334,
        -0.022, -0.022666666666666, -0.022, -0.022333333333334, -0.024, -0.023,
        -0.022666666666667, -0.023666666666667, -0.023333333333332,
        -0.024000000000001, -0.022666666666667, -0.024333333333333,
        -0.023333333333334, -0.025, -0.025, -0.021666666666667,
        -0.024333333333333, -0.024333333333334, -0.024, -0.024333333333333,
        -0.025666666666667, -0.023666666666666, -0.024333333333334,
        -0.024333333333333, -0.025, -0.025, -0.025333333333334, -0.025, -0.026,
        -0.026, -0.025333333333334, -0.026666666666666, -0.026666666666667,
        -0.027, -0.027333333333333, -0.025666666666667, -0.027666666666667,
        -0.027000000000001, -0.025999999999999, -0.026666666666667,
        -0.027666666666667, -0.026333333333333, -0.026333333333334,
        -0.026333333333333, -0.026666666666667, -0.026, -0.026333333333334,
        -0.027666666666667, -0.027, -0.028666666666665, -0.027666666666668,
        -0.028666666666666, -0.028666666666667, -0.028000000000001,
        -0.028666666666667, -0.029333333333333, -0.027666666666666,
        -0.029666666666667, -0.029333333333334, -0.029666666666667,
        -0.029999999999999, -0.030333333333334, -0.029, -0.030666666666667,
        -0.029333333333333, -0.028333333333334, -0.028666666666665,
        -0.029000000000001;

    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();

    // Define the steps
    const int iter = 262;
    std::ofstream fStress, fVoid;
    fStress.open("mcc_drained_hardening_p_q.txt");
    fVoid.open("mcc_drained_hardening_void_ratio.txt");
    // Compute updated stress
    for (int i = 0; i < iter; i++) {
      // Get dstrain
      dstrain(0) = dstrain(1) = -dstrainh(i) / 100;
      dstrain(2) = -dstrainv(i) / 100;
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
    // REQUIRE(stress(0) == Approx(-300000).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(-300000).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(-300000).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));
  }
}
