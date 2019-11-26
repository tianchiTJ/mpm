#include <fstream>
#include <limits>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/cam_clay.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

//! Check Exxon data for 12-5
TEST_CASE("Exxon data 12-5", "[material][cam_clay][3D]") {
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
  jmaterial["p_ref"] = 1000;
  jmaterial["e_ref"] = 1.12;
  jmaterial["pc0"] = 1000;
  jmaterial["ocr"] = 1.;
  jmaterial["m"] = 1.2;
  jmaterial["lambda"] = 0.4;
  jmaterial["kappa"] = 0.05;
  jmaterial["three_invariants"] = false;
  jmaterial["bonding"] = true;
  // Bonding parameters
  jmaterial["s_h"] = 0.4;
  jmaterial["mc_a"] = 300;
  jmaterial["mc_b"] = 1.6;
  jmaterial["mc_c"] = 100;
  jmaterial["mc_d"] = 1.6;
  jmaterial["degradation"] = 1.;

  unsigned id = 0;
  auto material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "CamClay3D", std::move(id), jmaterial);

  auto cam_clay = std::make_shared<mpm::CamClay<Dim>>(id, jmaterial);

  // Initialise dstrain
  Eigen::Matrix<double, 3406, 1> dstrainv;
  Eigen::Matrix<double, 3406, 1> dstrainh;

  // Number of coordinate lines
  unsigned nlines = 0;

  std::fstream file;
  file.open("./strain.txt", std::ios::in);

  std::string line;

  while (std::getline(file, line)) {
    boost::algorithm::trim(line);
    std::istringstream istream(line);

    istream >> dstrainv[nlines] >> dstrainh[nlines];

    ++nlines;
  }
  file.close();

  // Initialise stress
  mpm::Material<Dim>::Vector6d stress;
  stress.setZero();
  stress(0) = -material->property("pc0") / material->property("ocr");
  stress(1) = -material->property("pc0") / material->property("ocr");
  stress(2) = -material->property("pc0") / material->property("ocr");

  // Compute stress invariants
  mpm::dense_map state_vars = material->initialise_state_variables();
  mpm::Material<Dim>::Vector6d n;
  cam_clay->compute_stress_invariants(stress, n, &state_vars);

  mpm::Material<Dim>::Vector6d dstrain;
  dstrain.setZero();

  // Define the steps
  std::ofstream fPc, fPcc, fPcd, fStrain, fStress, fVoid, fChi;
  fPc.open("mcc_bonding_drained_hardening_pc.txt");
  fPcc.open("mcc_bonding_drained_hardening_pcc.txt");
  fPcd.open("mcc_bonding_drained_hardening_pcd.txt");
  fStrain.open("mcc_bonding_drained_hardening_strain.txt");
  fStress.open("mcc_bonding_drained_hardening_p_q.txt");
  fVoid.open("mcc_bonding_drained_hardening_void_ratio.txt");
  fChi.open("mcc_bonding_drained_hardening_chi.txt");

  const unsigned nsteps = 3406;
  const unsigned output_steps = 100;
  // Compute updated stress
  for (int i = 0; i < nsteps; i++) {
    // Get dstrain
    dstrain(0) = dstrain(1) = -dstrainh[i];
    dstrain(2) = -dstrainv[i];
    // Compute stress
    stress =
        cam_clay->compute_stress(stress, dstrain, particle.get(), &state_vars);

    //! Initialise writing of inputfiles
    if (i % output_steps == 0) {
      fPc << state_vars.at("pc") << "\n";
      fStrain << dstrain(0) << "\t" << dstrain(2) << "\t"
              << state_vars.at("dpvstrain") << "\t"
              << state_vars.at("dpdstrain") << "\n";
      fStress << state_vars.at("p") << "\t" << state_vars.at("q") << "\n";
      fVoid << state_vars.at("void_ratio") << "\n";
      fChi << state_vars.at("chi") << "\n";
      fPcc << state_vars.at("pcc") << "\n";
      fPcd << state_vars.at("pcd") << "\n";
    }
  }
  fStress.close();
  fVoid.close();
  fPc.close();
  fPcc.close();
  fPcd.close();
  fStrain.close();
  fStress.close();
  fVoid.close();
  fChi.close();
}
