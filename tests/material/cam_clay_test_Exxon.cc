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

//! Check drained condition with bonding for hardening in 3D
TEST_CASE("drained condition hardening with bonding is checked in 3D",
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
  jmaterial["bonding"] = true;
  // Bonding parameters
  jmaterial["s_h"] = 0.4;
  jmaterial["mc_a"] = 30000;
  jmaterial["mc_b"] = 1.6;
  jmaterial["mc_c"] = 100000;
  jmaterial["mc_d"] = 1.6;
  jmaterial["degradation"] = 1.;

  unsigned id = 0;
  auto material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "CamClay3D", std::move(id), jmaterial);

  auto cam_clay = std::make_shared<mpm::CamClay<Dim>>(id, jmaterial);

  // Number of steps
  const unsigned long long nsteps = 200000;
  // Interval of output step;
  const unsigned output_steps = 1000;
  // Initialise dstrain
  std::vector<double> dstrainv;
  std::vector<double> dstrainh;

  // Initialise dstrain
  const double dstrainq = 5.E-6;
  double dstrainp = 0;
  // Initialise stress invariances
  double p = material->property("pc0") / material->property("ocr");
  double q = 0;
  // Compute e0
  const double e0 =
      material->property("e_ref") -
      material->property("lambda") *
          log(material->property("pc0") / material->property("p_ref")) -
      material->property("kappa") * log(material->property("ocr"));
  // Initialise e
  double e = e0;
  // Initialise pc
  double pc = material->property("pc0");
  // Initialise bonded parameters
  double pcc = material->property("mc_a") *
               pow(material->property("s_h"), material->property("mc_b"));
  double pcd = material->property("mc_c") *
               pow(material->property("s_h"), material->property("mc_d"));
  double chi = 1.;

  for (int i = 0; i < nsteps; ++i) {
    double yield_function = q * q + pow(material->property("m"), 2) *
                                        (p + pcc) * (p - pc - pcc - pcd);
    // Initialise modulus
    double e_b = (1 + e) / material->property("kappa") * p;
    double e_s = 1.5 * (1 - 2 * material->property("poisson_ratio")) /
                 (1 + material->property("poisson_ratio")) * e_b;
    if (yield_function < 0) {
      // Compute dp, dq
      const double dq = 3 * e_s * dstrainq;
      const double dp = dq / 3;
      // Update dstrain
      dstrainp = dp / e_b;
      dstrainv.push_back(dstrainq + dstrainp / 3);
      dstrainh.push_back((dstrainp * 2 / 3 - dstrainq) / 2);
      // Update p, q
      p += dp;
      q += dq;
    } else {
      // Compute derivatives
      const double eta = q / p;
      const double df_dq = 2 * q;
      const double df_dp =
          pow(material->property("m"), 2) * (2 * p - (pc + pcd));
      const double hardening =
          pow(p, 2) * (pow(material->property("m"), 4) - pow(eta, 4)) *
              material->property("kappa") /
              (material->property("lambda") - material->property("kappa")) *
              e_b +
          (-2 * pcc - pc - pcd) *
              (-material->property("degradation") * material->property("mc_a") *
               material->property("mc_b") *
               pow(chi * material->property("s_h"),
                   material->property("mc_b"))) +
          (-p - pcc) *
              (-material->property("degradation") * material->property("mc_c") *
               material->property("mc_d") *
               pow(chi * material->property("s_h"),
                   material->property("mc_d")));
      const double delta_phi =
          (e_b * df_dp * dstrainp + 3 * e_s * df_dq * dstrainq) /
          (hardening + e_b * df_dp * df_dp + 3 * e_s * df_dq * df_dq);
      // Compute dp, dq
      const double dq = 3 * e_s * (dstrainq - delta_phi * df_dq);
      const double dp = dq / 3;
      // Update dstrain
      dstrainp = dp / e_b + delta_phi * df_dp;
      dstrainv.push_back(dstrainq + dstrainp / 3);
      dstrainh.push_back((dstrainp * 2 / 3 - dstrainq) / 2);
      // Update p, q
      p += dp;
      q += dq;
      // Update e
      e += -dstrainp * (1 + e0);
      // Update pc
      pc *= exp((1 + e) /
                (material->property("lambda") - material->property("kappa")) *
                delta_phi * (2 * p - pc));
      // Bonded parameters
      chi -= material->property("degradation") * chi *
             (delta_phi * (sqrt(6) * q / pow(material->property("m"), 2)));
      pcd = material->property("mc_a") *
            pow(chi * material->property("s_h"), material->property("mc_b"));
      pcc = material->property("mc_c") *
            pow(chi * material->property("s_h"), material->property("mc_d"));
    }
  }

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
