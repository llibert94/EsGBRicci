/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBH.hpp"
#include "CouplingAndPotential.hpp"
#include "FourDerivScalarTensor.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedPunctureGauge.hpp"

class SimulationParameters : public SimulationParametersBase {
public:
  SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp) {
    read_params(pp);
    check_params();
  }

  /// Read parameters
  void read_params(GRParmParse &pp) {

    // Do we want puncture tracking and constraint norm calculation?
    pp.load("track_punctures", track_punctures, false);
    pp.load("puncture_tracking_level", puncture_tracking_level, max_level);
    pp.load("calculate_constraint_norms", calculate_constraint_norms, false);

    // Initial scalar field data
    initial_params.center = center; // already read in SimulationParametersBase
    pp.load("G_Newton", G_Newton, 1.0);
    pp.load("scalar_amplitude", initial_params.amplitude, 0.);
    pp.load("scalar_width", initial_params.width, 1.0);
    pp.load("scalar_r0", initial_params.r0, 0.);

    // Coupling and potential
    pp.load("lambda_GB", coupling_and_potential_params.lambda_GB, 0.);
    pp.load("g2", coupling_and_potential_params.g2, 0.);
    pp.load("cutoff_GB", coupling_and_potential_params.cutoff_GB, 0.15);
    pp.load("factor_GB", coupling_and_potential_params.factor_GB, 100.);
    pp.load("scalar_mass", coupling_and_potential_params.scalar_mass, 0.);

    // Modified gauge
    pp.load("a0", modified_ccz4_params.a0, 0.);
    pp.load("b0", modified_ccz4_params.b0, 0.);
    pp.load("lapse_advec_coeff", modified_ccz4_params.lapse_advec_coeff, 0.);
    pp.load("lapse_power", modified_ccz4_params.lapse_power, 1.);
    pp.load("lapse_coeff", modified_ccz4_params.lapse_coeff, 2.);
    pp.load("shift_Gamma_coeff", modified_ccz4_params.shift_Gamma_coeff, 0.75);
    pp.load("shift_advec_coeff", modified_ccz4_params.shift_advec_coeff, 0.);
    pp.load("eta", modified_ccz4_params.eta, 1.);
    modified_ccz4_params.kappa1 = ccz4_base_params.kappa1;
    modified_ccz4_params.kappa2 = ccz4_base_params.kappa2;
    modified_ccz4_params.kappa3 = ccz4_base_params.kappa3;
    modified_ccz4_params.covariantZ4 = ccz4_base_params.covariantZ4;

    // Initial data
    pp.load("massA", bh1_params.mass);
    pp.load("momentumA", bh1_params.momentum);
    pp.load("massB", bh2_params.mass);
    pp.load("momentumB", bh2_params.momentum);

    // Get the centers of the BHs either explicitly or as
    // an offset (not both, or they will be offset from center
    // provided)
    std::array<double, CH_SPACEDIM> centerA, centerB;
    std::array<double, CH_SPACEDIM> offsetA, offsetB;
    pp.load("centerA", centerA, center);
    pp.load("centerB", centerB, center);
    pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
    pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
    FOR(idir) {
      bh1_params.center[idir] = centerA[idir] + offsetA[idir];
      bh2_params.center[idir] = centerB[idir] + offsetB[idir];
    }
  }

  void check_params() {
    check_parameter("a(x)", modified_ccz4_params.a0,
                    modified_ccz4_params.a0 > -1, "should be >-1");
    warn_parameter("a(x)", modified_ccz4_params.a0,
                   modified_ccz4_params.a0 != 0, "should be !=0");
    warn_parameter("b(x)", modified_ccz4_params.b0,
                   (modified_ccz4_params.b0 > 0) &&
                       (modified_ccz4_params.b0 != modified_ccz4_params.a0),
                   "should be >0 and !=a(x)");
    check_parameter("kappa1", modified_ccz4_params.kappa1,
                    modified_ccz4_params.kappa1 > 0, "should be > 0");
    check_parameter("kappa2", modified_ccz4_params.kappa2,
                    modified_ccz4_params.kappa2 >
                        -2. / (2. + modified_ccz4_params.b0),
                    "should be > -2/(2+b(x))");
    warn_parameter("massA", bh1_params.mass, bh1_params.mass >= 0,
                   "should be >= 0");
    warn_parameter("massB", bh2_params.mass, bh2_params.mass >= 0,
                   "should be >= 0");
    warn_array_parameter(
        "momentumA", bh1_params.momentum,
        std::sqrt(ArrayTools::norm2(bh1_params.momentum)) <
            0.3 * bh1_params.mass,
        "approximation used for boosted BH only valid for small boosts");
    warn_array_parameter(
        "momentumB", bh2_params.momentum,
        std::sqrt(ArrayTools::norm2(bh2_params.momentum)) <
            0.3 * bh1_params.mass,
        "approximation used for boosted BH only valid for small boosts");
    FOR(idir) {
      std::string nameA = "centerA[" + std::to_string(idir) + "]";
      std::string nameB = "centerB[" + std::to_string(idir) + "]";
      double center_A_dir = bh1_params.center[idir];
      double center_B_dir = bh2_params.center[idir];
      warn_parameter(nameA, center_A_dir,
                     (center_A_dir >= 0.0) &&
                         (center_A_dir <= (ivN[idir] + 1) * coarsest_dx),
                     "should be within the computational domain");
      warn_parameter(nameB, center_B_dir,
                     (center_B_dir >= 0.0) &&
                         (center_B_dir <= (ivN[idir] + 1) * coarsest_dx),
                     "should be within the computational domain");
    }
    check_parameter("puncture_tracking_level", puncture_tracking_level,
                    (puncture_tracking_level >= 0) &&
                        (puncture_tracking_level <= max_level),
                    "must be between 0 and max_level (inclusive)");
  }

  bool track_punctures, calculate_constraint_norms;
  int puncture_tracking_level;

  // Collection of parameters necessary for initial conditions
  // Set these even in the case of TwoPunctures as they are used elsewhere
  // e.g. for puncture tracking/tagging
  double G_Newton;
  InitialScalarData::params_t initial_params;
  CouplingAndPotential::params_t coupling_and_potential_params;
  BoostedBH::params_t bh2_params;
  BoostedBH::params_t bh1_params;
  ModifiedCCZ4RHS<
      FourDerivScalarTensor<CouplingAndPotential>, ModifiedPunctureGauge,
      FourthOrderDerivatives>::modified_params_t modified_ccz4_params;
};

#endif /* SIMULATIONPARAMETERS_HPP */
