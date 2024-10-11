/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "ModifiedGravitySimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBH.hpp"
#include "CouplingAndPotential.hpp"
#include "FourDerivScalarTensor.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedPunctureGauge.hpp"

class SimulationParameters : public ModifiedGravitySimulationParametersBase<
                                 FourDerivScalarTensor<CouplingAndPotential>>
{
  public:
    SimulationParameters(GRParmParse &pp)
        : ModifiedGravitySimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters
    void read_params(GRParmParse &pp)
    {

        // Do we want puncture tracking and constraint norm calculation?
        pp.load("track_punctures", track_punctures, false);
        pp.load("puncture_tracking_level", puncture_tracking_level, max_level);
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);

        // Coupling and potential
        pp.load("lambda_GB", coupling_and_potential_params.lambda_GB, 0.);
        pp.load("g2", coupling_and_potential_params.g2, 0.);
        pp.load("beta_ricci", coupling_and_potential_params.beta_ricci, 0.);
        pp.load("cutoff_GB", coupling_and_potential_params.cutoff_GB, 0.15);
        pp.load("factor_GB", coupling_and_potential_params.factor_GB, 100.);
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass, 0.);

        // Initial data
        pp.load("G_Newton", G_Newton, 1.0);

        pp.load("expand_matrix", expand_matrix, 0);

        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("scalar_amplitude", initial_params.amplitude, 0.);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_r0", initial_params.r0, 0.);

        // Initial BH data
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
        FOR(idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }

#ifdef USE_AHFINDER
        pp.load("AH_1_initial_guess", AH_1_initial_guess,
                0.5 * bh1_params.mass);
        pp.load("AH_2_initial_guess", AH_2_initial_guess,
                0.5 * bh2_params.mass);
        pp.load("AH_set_origins_to_punctures", AH_set_origins_to_punctures,
                false);
#endif
    }

    void check_params()
    {
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
        FOR(idir)
        {
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
    double G_Newton;
    int expand_matrix;
    InitialScalarData::params_t initial_params;
    CouplingAndPotential::params_t coupling_and_potential_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;

#ifdef USE_AHFINDER
    double AH_1_initial_guess;
    double AH_2_initial_guess;
    bool AH_set_origins_to_punctures;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP */
