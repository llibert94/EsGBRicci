/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBH4DSTLEVEL_HPP_
#define BINARYBH4DSTLEVEL_HPP_

#include "CouplingAndPotential.hpp"
#include "DefaultLevelFactory.hpp"
#include "FourDerivScalarTensor.hpp"
#include "GRAMRLevel.hpp"
#include "ModifiedPunctureGauge.hpp"
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp"

class BinaryBH4dSTLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BinaryBH4dSTLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);
#ifdef USE_TWOPUNCTURES
    TPAMR &m_tp_amr = dynamic_cast<TPAMR &>(m_gr_amr);
#endif /* USE_TWOPUNCTURES */

    // Typedef for 4dST
    typedef FourDerivScalarTensor<CouplingAndPotential>
        FourDerivScalarTensorWithCouplingAndPotential;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    /// Identify and tag the cells that need higher resolution
    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

#ifdef CH_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* CH_USE_HDF5 */
};

#endif /* BINARYBH4DSTLEVEL_HPP_ */
