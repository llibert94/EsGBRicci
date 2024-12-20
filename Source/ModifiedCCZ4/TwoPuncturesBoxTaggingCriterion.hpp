/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TWOPUNCTURESBOXTAGGINGCRITERION_HPP_
#define TWOPUNCTURESBOXTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

class TwoPuncturesBoxTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const int m_max_level;
    const std::vector<double> m_puncture_masses;
    const std::vector<std::array<double, CH_SPACEDIM>> &m_puncture_coords;

  public:
    TwoPuncturesBoxTaggingCriterion(
        const double dx, const int a_level, const int a_max_level,
        const std::vector<std::array<double, CH_SPACEDIM>> &a_puncture_coords,
        const std::vector<double> a_puncture_masses = {0.5, 0.5})
        : m_dx(dx), m_level(a_level), m_max_level(a_max_level),
          m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0;

        double puncture_separation_squared = 0.0;
        FOR(i)
        {
            double displacement =
                m_puncture_coords[0][i] - m_puncture_coords[1][i];
            puncture_separation_squared += displacement * displacement;
        }
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        // we want each level to be double the innermost one in size
        const double factor = pow(2.0, m_max_level - m_level - 1);
        double sum_masses = m_puncture_masses[0] + m_puncture_masses[1];

        if (puncture_separation_squared > sum_masses * sum_masses)
        {
            // loop over puncture masses
            for (int ipuncture = 0; ipuncture < m_puncture_masses.size();
                 ++ipuncture)
            {
                // where am i?
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_puncture_coords[ipuncture]);
                const data_t max_abs_xy =
                    simd_max(abs(coords.x), abs(coords.y));
                const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
                auto regrid = simd_compare_lt(
                    max_abs_xyz, 1.5 * factor * m_puncture_masses[ipuncture]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }
        else
        {
            // where am i?
            const std::array<double, CH_SPACEDIM> puncture_centre;
            FOR(i)
            puncture_centre[i] =
                0.5 * (m_puncture_coords[0][i] + m_puncture_coords[1][i]);
            const Coordinates<data_t> coords(current_cell, m_dx,
                                             puncture_centre);
            const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
            const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
            auto regrid =
                simd_compare_lt(max_abs_xyz, 1.5 * factor * sum_masses);
            criterion = simd_conditional(regrid, 100.0, criterion);
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* TWOPUNCTURESBOXTAGGINGCRITERION_HPP_ */
