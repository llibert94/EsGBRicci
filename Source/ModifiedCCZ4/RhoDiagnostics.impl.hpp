/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(RHODIAGNOSTICS_HPP_)
#error "This file should only be included through RhoDiagnostics.hpp"
#endif

#ifndef RHODIAGNOSTICS_IMPL_HPP_
#define RHODIAGNOSTICS_IMPL_HPP_

template <class theory_t>
RhoDiagnostics<theory_t>::RhoDiagnostics(
    const theory_t a_theory, double dx,
    const std::array<double, CH_SPACEDIM> a_center)
    : my_theory(a_theory), m_center(a_center), m_deriv(dx)
{
}

template <class theory_t>
template <class data_t>
void RhoDiagnostics<theory_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<BSSNTheoryVars>();
    const auto d1 = m_deriv.template diff1<BSSNTheoryVars>(current_cell);
    const auto d2 = m_deriv.template diff2<BSSNTheoryVars>(current_cell);

    // Coordinates
    Coordinates<data_t> coords{current_cell, m_deriv.m_dx, m_center};

    // Compute diagnostics
    AllRhos<data_t> all_rhos = my_theory.compute_all_rhos(vars, d1, d2, coords);

    // Write the constraints into the output FArrayBox
    current_cell.store_vars(all_rhos.phi, c_rho_phi);
    current_cell.store_vars(all_rhos.g2, c_rho_g2);
    current_cell.store_vars(all_rhos.g3, c_rho_g3);
    current_cell.store_vars(all_rhos.GB, c_rho_GB);
    current_cell.store_vars(all_rhos.total, c_rho_total);
}

#endif /* RHODIAGNOSTICS_IMPL_HPP_ */
