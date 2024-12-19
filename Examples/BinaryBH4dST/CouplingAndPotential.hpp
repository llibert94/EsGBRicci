/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "Coordinates.hpp"
#include "simd.hpp"

class CouplingAndPotential
{
  public:
    struct params_t
    {
        double lambda_GB; // Gauss-Bonnet coupling
        double g2;        // coupling to the square of the kinetic term
        double beta_ricci;
        double cutoff_GB;   // cutoff for switching off the Gauss-Bonnet terms
                            // inside the BH
        double factor_GB;   // factor for the function smoothening the GB cutoff
	double scalar_mass; // mass in the potential
	double quadratic_factor; // see arxiv 2306.14966 eq 65
	double quartic_factor;   // see arxiv 2306.14966 eq 65 
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    CouplingAndPotential(params_t a_params) : m_params(a_params) {}

    //! Set the EsGB coupling function fhere
    template <class data_t, template <typename> class vars_t>
    void compute_coupling_and_potential(data_t &dfdphi, data_t &d2fdphi2,
                                        data_t &g2, data_t &dg2dphi,
                                        data_t &beta, data_t &dbetadphi,
                                        data_t &d2betadphi2, data_t &V_of_phi,
                                        data_t &dVdphi,
                                        const vars_t<data_t> &vars,
                                        const Coordinates<data_t> &coords) const
    {
        // excision setting the coupling to 0 in the interior of the BH with a
        // smooth function
        data_t cutoff_factor =
            1. + exp(-m_params.factor_GB * (vars.chi - m_params.cutoff_GB));

        // Shift-symmetric coupling: f(\phi) = \lambda^{GB}\phi

        // The first derivative of the GB coupling function
        //dfdphi = m_params.lambda_GB / cutoff_factor * vars.phi;
	  dfdphi = m_params.lambda_GB / cutoff_factor *
                 exp(-m_params.quadratic_factor * vars.phi * vars.phi *
                     (1. + m_params.quartic_factor * vars.phi * vars.phi)) *
                 vars.phi *
                 (1. + 2. * m_params.quartic_factor * vars.phi * vars.phi);
        // The second derivative of the GB coupling function
        //d2fdphi2 = m_params.lambda_GB / cutoff_factor;
	  d2fdphi2  = m_params.lambda_GB / cutoff_factor *
            exp(-m_params.quadratic_factor * vars.phi * vars.phi *
                (1. + m_params.quartic_factor * vars.phi * vars.phi)) *
            (1. + 3. * m_params.quartic_factor * vars.phi * vars.phi -
             2. * m_params.quadratic_factor * vars.phi * vars.phi *
                 (1. + 2. * m_params.quartic_factor * vars.phi * vars.phi) *
                 (1. + 2. * m_params.quartic_factor * vars.phi * vars.phi));

        beta = 0.5 * m_params.beta_ricci * vars.phi * vars.phi;
        dbetadphi = m_params.beta_ricci * vars.phi;
        d2betadphi2 = m_params.beta_ricci;

        // The coupling to the square of the kinetic term
        g2 = m_params.g2;
        // The first derivative of the g2 coupling
        dg2dphi = 0.;
        // The potential of the scalar field
        V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);
        // The first derivative of the potential
        dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
    }
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
