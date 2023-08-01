/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ESGB_HPP_
#define ESGB_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultCoupling.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is an example of a matter_t object which calculates the
     matter type specific elements for the RHS update and the evaluation
     of the constraints. This includes the Energy Momentum Tensor, and
     the matter evolution terms. In this case, a scalar field coupled to
     the Gauss-Bonnet invariant, the matter elements are phi and its
     curvature, Kphi.
     It is templated over a coupling function coupling_t which the
     user must specify in a class, although a default is provided which
     sets dfGB and df2GB to zero.
     It assumes a scalar field without potential, coming from the action
     S=\frac{1}{8\pi G}\int d^4x\left(R - 1/2\nabla_{\mu}\nabla^{\mu}\phi
     +f(\phi)/4{\mathcal L}_{GB}\right)
     \sa MatterCCZ4(), ConstraintsMatter()
*/
template <class data_t> struct rho_Si_t
{
    Tensor<1, data_t> Si; //!< S_i = T_ia_n^a
    data_t rho;           //!< rho = T_ab n^a n^b
};

template <class coupling_t = DefaultCoupling> class EsGB
{
  protected:
    //! The local copy of the coupling
    coupling_t my_coupling;

  public:
    //!  Constructor of class EsGB, inputs are the matter parameters.
    EsGB(const coupling_t a_coupling) : my_coupling(a_coupling) {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t phi;  // the scalar field
        data_t Kphi; // the curvature of the scalar field

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
            VarsTools::define_enum_mapping(mapping_function, c_Kphi, Kphi);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
        }
    };

    //! The function which calculates the rho and Si elements of the EM
    //! Tensor, given the vars and derivatives, for the given coupling
    //! function
    //! NOTE: this is computed in a separate function from the other
    //! elements of the EM, namely S and Sij, so that we can call directly
    //! this function to calculate the constraint without needing to
    //! specify the gauge variables
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    rho_Si_t<data_t> compute_rho_Si(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2 //!< the value of the 2nd derivs
    ) const;

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, for the given coupling function
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2, //!< the value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

    //! The function which adds in the RHS for the matter field vars,
    //! given the coupling function
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &rhs,    //!< the value of the RHS for all vars
        const vars_t<data_t> &vars, //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2, //!< the value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

    //! The function which computes the LHS matrix for some of the vars, which
    //! are A, K and Kphi for this EsGB example
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void compute_lhs(
        const int N,
        data_t *LHS, //!< the LHS matrix for some vars (namely A, K and Kphi)
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2, //!< the value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

    //! The function which solves the linear system using the LHS matrix
    //! computed in compute_lhs and the RHS calculated before
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void solve_lhs(
        rhs_vars_t<data_t> &rhs,             //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms
};

//! Class for the function that solves the linear system (included in
//! lapack.cpp) A separate class has been needed to be defined because of the
//! explicit specialization of this function, which does not work for a class
//! with an extra template argument (namely the coupling for the EsGB class)
class EsGB_lapack
{
  public:
    EsGB_lapack() {}
    template <class data_t>
    static void solve_system_lapack(const int N, data_t *LHS, data_t *RHS);
};

#include "EsGB.impl.hpp"

#endif /* ESGB_HPP_ */
