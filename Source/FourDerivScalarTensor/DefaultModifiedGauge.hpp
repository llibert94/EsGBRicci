/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTMODIFIEDGAUGE_HPP_
#define DEFAULTMODIFIEDGAUGE_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultModifiedGauge
{
  public:
    //! The constructor
    DefaultModifiedGauge() {}

    //! Set the modified gauge functions here to zero
    template <class data_t, template <typename> class coords_t>
    void compute_modified_gauge(data_t &a_of_x, data_t &b_of_x,
                           const coords_t<data_t> &coords) const
    {
        // a(x) from \tilde{g}^{\mu\nu} = g^{\mu\nu} - a(x)n^{\mu}n^{\nu}
        a_of_x = 0.;

        // b(x) from \hat{g}^{\mu\nu} = g^{\mu\nu} - b(x)n^{\mu}n^{\nu}
        b_of_x = 0.;
    }
};

#endif /* DEFAULTMODIFIEDGAUGE_HPP_ */