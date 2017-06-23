/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_EOS_FORM_FACTORS_HQET_B_TO_C_HH
#define EOS_GUARD_EOS_FORM_FACTORS_HQET_B_TO_C_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class HQETBToC :
        public ParameterUser,
        public PrivateImplementationPattern<HQETBToC>
    {
        public:
            HQETBToC(const Parameters &, const Options &);
            ~HQETBToC();

            /*!
             * Matching coefficients for heavy-to-heavy b -> c vector currents, as
             * functions of the cusp omega.
             */
            // @{
            double c_1_vector(const double & omega) const;
            double c_2_vector(const double & omega) const;
            double c_3_vector(const double & omega) const;
            // @}

            /*!
             * Matching coefficients for heavy-to-heavy b -> c axial vector currents,
             * as functions of the cusp omega.
             */
            // @{
            double c_1_axialvector(const double & omega) const;
            double c_2_axialvector(const double & omega) const;
            double c_3_axialvector(const double & omega) const;
            // @}
    };
}

#endif
