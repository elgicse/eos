/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Elena Graverini
 * Copyright (c) 2017 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_B_DECAYS_LAMBDAB_TO_LAMBDAC2625_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_LAMBDAB_TO_LAMBDAC2625_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: Lambda_b -> Lambda_c(2625) l nu
     */
    class LambdabToLambdac2625LeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<LambdabToLambdac2625LeptonNeutrino>
    {
        public:
            LambdabToLambdac2625LeptonNeutrino(const Parameters & parameters, const Options & options);
            ~LambdabToLambdac2625LeptonNeutrino();

            // [BBGIOvD] parametrization for the differential decay width
            double F12T(const double & s) const;
            double F120(const double & s) const;
            double F12P(const double & s) const;
            double F32P(const double & s) const;
            double G12T(const double & s) const;
            double G120(const double & s) const;
            double G12P(const double & s) const;
            double G32P(const double & s) const;
            double s_plus(const double & s) const;
            double s_minus(const double & s) const;
            
            double a_l(const double & s) const;
            double b_l(const double & s) const;
            double c_l(const double & s) const;
            double gamma_0(const double & s) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double double_differential_branching_ratio(const double & s, const double & theta_l) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;

            // R_Lambda_c_2625
            double differential_r_lambdac2625(const double & s) const;
            double integrated_r_lambdac2625() const;

            // Zero-Recoil Sum Rule
            double f_inel() const;
            double g_inel() const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_theta_l;
    };
}

#endif
