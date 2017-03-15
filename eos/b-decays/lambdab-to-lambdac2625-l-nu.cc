/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Elena Graverini
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

#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

#include <math.h>
// Only for debug purposes
#include <iostream>

using namespace std;

namespace eos
{
    template <>
    struct Implementation<LambdabToLambdac2625LeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>> form_factors;

        Parameters parameters;

        UsedParameter m_Lambdab;

        UsedParameter tau_Lambdab;

        UsedParameter m_Lambdac2625;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            m_Lambdab(p["mass::Lambda_b"], u),
            tau_Lambdab(p["life_time::Lambda_b"], u),
            m_Lambdac2625(p["mass::Lambda_c(2625)"], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {

            form_factors = FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda_c(2625)@BBGIOvD2017", p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);
        }

        // [BBGIOvD] parametrization for the differential decay width
        double a_l(const double & s) const
        {
            const auto F12T = form_factors->f_time12_v(s);
            const auto F120 = form_factors->f_long12_v(s);
            const auto F12P = form_factors->f_perp12_v(s);
            const auto F32P = form_factors->f_perp32_v(s);
            const auto G12T = form_factors->f_time12_a(s);
            const auto G120 = form_factors->f_long12_a(s);
            const auto G12P = form_factors->f_perp12_a(s);
            const auto G32P = form_factors->f_perp32_a(s);
            const double s_plus = pow((m_Lambdab + m_Lambdac2625), 2) - pow(s, 2);
            const double s_minus = pow((m_Lambdab - m_Lambdac2625), 2) - pow(s, 2);

            std::cout   << "Elena " << F12T
                        << " " << F120
                        << " " << F12P
                        << " " << F32P
                        << " " << G12T
                        << " " << G120
                        << " " << G12P
                        << " " << G32P
                        << std::endl;

            double val = pow(F12T, 2) * (pow(m_l, 2) / pow(s, 2)) * pow((m_Lambdab - m_Lambdac2625), 2) * s_plus;
            val += (pow(F120, 2) * pow((m_Lambdab + m_Lambdac2625), 2) + (pow(F12P, 2) + 3 * pow(F32P, 2)) * (pow(m_l, 2) + pow(s, 2))) * s_minus;
            val += pow(G12T, 2) * (pow(m_l, 2) / pow(s, 2)) * pow((m_Lambdab + m_Lambdac2625), 2) * s_minus;
            val += (pow(G120, 2) * pow((m_Lambdab - m_Lambdac2625), 2) + (pow(G12P, 2) + 3 * pow(G32P, 2)) * (pow(m_l, 2) + pow(s, 2))) * s_plus;

            return val;
        }

        double b_l(const double & s) const
        {
            const auto F12T = form_factors->f_time12_v(s);
            const auto F120 = form_factors->f_long12_v(s);
            const auto F12P = form_factors->f_perp12_v(s);
            const auto F32P = form_factors->f_perp32_v(s);
            const auto G12T = form_factors->f_time12_a(s);
            const auto G120 = form_factors->f_long12_a(s);
            const auto G12P = form_factors->f_perp12_a(s);
            const auto G32P = form_factors->f_perp32_a(s);
            const double s_plus = pow((m_Lambdab + m_Lambdac2625), 2) - pow(s, 2);
            const double s_minus = pow((m_Lambdab - m_Lambdac2625), 2) - pow(s, 2);

            double val = 2 * (F12T * F120 + G12T * G120) * pow(m_l, 2) / pow(s, 2);
            val *= (pow(m_Lambdab, 2) - pow(m_Lambdac2625, 2)) * sqrt(s_plus * s_minus);
            val += (-4) * pow(s, 2) * sqrt(s_plus * s_minus) * (F12P * G12P + 3 * F32P * G32P);

            return val;
        }

        double c_l(const double & s) const
        {
            const auto F12T = form_factors->f_time12_v(s);
            const auto F120 = form_factors->f_long12_v(s);
            const auto F12P = form_factors->f_perp12_v(s);
            const auto F32P = form_factors->f_perp32_v(s);
            const auto G12T = form_factors->f_time12_a(s);
            const auto G120 = form_factors->f_long12_a(s);
            const auto G12P = form_factors->f_perp12_a(s);
            const auto G32P = form_factors->f_perp32_a(s);
            const double s_plus = pow((m_Lambdab + m_Lambdac2625), 2) - pow(s, 2);
            const double s_minus = pow((m_Lambdab - m_Lambdac2625), 2) - pow(s, 2);

            double val = pow(F120, 2) * (pow(m_Lambdab, 2) + pow(m_Lambdac2625, 2));
            val += (-1) * pow(s, 2) * (pow(F12P, 2) * 3 * pow(F32P, 2));
            val += pow(G120, 2) * (pow(m_Lambdab, 2) - pow(m_Lambdac2625, 2));
            val += (-1) * pow(s, 2) * (pow(G12P, 2) * 3 * pow(G32P, 2));
            val *= (-1) * (1 - pow(m_l, 2) / pow(s, 2));

            return val;
        }

        double gamma_0(const double & s) const
        {
            const double s_plus = pow((m_Lambdab + m_Lambdac2625), 2) - pow(s, 2);
            const double s_minus = pow((m_Lambdab - m_Lambdac2625), 2) - pow(s, 2);
            const double pi = acos(-1.0);

            double val = pow(g_fermi, 2) * sqrt(s_plus * s_minus);
            val *= (1 / (384 * pow(pi, 2) * pow(m_Lambdab, 3)));
            val *= pow(1 - pow(m_l, 2) / pow(s, 2), 2);

            return val;
        }

        // normalized to V_cb = 1
        double normalized_differential_decay_width(const double & s) const
        {
            return gamma_0(s) * 2 * (a_l(s) + c_l(s) / 3);
        }

        double differential_decay_width(const double & s) const
        {
            return normalized_differential_decay_width(s) * std::norm(model->ckm_cb());
        }

        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_Lambdab / hbar;
        }
    };

    LambdabToLambdac2625LeptonNeutrino::LambdabToLambdac2625LeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<LambdabToLambdac2625LeptonNeutrino>(new Implementation<LambdabToLambdac2625LeptonNeutrino>(parameters, options, *this))
    {
    }

    LambdabToLambdac2625LeptonNeutrino::~LambdabToLambdac2625LeptonNeutrino()
    {
    }

    double
    LambdabToLambdac2625LeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    LambdabToLambdac2625LeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<LambdabToLambdac2625LeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate(f, 128, s_min, s_max);
    }

    double
    LambdabToLambdac2625LeptonNeutrino::differential_r_lambdac2625(const double & s) const
    {
        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, 0 /*_imp->parameters["mass::mu"]()*/);
            br_muons = _imp->differential_branching_ratio(s);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            br_taus = _imp->differential_branching_ratio(s);
        }

        return br_taus / br_muons;
    }

    double
    LambdabToLambdac2625LeptonNeutrino::integrated_r_lambdac2625() const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<LambdabToLambdac2625LeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, 0 /*_imp->parameters["mass::mu"]()*/);
            br_muons = integrate(f, 128, 0.02, 11.62);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            br_taus = integrate(f, 128, 3.16, 11.62);
        }

        return br_taus / br_muons;
    }

    const std::string
    LambdabToLambdac2625LeptonNeutrino::description = "\
The decay Lambda_b -> Lambda_c(2625) l nu, where l=e,mu,tau is a lepton.";

    const std::string
    LambdabToLambdac2625LeptonNeutrino::kinematics_description_s = "\
The invariant mass of the l-nubar pair in GeV^2.";

}
