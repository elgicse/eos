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

        // normalized to V_cb = 1
        double normalized_differential_decay_width(const double & s) const
        {
            //double fp = form_factors->f_p(s);
            //double f0 = form_factors->f_0(s);
            //double lam = lambda(m_B * m_B, m_D * m_D, s);
            //double norm = power_of<2>(g_fermi())
            //    / (384.0 * power_of<3>(M_PI * m_B));
            //// make sure we return NaN if s < m_l^2
            //double sqrtv = sqrt(1.0 - m_l * m_l / s);
            //double v = sqrtv * sqrtv, v2 = v * v;
	    //
            //// correct result
            //return norm * sqrt(lam) * v2 * ((3.0 - v) * fp * fp * lam + 3.0 * f0 * f0 * (1.0 - v) * power_of<2>(m_B * m_B - m_D * m_D));
	  return 1.0;
        }

        double differential_decay_width(const double & s) const
        {
	  return 1.0; //normalized_differential_decay_width(s) * std::norm(model->ckm_cb());
        }

        double differential_branching_ratio(const double & s) const
        {
	  return 1.0; //differential_decay_width(s) * tau_B / hbar;
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
