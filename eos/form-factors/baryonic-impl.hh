/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015, 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_IMPL_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>

namespace eos
{
    /* Form Factors according to [MvD2016] for J=1/2^+ -> 1/2^+ transitions */
    template <typename Process_> class MvD2016FormFactors;

    /*
     * J=1/2^+ -> J=1/2^+ transitions
     */
    struct LambdaBToLambda {
        static constexpr const char * label = "Lambda_b->Lambda";
        // initial state mass
        static constexpr double m1 = 5.61951;
        // final state mass
        static constexpr double m2 = 1.115683;
        // semileptonic kinematic endpoint
        static constexpr double tm = (m1 - m2) * (m1 - m2);
        // pair production threshold: B + K
        static constexpr double tp = (5.279 + 0.494) * (5.279 + 0.494);
        // first resonances sorted by spin/parity
        static constexpr double mR2_0m = 5.367 * 5.367;
        static constexpr double mR2_0p = 5.711 * 5.711;
        static constexpr double mR2_1m = 5.416 * 5.416;
        static constexpr double mR2_1p = 5.750 * 5.750;
    };

    template <typename Process_> class DM2016FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            UsedParameter _alpha_0_time_v, _alpha_1_time_v, _alpha_2_time_v;
            UsedParameter _alpha_0_time_a, _alpha_1_time_a, _alpha_2_time_a;

            UsedParameter _alpha_0_long_v, _alpha_1_long_v, _alpha_2_long_v;
            UsedParameter _alpha_0_long_a, _alpha_1_long_a, _alpha_2_long_a;
            UsedParameter _alpha_0_perp_v, _alpha_1_perp_v, _alpha_2_perp_v;
            UsedParameter                  _alpha_1_perp_a, _alpha_2_perp_a;

            UsedParameter _alpha_0_long_t,  _alpha_1_long_t,  _alpha_2_long_t;
            UsedParameter _alpha_0_long_t5, _alpha_1_long_t5, _alpha_2_long_t5;
            UsedParameter _alpha_0_perp_t,  _alpha_1_perp_t,  _alpha_2_perp_t;
            UsedParameter                   _alpha_1_perp_t5, _alpha_2_perp_t5;

            static constexpr double _z(const double & t, const double & tp, const double & t0)
            {
                return (std::sqrt(tp - t) - std::sqrt(tp - t0)) / (std::sqrt(tp - t) + std::sqrt(tp - t0));
            }

        public:
            DM2016FormFactors(const Parameters & p) :
                // time, V
                _alpha_0_time_v(p["Lambda_b->Lambda::a_0_time^V@DM2016"], *this),
                _alpha_1_time_v(p["Lambda_b->Lambda::a_1_time^V@DM2016"], *this),
                _alpha_2_time_v(p["Lambda_b->Lambda::a_2_time^V@DM2016"], *this),
                // time, A
                _alpha_0_time_a(p["Lambda_b->Lambda::a_0_time^A@DM2016"], *this),
                _alpha_1_time_a(p["Lambda_b->Lambda::a_1_time^A@DM2016"], *this),
                _alpha_2_time_a(p["Lambda_b->Lambda::a_2_time^A@DM2016"], *this),

                // long, V
                _alpha_0_long_v(p["Lambda_b->Lambda::a_0_long^V@DM2016"], *this),
                _alpha_1_long_v(p["Lambda_b->Lambda::a_1_long^V@DM2016"], *this),
                _alpha_2_long_v(p["Lambda_b->Lambda::a_2_long^V@DM2016"], *this),
                // long, A
                _alpha_0_long_a(p["Lambda_b->Lambda::a_0_long^A@DM2016"], *this),
                _alpha_1_long_a(p["Lambda_b->Lambda::a_1_long^A@DM2016"], *this),
                _alpha_2_long_a(p["Lambda_b->Lambda::a_2_long^A@DM2016"], *this),
                // perp, V
                _alpha_0_perp_v(p["Lambda_b->Lambda::a_0_perp^V@DM2016"], *this),
                _alpha_1_perp_v(p["Lambda_b->Lambda::a_1_perp^V@DM2016"], *this),
                _alpha_2_perp_v(p["Lambda_b->Lambda::a_2_perp^V@DM2016"], *this),
                // perp, A
                _alpha_1_perp_a(p["Lambda_b->Lambda::a_1_perp^A@DM2016"], *this),
                _alpha_2_perp_a(p["Lambda_b->Lambda::a_2_perp^A@DM2016"], *this),

                // long, T
                _alpha_0_long_t(p["Lambda_b->Lambda::a_0_long^T@DM2016"], *this),
                _alpha_1_long_t(p["Lambda_b->Lambda::a_1_long^T@DM2016"], *this),
                _alpha_2_long_t(p["Lambda_b->Lambda::a_2_long^T@DM2016"], *this),
                // long, T5
                _alpha_0_long_t5(p["Lambda_b->Lambda::a_0_long^T5@DM2016"], *this),
                _alpha_1_long_t5(p["Lambda_b->Lambda::a_1_long^T5@DM2016"], *this),
                _alpha_2_long_t5(p["Lambda_b->Lambda::a_2_long^T5@DM2016"], *this),
                // perp, T
                _alpha_0_perp_t(p["Lambda_b->Lambda::a_0_perp^T@DM2016"], *this),
                _alpha_1_perp_t(p["Lambda_b->Lambda::a_1_perp^T@DM2016"], *this),
                _alpha_2_perp_t(p["Lambda_b->Lambda::a_2_perp^T@DM2016"], *this),
                // perp, T5
                _alpha_1_perp_t5(p["Lambda_b->Lambda::a_1_perp^T5@DM2016"], *this),
                _alpha_2_perp_t5(p["Lambda_b->Lambda::a_2_perp^T5@DM2016"], *this)
            {
            }

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, unsigned)
            {
                return new DM2016FormFactors(parameters);
            }

            // vector current
            virtual double f_time_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_0p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_v() + _alpha_1_time_v() * z + _alpha_2_time_v() * z2);
            }

            virtual double f_long_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_v() + _alpha_1_long_v() * z + _alpha_2_long_v() * z2);
            }

            virtual double f_perp_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_v() + _alpha_1_perp_v() * z + _alpha_2_perp_v() * z2);
            }

            // axial vector current
            virtual double f_time_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_0m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_a() + _alpha_1_time_a() * z + _alpha_2_time_a() * z2);
            }

            virtual double f_long_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_long_a() * z + _alpha_2_long_a() * z2);
            }

            virtual double f_perp_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                // Using alpha_0_long_a instead of alpha_0_perp_a, in order to
                // fulfill relation eq. (7), [DM2016], p. 3.
                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_perp_a() * z + _alpha_2_perp_a() * z2);
            }

            // tensor current
            virtual double f_long_t(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t() + _alpha_1_long_t() * z + _alpha_2_long_t() * z2);
            }

            virtual double f_perp_t(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_t() + _alpha_1_perp_t() * z + _alpha_2_perp_t() * z2);
            }

            // axial tensor current
            virtual double f_long_t5(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_long_t5() * z + _alpha_2_long_t5() * z2);
            }

            virtual double f_perp_t5(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                // Using alpha_0_long_t5 instead of alpha_0_perp_t5, in order to
                // fulfill relation eq. (8), [DM2016], p. 3.
                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_perp_t5() * z + _alpha_2_perp_t5() * z2);
            }
    };

    /*
     * J=1/2^+ -> J=3/2^- transitions
     */
    struct LambdaBToLambdaC2625 {
        static constexpr const char * label = "Lambda_b->Lambda_c(2625)";
        // initial state mass
        static constexpr double m1 = 5.61951;
        // final state mass
        static constexpr double m2 = 2.62811;
        // semileptonic kinematic endpoint
        static constexpr double tm = (m1 - m2) * (m1 - m2);
        // pair production threshold: Lambda_b + Lambda_c(2625)
        static constexpr double tp = (m1 + m2) * (m1 + m2);
        // first resonances sorted by spin/parity
        // we use the shifts from [DLM2015], table VII.
        static constexpr double mBc = 6.2751;
        static constexpr double mR2_0m = std::pow(mBc + 0.000, 2);
        static constexpr double mR2_0p = std::pow(mBc + 0.449, 2);
        static constexpr double mR2_1m = std::pow(mBc + 0.056, 2);
        static constexpr double mR2_1p = std::pow(mBc + 0.492, 2);
    };

    template <typename Process_> class BBGIOvD2017FormFactors_constant :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:

        public:
      BBGIOvD2017FormFactors_constant(const Parameters & p)
            {
            }

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, unsigned)
            {
                return new BBGIOvD2017FormFactors_constant(parameters);
            }

            // vector current
            virtual double f_time12_v(const double & s) const
            {
                return 1.0;
            }

            virtual double f_long12_v(const double & s) const
            {
                return 1.0;
            }

            virtual double f_perp12_v(const double & s) const
            {
                return 1.0;
            }

            virtual double f_perp32_v(const double & s) const
            {
                return 0.0;
            }

            // axial vector current
            virtual double f_time12_a(const double & s) const
            {
                return 1.0;
            }

            virtual double f_long12_a(const double & s) const
            {
                return 1.0;
            }

            virtual double f_perp12_a(const double & s) const
            {
                return 1.0;
            }

            virtual double f_perp32_a(const double & s) const
            {
                return 0.0;
            }
    };


    template <typename Process_> class BBGIOvD2017FormFactors_1pole :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:

        public:
      BBGIOvD2017FormFactors_1pole(const Parameters & p)
            {
            }

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, unsigned)
            {
                return new BBGIOvD2017FormFactors_1pole(parameters);
            }

            // vector current
            virtual double f_time12_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_0m;
                return 1.0 / (1.0 - s / mR2);
            }

            virtual double f_long12_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;
                return 1.0 / (1.0 - s / mR2);
            }

            virtual double f_perp12_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;
                return 1.0 / (1.0 - s / mR2);
            }

            virtual double f_perp32_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;
                return 0.2 / (1.0 - s / mR2);
            }

            // axial vector current
            virtual double f_time12_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_0p;
                return 1.0 / (1.0 - s / mR2);
            }

            virtual double f_long12_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;
                return 1.0 / (1.0 - s / mR2);
            }

            virtual double f_perp12_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;
                return 1.0 / (1.0 - s / mR2);
            }

            virtual double f_perp32_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;
                return 0.2 / (1.0 - s / mR2);
            }
    };


    template <typename Process_> class BBGIOvD2017FormFactors_1polelambda :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:
            static constexpr double _lam(const double & mLb2, const double & mLcs2, const double & s)
            {
                return std::sqrt(eos::lambda(mLb2, mLcs2, s)) / (mLb2 - mLcs2);
            }

        public:
      BBGIOvD2017FormFactors_1polelambda(const Parameters & p)
            {
            }

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, unsigned)
            {
                return new BBGIOvD2017FormFactors_1polelambda(parameters);
            }

            // vector current
            virtual double f_time12_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_0m;
                static const double mLb2 = Process_::m1 * Process_::m1;
                static const double mLcs2 = Process_::m2 * Process_::m2;
                return 1.0 / (1.0 - s / mR2) * _lam(mLb2, mLcs2, s);
            }

            virtual double f_long12_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;
                static const double mLb2 = Process_::m1 * Process_::m1;
                static const double mLcs2 = Process_::m2 * Process_::m2;
                return 1.0 / (1.0 - s / mR2) * _lam(mLb2, mLcs2, s);
            }

            virtual double f_perp12_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;
                static const double mLb2 = Process_::m1 * Process_::m1;
                static const double mLcs2 = Process_::m2 * Process_::m2;
                return 1.0 / (1.0 - s / mR2) * _lam(mLb2, mLcs2, s);
            }

            virtual double f_perp32_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;
                return 0.2 / (1.0 - s / mR2);
            }

            // axial vector current
            virtual double f_time12_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_0p;
                static const double mLb2 = Process_::m1 * Process_::m1;
                static const double mLcs2 = Process_::m2 * Process_::m2;
                return 1.0 / (1.0 - s / mR2) * _lam(mLb2, mLcs2, s);
            }

            virtual double f_long12_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;
                static const double mLb2 = Process_::m1 * Process_::m1;
                static const double mLcs2 = Process_::m2 * Process_::m2;
                return 1.0 / (1.0 - s / mR2) * _lam(mLb2, mLcs2, s);
            }

            virtual double f_perp12_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;
                static const double mLb2 = Process_::m1 * Process_::m1;
                static const double mLcs2 = Process_::m2 * Process_::m2;
                return 1.0 / (1.0 - s / mR2) * _lam(mLb2, mLcs2, s);
            }

            virtual double f_perp32_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;
                return 0.2 / (1.0 - s / mR2);
            }
    };

    template <typename Process_> class BBGIOvD2017FormFactors_linear :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:
            static constexpr double mLb = Process_::m1;
            static constexpr double mLcs = Process_::m2;
            static constexpr double mLb2 = mLb * mLb;
            static constexpr double mLcs2 = mLcs * mLcs;
            static constexpr double s_max = (mLb - mLcs) * (mLb - mLcs);
            static constexpr double _lam(const double & mLb2, const double & mLcs2, const double & s)
            {
                return std::sqrt(eos::lambda(mLb2, mLcs2, s)) / (mLb2 - mLcs2);
            }
            UsedParameter _rho;

        public:
      BBGIOvD2017FormFactors_linear(const Parameters & p) :
        _rho(p["Lambda_b->Lambda_c(2625)::rho@BBGIOvD2017-linear"], *this)
            {

            }

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, unsigned)
            {
                return new BBGIOvD2017FormFactors_linear(parameters);
            }

            // vector current
            virtual double f_time12_v(const double & s) const
            {
                return 1.0 + _rho * (1.0 - s / s_max);
            }

            virtual double f_long12_v(const double & s) const
            {
                return 1.0 + _rho * (1.0 - s / s_max);
            }

            virtual double f_perp12_v(const double & s) const
            {
                return 1.0 + _rho * (1.0 - s / s_max);
            }

            virtual double f_perp32_v(const double & s) const
            {
                return 0.0;
            }

            // axial vector current
            virtual double f_time12_a(const double & s) const
            {
                return 1.0 + _rho * (1.0 - s / s_max);
            }

            virtual double f_long12_a(const double & s) const
            {
                return 1.0 + _rho * (1.0 - s / s_max);
            }

            virtual double f_perp12_a(const double & s) const
            {
                return 1.0 + _rho * (1.0 - s / s_max);
            }

            virtual double f_perp32_a(const double & s) const
            {
                return 0.0;
            }
    };

    template <typename Process_> class BBGIOvD2017FormFactors_HQT :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:
            static constexpr double mLb = Process_::m1;
            static constexpr double mLcs = Process_::m2;
            static constexpr double mLb2 = mLb * mLb;
            static constexpr double mLcs2 = mLcs * mLcs;
            static constexpr double lambdabar = 1.0;
            static constexpr double lambdabarprime = 1.0;
            static constexpr double s_max = (mLb - mLcs) * (mLb - mLcs);
            static constexpr double _zeta(const double & s, const double & s_max, const double & c, const double & rho)
            {
                return c + rho * (1.0 - s / s_max);
            }
            const double s_plus(const double & s)
            {
                return std::pow((mLb + mLcs), 2) - s;
            }
            const double s_minus(const double & s)
            {
                return std::pow((mLb - mLcs), 2) - s;
            }

            UsedParameter _rho, _rho3b, _c, _c3b;

        public:
      BBGIOvD2017FormFactors_HQT(const Parameters & p) :
        _rho(   p["Lambda_b->Lambda_c(2625)::rho@BBGIOvD2017-HQT"],   *this)
        _rho3b( p["Lambda_b->Lambda_c(2625)::rho3b@BBGIOvD2017-HQT"], *this)
        _c(     p["Lambda_b->Lambda_c(2625)::c@BBGIOvD2017-HQT"],     *this)
        _c3b(   p["Lambda_b->Lambda_c(2625)::c3b@BBGIOvD2017-HQT"],   *this)
            {

            }

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, unsigned)
            {
                return new BBGIOvD2017FormFactors_HQT(parameters);
            }

            // Isgur-Wise functions
            virtual double _z(const double & s) const
            {
                return _zeta(s, s_max, _c, _rho);
            }
            virtual double _z3b(const double & s) const
            {
                return _zeta(s, s_max, _c3b, _rho3b);
            }

            // vector current
            virtual double f_time12_v(const double & s) const
            {
                double res = std::sqrt(s_minus(s)) / (2.0 * std::pow(mLb * mLcs), 3.0 / 2.0);
                res *= ( ( s_plus(s) + (lambdabar/(2.0*mLb) - lambdabarprime/(2.0*mLcs)) * (mLb+mLcs) * (mLb2-mLcs2-s) / (mLb-mLcs) ) * _z(s) 
                        + _z3b(s) * (s_plus(s) + 2.0*s) / (mLb-mLcs) );
                return res;
            }

            virtual double f_long12_v(const double & s) const
            {
                double res = std::sqrt(s_plus(s)) / (2.0 * std::pow(mLb * mLcs), 3.0 / 2.0);
                res *= ( ( s_minus(s) + (lambdabar/(2.0*mLb) - lambdabarprime/(2.0*mLcs)) * (mLb-mLcs) * (mLb2-mLcs2-s) / (mLb+mLcs) ) * _z(s) 
                        + _z3b(s) * (mLb2-mLcs2-s) / (mLb+mLcs) );
                return res;
            }

            virtual double f_perp12_v(const double & s) const
            {
                double res = std::sqrt(s_plus(s)) / (2.0 * std::pow(mLb * mLcs), 3.0 / 2.0);
                res *= ( ( s_minus(s) + lambdabar * (3.0*mLb2+mLcs2-s) / (2.0*mLb) - (mLb2+3.0*mLcs2-s) * lambdabarprime / (2.0*mLcs2) ) * _z(s)
                        - mLcs * _z3b(s) );
                return res;
            }

            virtual double f_perp32_v(const double & s) const
            {
                double res = -1.0 * std::sqrt(s_plus(s)) * _z3b(s) / (2.0 * std::pow(mLb, 3.0/2.0) * std::sqrt(mLcs))
                return res;
            }

            // axial vector current
            virtual double f_time12_a(const double & s) const
            {
                double res = std::sqrt(s_plus(s)) / (2.0 * std::pow(mLb * mLcs), 3.0 / 2.0);
                res *= ( ( s_minus(s) + ((mLb-mLcs) / (mLb+mLcs)) * (lambdabar*(mLb2-mLcs2+s)/(2.0*mLb) - lambdabarprime*(mLb2-mLcs2-s)/(2.0*mLcs)) ) * _z(s) 
                        + _z3b(s) * s_minus(s) / (mLb+mLcs) );
                return res;
            }

            virtual double f_long12_a(const double & s) const
            {
                double res = std::sqrt(s_minus(s)) / (2.0 * std::pow(mLb * mLcs), 3.0 / 2.0);
                res *= ( ( s_plus(s) + ((mLb-mLcs) / (mLb+mLcs)) * (lambdabar*(mLb2-mLcs2+s)/(2.0*mLb) - lambdabarprime*(mLb2-mLcs2-s)/(2.0*mLcs)) ) * _z(s) 
                        + _z3b(s) * (mLb2-mLcs2+s) / (mLb-mLcs) );
                return res;
            }

            virtual double f_perp12_a(const double & s) const
            {
                double res = -1.0 * std::sqrt(s_plus(s)) / (2.0 * std::pow(mLb * mLcs), 3.0 / 2.0);
                res *= ( ( s_minus(s) + (lambdabar*(3.0*mLb2+mLcs2-s)/(2.0*mLb) - lambdabarprime*(mLb2+3.0*mLcs2-s)/(2.0*mLcs)) ) * _z(s) 
                        + _z3b(s) * mLcs );
                return res;
            }

            virtual double f_perp32_a(const double & s) const
            {
                double res = -1.0 * std::sqrt(s_minus(s)) * _z3b(s) / (2.0 * std::pow(mLb, 3.0/2.0) * std::sqrt(mLcs))
                return res;
            }
    };






}

#endif
