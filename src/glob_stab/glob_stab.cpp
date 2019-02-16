/*
    Copyright (C) 2017 Thomas Schauss

    This file is part of glob_stab.

    glob_stab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    glob_stab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with glob_stab. If not, see <http://www.gnu.org/licenses/>.
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "bernstein_exponent.h"
#include "config.h"
#include "helpers.h"
#include "iv.h"
#include "parameters.h"
#include "taylor_model.h"

#include "optimization.h"
#include "thread_pool.hpp"

boost::program_options::options_description options()
{
    // Declare generic options.
    boost::program_options::options_description opt("Generic Options");
    opt.add_options()
        ("threads,t",    boost::program_options::value<int>()->default_value(4), "maximum number of parallel threads")
        ("parameters,p", boost::program_options::value<std::string>()->default_value("cfg/FaFa.p"), "parameter configuration (file or inline)")
        ("chEq,c",       boost::program_options::value<std::string>()->default_value("cfg/FaFa.chEq"), "characteristic equation (file or inline)")
        ("delta,d",      boost::program_options::value<double>()->default_value(0.0), "required system damping")
        ("output,o",     boost::program_options::value<std::string>()                            , "file to store results in")
        ("asy,a"                                                                                 , "store results as asy (instead of matlab)")
        ("stability,s",  boost::program_options::value<int>()->default_value(1), "stability check (0:none, 1:standard (boundary), 2:standard (right-half plane), 3:value-set only (boundary), 4:value-set only (omega+sigma))")
        ("complex,j"                                                                             , "complex equation (set s=sigma+j*omega)")
        ("no-inverse-valueset,i"                                                                 , "do not use inverse to examine unbounded omega/sigma (with s1/s2 and finalization)")
        ("finalize,f"                                                                            , "for -s1, perform finalizing step")
        ("omega_max",    boost::program_options::value<double>()->default_value(1.0e3), "Upper bound for omega in value-set only approach")
        ("omega_min",    boost::program_options::value<double>()->default_value(-0.0001), "Lower bound for omega in value-set only approach")
        ("omega_eps",    boost::program_options::value<double>()->default_value(1.0e-10), "Epsilon value for omega in standard-approach")
        ("sigma_max",    boost::program_options::value<double>()->default_value(1.0e3), "Upper bound for sigma in value-set only approach")
        ("sigma_min",    boost::program_options::value<double>()->default_value(0.0), "Lower bound for sigma in value-set only approach")
        ("sigma_eps",    boost::program_options::value<double>()->default_value(1.0e-10), "Epsilon value for sigma in standard-approach")
        ("verbosity,v",  boost::program_options::value<int>()->default_value(0),        "verbosity of output")
    ;
    return opt;
}

std::string fileOrInline(const std::string &in)
{
    std::string out;
    std::ifstream hnd(in.c_str());
    if (hnd.good()) {
        std::cout << "Loading file " << in << std::endl;
        std::getline(hnd, out, '\0');
        hnd.close();
    }
    else
        out = in;

    return out;
}

int main(int argc, char *argv[])
{
    try {
        double t1, t2, t3, t4;

        CONFIG->addOptions(options());
        CONFIG->addOptions(TaylorModel<Iv>::options());
        CONFIG->parseOptions(argc,argv);
        CONFIG->printOptions();

        Parameters::parse(fileOrInline(CONFIG_VAR(parameters,std::string)));

        double delta = CONFIG_VAR(delta,double);
        // for value-set tests
        if (CONFIG_VAR(stability,int) > 0) {
            if ( (CONFIG_VAR(stability,int) > 2) || CONFIG_ISSET(no-inverse-valueset) ) {
                Parameters::add("sigma",Iv(CONFIG_VAR(sigma_min,double),CONFIG_VAR(sigma_max,double)));
                Parameters::add("omega",Iv(CONFIG_VAR(omega_min,double),CONFIG_VAR(omega_max,double)));
            } else {
                /* we here use intervals up to 1.0001 so that normal and inverse check overlap!
                 * we use an interval for omega which includes zero, only the positive part is used for some checks (see below).
                 * unfortunately there is currently no function which allows adjusting a parameter which has been set, so this is done by splitting! */
                Parameters::add("sigma",Iv(0.0,1.0001+delta));
                Parameters::add("omega",Iv(-0.0001,1.0001));
            }
        }

        std::cout << Parameters().print();
        TaylorModel<Iv>::init(Parameters().names());

        Optimization opt(CONFIG_VAR(threads,int));
        ThreadPool tp;
        bool marginally_stable = false;

        t1 = sec();
        if (CONFIG_VAR(stability,int) == 0) {
            // Find region with all chEq > 0
            Equation ch_eq0(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq");
            opt.constrain(ch_eq0);

        } else if (CONFIG_VAR(stability,int) == 1) {
            // standard stability check

            // evaluate real root boundary (RRB) using chEq[0] for s=0
            Equation ch_eq0(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq[0]");
            if (CONFIG_ISSET(complex))
                ch_eq0.set("s",0.0);
            opt.constrain(ch_eq0);
            marginally_stable = opt.result()->hasMarginal();

            // if marginal stability is determined from RRB, evaluate chEq[1] for s=0
            if (marginally_stable) {
                Equation ch_eq1(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq[1]");
                if (CONFIG_ISSET(complex))
                    ch_eq1.set("s",0.0);
                Function fcn_eq1 = Function(ch_eq1);
                fcn_eq1.setMarginal(true); // by setting this a root at s=0 is considered as instability in Optimization::doConstrain()
                opt.constrain(fcn_eq1);
            }

            // evaluate value set for boundary
            Equation ch_eq(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq");
            ValueSet vs(ch_eq,delta);
            vs.doBoundary(true);
            if (CONFIG_ISSET(complex))
                vs.doComplex(true);

            /* it is sufficient to use positive omega as omega==0.0 is evaluated by ch_eq0,
               if system is marginally stable due to RRB we exclude omega == 0.0 and start with omega_eps!*/
            Parameters p;
            Vec<Parameters> p_vec;
            if (marginally_stable)
                p_vec = p.splitAt(p.id("omega"),CONFIG_VAR(omega_eps,double));
            else
                p_vec = p.splitAt(p.id("omega"),CONFIG_VAR(omega_min,double)<0?0:CONFIG_VAR(omega_min,double));
            Function fcn1(vs,p_vec[1]);
            tp.schedule(boost::bind(&Optimization::constrain,opt,fcn1));

            // evaluate value set for transformed characteristic function
            if (!CONFIG_ISSET(no-inverse-valueset)) {
                vs.doInverse(true);

                /* it is necessary to use positive omega for inverse check, as otherwise the exponential becomes infinite */
                p_vec = p.splitAt(p.id("omega"),0.0);
                Function fcn2(vs,p_vec[1]);
                tp.schedule(boost::bind(&Optimization::constrain,opt,fcn2));
            }

        } else if (CONFIG_VAR(stability,int) == 2) {
            // stability check where value set is mapped for complete right-half plane (RHP)

            TaylorModel<Iv>::correctTrigonometric(true); // when mapping the RHP this should be set

            // evaluate real root boundary (RRB) using chEq[0] for s=0
            Equation ch_eq0(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq[0]");
            if (CONFIG_ISSET(complex))
                ch_eq0.set("s",0.0);
            opt.constrain(ch_eq0);
            marginally_stable = opt.result()->hasMarginal();

            // if marginal stability is determined from RRB, evaluate chEq[1] for s=0
            if (marginally_stable) {
                Equation ch_eq1(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq[1]");
                if (CONFIG_ISSET(complex))
                    ch_eq1.set("s",0.0);
                Function fcn_eq1 = Function(ch_eq1);
                fcn_eq1.setMarginal(true); // by setting this a root at s=0 is considered as instability in Optimization::doConstrain()
                opt.constrain(fcn_eq1);
            }

            // evaluate value set for right-half plane (RHP)
            Equation ch_eq(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq");
            ValueSet vs(ch_eq,delta);
            if (CONFIG_ISSET(complex))
                vs.doComplex(true);

            /* in this case we use the interval containing zero for the normal check as sigma!=0 is not checked using ch_eq0,
               if system is marginally stable due to RRB we excluse omega == 0.0 and sigma == 0.0 and start with small omega and sigma!*/
            Parameters p;
            if (marginally_stable) {
                /* omega = [omega_eps, 1], sigma = [0,1] */
                Vec<Parameters> p_vec1 = p.splitAt(p.id("omega"),CONFIG_VAR(omega_eps,double));
                Function fcn1(vs,p_vec1[1]);
                tp.schedule(boost::bind(&Optimization::constrain,opt,fcn1));
                /* omega = [0, omega_eps], sigma = [sigma_eps,1] */
                Vec<Parameters> p_vec2 = p_vec1[0].splitAt(p.id("sigma"),CONFIG_VAR(sigma_eps,double));
                Function fcn2(vs,p_vec2[1]);
                tp.schedule(boost::bind(&Optimization::constrain,opt,fcn2));
            } else
                tp.schedule(boost::bind(&Optimization::constrain,opt,vs));

            // evaluate value set for transformed characteristic function
            if (!CONFIG_ISSET(no-inverse-valueset)) {
                vs.doInverse(true);

                /* it is necessary to use positive omega for inverse check, as otherwise the exponential becomes infinite */
                Vec<Parameters> p_vec = p.splitAt(p.id("omega"),0.0);
                Function fcn3(vs,p_vec[1]);
                tp.schedule(boost::bind(&Optimization::constrain,opt,fcn3));
            }

        } else if (CONFIG_VAR(stability,int) == 3) {
            // evaluate value set for boundary
            Equation ch_eq(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq");
            ValueSet vs(ch_eq,delta);
            vs.doBoundary(true);
            if (CONFIG_ISSET(complex))
                vs.doComplex(true);

            opt.constrain(vs);

        } else if (CONFIG_VAR(stability,int) == 4) {
            // evaluate value set for right-half plane (RHP)
            TaylorModel<Iv>::correctTrigonometric(true); // when mapping the RHP this should be set
            Equation ch_eq(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq");
            ValueSet vs(ch_eq,delta);
            if (CONFIG_ISSET(complex))
                vs.doComplex(true);

            opt.constrain(vs);

        }
        tp.waitDone();

        bool one_marginally_stable_finalize = false; // one of the disjoint regions is marginally stable
        t2 = sec();
        if (CONFIG_ISSET(finalize) && ( (CONFIG_VAR(stability,int) == 1) || (CONFIG_VAR(stability,int) == 3) ) ){
            // check stability of disjoint regions
            TaylorModel<Iv>::correctTrigonometric(true); // when mapping the RHP this should be set
            Vec<Parameters> final_test(opt.result()->finalTestParameters());
            for (size_t i=0; i<final_test.size(); i++) {
                // using one fixed parametrization for all following stability checks
                std::cout << "Doing final test " << i+1 << "/" << final_test.size() << " for " << final_test[i].print() << std::endl;

                // evaluate real root boundary (RRB) using chEq[0] for s=0
                Equation ch_eq0(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq[0]");
                if (CONFIG_ISSET(complex))
                    ch_eq0.set("s",0.0);
                opt.constrain(Function(ch_eq0,final_test[i]));
                bool marginally_stable_finalize = opt.result()->hasMarginal();

                // if marginal stability is determined from RRB, evaluate chEq[1] for s=0
                if (marginally_stable_finalize) {
                    one_marginally_stable_finalize = true;
                    Equation ch_eq1(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq[1]");
                    if (CONFIG_ISSET(complex))
                        ch_eq1.set("s",0.0);
                    Function fcn_eq1 = Function(ch_eq1,final_test[i]);
                    fcn_eq1.setMarginal(true); // by setting this a root at s=0 is considered as instability in Optimization::doConstrain()
                    opt.constrain(fcn_eq1);
                }

                // evaluate value set for right-half plane (RHP)
                Equation ch_eq(fileOrInline(CONFIG_VAR(chEq,std::string)),"chEq");
                ValueSet vs(ch_eq,delta);
                if (CONFIG_ISSET(complex))
                    vs.doComplex(true);

                /* in this case we use the interval containing zero for the normal check as sigma!=0 is not checked using ch_eq0,
                   if system is marginally stable due to RRB we exclude omega == 0.0, sigma == 0.0 and start with small omega or sigma!*/
                if (marginally_stable_finalize) {
                    /* omega = [omega_eps, 1], sigma = [0,1] */
                    Vec<Parameters> p_vec1 = final_test[i].splitAt(final_test[i].id("omega"),CONFIG_VAR(omega_eps,double));
                    Function fcn1(vs,p_vec1[1]);
                    tp.schedule(boost::bind(&Optimization::constrain,opt,fcn1));
                    /* omega = [-eps, omega_eps], sigma = [sigma_eps,1] */
                    Vec<Parameters> p_vec2 = p_vec1[0].splitAt(p_vec1[0].id("sigma"),CONFIG_VAR(sigma_eps,double));
                    Function fcn2(vs,p_vec2[1]);
                    tp.schedule(boost::bind(&Optimization::constrain,opt,fcn2));
                } else {
                    Function fcn1(vs,final_test[i]);
                    tp.schedule(boost::bind(&Optimization::constrain,opt,fcn1));
                }

                // evaluate value set for transformed characteristic function
                if ( (CONFIG_VAR(stability,int) == 1) && !CONFIG_ISSET(no-inverse-valueset)) {
                    vs.doInverse(true);
                    //opt.constrain(Function(vs,final_test[i]));

                    /* it is necessary to use positive omega for inverse check, as otherwise the exponential becomes infinite */
                    Vec<Parameters> p_vec = final_test[i].splitAt(final_test[i].id("omega"),0.0);
                    Function fcn3(vs,p_vec[1]);
                    tp.schedule(boost::bind(&Optimization::constrain,opt,fcn3));
                }

                tp.waitDone();
            }
        }
        t3 = sec();

        if ( CONFIG_ISSET(output) ) {
            // store results to matlab-format or Asymptote-format
            if ( CONFIG_ISSET(asy) ) {
                if (!opt.result()->storeAsy(CONFIG_VAR(output,std::string))) {
                    std::cout << std::endl <<*(opt.result()) << std::endl;
                }
            } else {
                opt.result()->storeMatlab(CONFIG_VAR(output,std::string));
            }
        } else {
            // output results to stdout
            std::cout << std::endl <<*(opt.result()) << std::endl;
        }
        t4 = sec();

        // Display statistics and other usefull information

        std::cout << std::endl << opt << std::endl;

        if (marginally_stable) {
            std::cout << std::endl;
            std::cout << "In 'normal' evaluation:" << std::endl;
            if ( (CONFIG_VAR(stability,int) == 1) || (CONFIG_VAR(stability,int) == 3) ) {
                std::cout << "System was marginally stable => excluded omega == 0.0 from value-set check!" << std::endl;
                std::cout << "Excluded omega < " << CONFIG_VAR(omega_eps,double) << std::endl << std::endl;
            } else {
                std::cout << "System was marginally stable => excluded omega == 0.0 and sigma == 0.0 from value-set check!" << std::endl;
                std::cout << "Excluded omega < " << CONFIG_VAR(omega_eps,double) << " and sigma < " << CONFIG_VAR(sigma_eps,double) << std::endl << std::endl;
            }
        }

        if (one_marginally_stable_finalize) {
            std::cout << std::endl;
            std::cout << "In 'finalization':" << std::endl;
            std::cout << "System was marginally stable => excluded omega == 0.0 and sigma == 0.0 from value-set check!" << std::endl;
            std::cout << "Excluded omega < " << CONFIG_VAR(omega_eps,double) << " and sigma < " << CONFIG_VAR(sigma_eps,double) << std::endl << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Bounding time:    " << t2-t1 << " s" << std::endl;
        std::cout << "Finalizing time:  " << t3-t2 << " s" << std::endl;
        std::cout << "Store/Print time: " << t4-t3 << " s" << std::endl << std::endl;
        std::cout << "Total time:       " << t4-t1 << " s" << std::endl << std::endl;
    }
    catch (ConfigException &err) {
        std::cerr << "Config error: " << err.what() << std::endl;
        CONFIG->printDescription();
        return 1;
    }
    catch (std::exception &err) {
        std::cerr << "std::exception: " << err.what() << std::endl;
        return 1;
    }
    catch (const char* err) {
        std::cerr << "Error: " << err << std::endl;
        return 1;
    }

    Config::clear();
    BernsteinExponent::clear();

    return 0;
}


