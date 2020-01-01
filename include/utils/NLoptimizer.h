#ifndef NLOPTIMIZER_H
#define NLOPTIMIZER_H

#include <nlopt.hpp>
#include <vector>
#include <iostream>

struct NLoptimizer
{
    int operator()(
        const std::vector<double>& lb,
        const std::vector<double>& ub,
        nlopt::vfunc& optfunc, void* object,
        const double& tor, const int& maxIter,
        std::vector<double>& sol, double& energy)
    {
        nlopt::result result;
        nlopt::opt opt(nlopt::LD_LBFGS, static_cast<int>(lb.size()));
        opt.set_ftol_rel(tor);
        opt.set_maxeval(maxIter);
        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);
        opt.set_min_objective(optfunc, object);
        try
        {
            result = opt.optimize(sol, energy);
        }
        catch(std::exception& e)
        {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
        return result;
    }
};

#endif
