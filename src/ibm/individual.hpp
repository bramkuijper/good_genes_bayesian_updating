#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        double t[2]{0.0,0.0};
        double p[2]{0.0,0.0};
        double prior_mean_x[2]{0.0,0.0};
        double prior_sigma_x[2]{0.0,0.0};
        double v[2]{0.0,0.0};

        unsigned mate;

        // realized ornament
        double x{0.0};

        // male survival prob

        // standard constructor
        Individual(Parameters const &param);

        // copy ctor
        Individual(Individual const &other);

        // birth constructor
        Individual(Individual const &mother,
                Individual const&father,
                std::mt19937 &rng_r,
                Parameters const &params);

        void operator=(Individual const &other);
};

#endif 
