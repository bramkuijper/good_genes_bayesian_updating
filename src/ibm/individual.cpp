#include <cmath>
#include "individual.hpp"
#include "parameters.hpp"

// standard constructor
Individual::Individual(Parameters const &params) :
    t{params.init_t, params.init_t},
    p{params.init_p, params.init_p},
    prior_mean_x{params.init_prior_mean_x,params.init_prior_mean_x},
    prior_sigma_x{params.init_prior_sigma_x,params.init_prior_sigma_x},
    v{params.init_v, params.init_v},
    x{params.init_t * std::exp(-std::fabs(params.v_opt - params.init_v))}
{
}

Individual::Individual(Individual const &other) :
    t{other.t[0], other.t[1]},
    p{other.p[0], other.p[1]},
    prior_mean_x{other.prior_mean_x[0],other.prior_mean_x[1]},
    prior_sigma_x{other.prior_sigma_x[0],other.prior_sigma_x[1]},
    v{other.v[0], other.v[1]},
    x{other.x}
{}


Individual::Individual(
                Individual const &mother,
                Individual const &father,
                std::mt19937 &rng_r,
                Parameters const &params) :
    x{0.0}
{
    std::bernoulli_distribution segregator{0.5};
    std::normal_distribution standard_normal{};
    std::uniform_real_distribution<double> uniform{0.0,1.0};

    p[0] = mother.p[segregator(rng_r)];
    p[1] = father.p[segregator(rng_r)];

    t[0] = mother.t[segregator(rng_r)];
    t[1] = father.t[segregator(rng_r)];

    v[0] = mother.v[segregator(rng_r)];
    v[1] = father.v[segregator(rng_r)];

    prior_mean_x[0] = mother.prior_mean_x[segregator(rng_r)];
    prior_mean_x[1] = father.prior_mean_x[segregator(rng_r)];

    prior_sigma_x[0] = mother.prior_sigma_x[segregator(rng_r)];
    prior_sigma_x[1] = father.prior_sigma_x[segregator(rng_r)];


    double incr;

    for (unsigned allele_idx{0}; allele_idx < 2; ++allele_idx)
    {
        if (uniform(rng_r) < params.mu_p)
        {
            incr = uniform(rng_r) * params.max_mut_p;

            p[allele_idx] = p[allele_idx] + (uniform(rng_r) < 0.5 ? -incr : incr);
        }

        if (uniform(rng_r) < params.mu_t)
        {
            incr = uniform(rng_r) * params.max_mut_t;
            t[allele_idx] = t[allele_idx] + (uniform(rng_r) < 0.5 ? -incr : incr);
        }
        
        if (uniform(rng_r) < params.mu_v)
        {
            incr = uniform(rng_r) * params.max_mut_v;
            v[allele_idx] = v[allele_idx] + (uniform(rng_r) < params.biasv ? -incr : incr);
        }

        if (uniform(rng_r) < params.mu_prior_mean_x)
        {
            incr = standard_normal(rng_r) * params.mu_stddev;
            prior_mean_x[allele_idx] = prior_mean_x[allele_idx] + incr;
        }

        if (uniform(rng_r) < params.mu_prior_sigma_x)
        {
            incr = standard_normal(rng_r) * params.mu_stddev;
            prior_sigma_x[allele_idx] = prior_sigma_x[allele_idx] + incr;
            
            if (prior_sigma_x[allele_idx] < 0)
            {
                prior_sigma_x[allele_idx] = 0.0;
            }
        }
    } // end for allele_idx
} // end Individual()

void Individual::operator=(Individual const &other)
{
    for (unsigned allele_idx{0}; allele_idx < 2; ++allele_idx)
    {
        p[allele_idx] = other.p[allele_idx];
        t[allele_idx] = other.t[allele_idx];
        v[allele_idx] = other.v[allele_idx];
        prior_mean_x[allele_idx] = other.prior_mean_x[allele_idx];
        prior_sigma_x[allele_idx] = other.prior_sigma_x[allele_idx];
    }

    x = other.x;
}//end operator=
