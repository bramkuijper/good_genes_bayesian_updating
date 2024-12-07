#include <cassert>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include "good_genes.hpp"
#include "parameters.hpp"

// the constructor, which initializes 
// and then runs the whole simulation
GoodGenes::GoodGenes(Parameters const &params) :
    par{params}, // initialize the parameter object
    data_file{par.file_name}, // start the output file
    uniform{0.0,1.0}, // initialize a standard uniform distribution
    rd{}, // rd() generates the seed fed to the random number generator
    seed{rd()}, // save the seed
    rng_r{seed}, // initialize the random number generator
    males(par.n/2,Individual(par)), // initialize males
    females(par.n/2,Individual(par)) // initialize females
{
    // add headers to the data file
    write_data_headers();

    // set phenotypes for the first generation
    phenotypes();

    for (time_step = 0; 
            time_step <= par.max_num_gen; ++time_step)
    {
        survival();
        reproduction();
        phenotypes();

        if (time_step % par.numoutgen == 0)
        {
            write_data();
        }
    }

    write_parameters();

} // end GoodGenes() constructor


void GoodGenes::phenotypes()
{
    double t,v;

    for (auto male_iter{males.begin()}; male_iter != males.end(); ++male_iter)
    {
        t = 0.5 * (male_iter->t[0] + male_iter->t[1]);
        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);

        male_iter->x = t * std::exp(-std::fabs(par.v_opt - v));
    }
}

void GoodGenes::survival()
{
    // aux variables to store trait values
    double p,v,surv;

    unsigned nm = males.size();
    unsigned nf = females.size();

    // initialize survival probabilities
    mean_p_survive_f = 0.0;
    mean_p_survive_m = 0.0;
    
    for (auto female_iter{females.begin()}; female_iter != females.end(); )
    {
        p = 0.5 * (female_iter->p[0] + female_iter->p[1]);
        v = 0.5 * (female_iter->v[0] + female_iter->v[1]);

        surv = std::exp(-std::fabs(par.v_opt - v)) ;

        // update survival prob
        mean_p_survive_f += surv;

        // individual dies
        if (uniform(rng_r) > surv)
        {
            std::swap(*female_iter, females.back()); // swap this female with final elmt
            females.pop_back(); // remove final elmt
                                // this is way faster than std::erase() 
        }
        else
        {
            ++female_iter;
        }
    } // end for females

    double x;

    for (auto male_iter{males.begin()}; male_iter != males.end(); )
    {
        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);

        x = male_iter->x;

        surv = std::exp(-par.c * x * x - std::fabs(par.v_opt - v)) ;
        
        mean_p_survive_m += surv;

        // individual dies
        if (uniform(rng_r) > surv)
        {
            std::swap(*male_iter, males.back()); // swap this male with final elmt
            males.pop_back(); // remove final elmt
                                // this is way faster than std::erase() 
        }
        else
        {
            ++male_iter;
        }
    }

    mean_p_survive_f /= nf;
    mean_p_survive_m /= nm;
} // end survival


void GoodGenes::reproduction()
{
    long unsigned female_idx,male_idx;

    std::vector<Individual> nextgen{};

    std::uniform_int_distribution<unsigned> female_sampler(0, females.size() - 1);

    for (unsigned newborn_idx{0}; newborn_idx < par.n; ++newborn_idx)
    {
        female_idx = female_sampler(rng_r);

        assert(female_idx >= 0);
        assert(female_idx < females.size());

        male_idx = choose(females[female_idx]);
        
        assert(male_idx >= 0);
        assert(male_idx < males.size());

        Individual Kid(
                females[female_idx],
                males[male_idx],
                rng_r,
                par);

        nextgen.push_back(Kid);
    }

    assert(nextgen.size() == par.n);

    females.clear();
    males.clear();

    for (unsigned newborn_idx{0}; newborn_idx < par.n; ++newborn_idx)
    {
        if (uniform(rng_r) < 0.5)
        {
            males.push_back(nextgen[newborn_idx]);
        } else
        {
            females.push_back(nextgen[newborn_idx]);
        }
    }
} // end reproduction

void GoodGenes::write_parameters()
{
    data_file 
        << std::endl 
        << std::endl
        << "seed;" << seed << ";" << std::endl
        << "n;" << par.n << ";" << std::endl
        << "mu_p;" << par.mu_p << ";" << std::endl
        << "mu_t;" << par.mu_t << ";" << std::endl
        << "mu_v;" << par.mu_v << ";" << std::endl
        << "mu_prior_mean_x;" << par.mu_prior_mean_x << ";" << std::endl
        << "mu_prior_sigma_x;" << par.mu_prior_sigma_x << ";" << std::endl
        << "max_mut_p;" << par.max_mut_p << ";" << std::endl
        << "max_mut_t;" << par.max_mut_t << ";" << std::endl
        << "max_mut_v;" << par.max_mut_v << ";" << std::endl
        << "mu_stddev;" << par.mu_stddev << ";" << std::endl
        << "biasv;" << par.biasv << ";" << std::endl
        << "a;" << par.a << ";" << std::endl
        << "b;" << par.b << ";" << std::endl
        << "c;" << par.c << ";" << std::endl
        << "choice_sample_size;" << par.choice_sample_size << ";" << std::endl
        << "init_t;" << par.init_t << ";" << std::endl
        << "init_p;" << par.init_p << ";" << std::endl
        << "init_v;" << par.init_v << ";" << std::endl
        << "v_opt;" << par.v_opt << ";" << std::endl
        << "max_num_gen;" << par.max_num_gen << ";" << std::endl
        << "numoutgen;" << par.numoutgen << ";" << std::endl;
}

void GoodGenes::write_data()
{
    double meanp{0.0};
    double ssp{0.0};
    double meant{0.0};
    double sst{0.0};
    double meanv{0.0};
    double ssv{0.0};
    double meanx{0.0};
    double ssx{0.0};
    
    //  covariances
    double stv{0.0};
    double stp{0.0};
    double spv{0.0};


    // keep track of population sizes
    unsigned long nf{females.size()};
    unsigned long nm{males.size()};

    // aux variables to store trait values
    double p,t,v,x;

    for (auto female_iter{females.begin()};
            female_iter != females.end();
            ++female_iter)
    {
        p = 0.5 * (female_iter->p[0] + female_iter->p[1]);
        meanp += p;
        ssp += p*p;
        
        t = 0.5 * (female_iter->t[0] + female_iter->t[1]);
        meant += t;
        sst += t*t;

        v = 0.5 * (female_iter->v[0] + female_iter->v[1]);
        meanv += v;
        ssv += v*v;

        stp += t * p;
        spv += p * v;
        stv += t * v;
    }
    
    for (auto male_iter{males.begin()};
            male_iter != males.end();
            ++male_iter)
    {
        p = 0.5 * (male_iter->p[0] + male_iter->p[1]);
        meanp += p;
        ssp += p*p;
        
        t = 0.5 * (male_iter->t[0] + male_iter->t[1]);
        meant += t;
        sst += t*t;

        v = 0.5 * (male_iter->v[0] + male_iter->v[1]);
        meanv += v;
        ssv += v*v;
        
        x = male_iter->x;
        meanx += x;
        ssx += x*x;
        
        stp += t * p;
        spv += p * v;
        stv += t * v;
    }


    meanp /= (nf + nm);
    meant /= (nf + nm);
    meanv /= (nf + nm);
    meanx /= nm;
    
    double varp = ssp / (nf + nm) - meanp * meanp;
    double vart = sst / (nf + nm) - meant * meant;
    double varv = ssv / (nf + nm) - meanv * meanv;
    double varx = ssx / nm - meanx * meanx;

    double covtp = stp / (nf + nm) - meant * meanp;
    double covtv = stv / (nf + nm) - meant * meanv;
    double covpv = spv / (nf + nm) - meanp * meanv;

    data_file << time_step << ";"
        << meanp << ";"
        << meant << ";"
        << meanv << ";"
        << meanx << ";"
        << varp << ";"
        << vart << ";"
        << varv << ";"
        << varx << ";"
        << covtp << ";"
        << covtv << ";"
        << covpv << ";"
        << mean_p_survive_f << ";"
        << mean_p_survive_m << ";"
        << nf << ";"
        << nm << ";" << std::endl;

} // write_data()

void GoodGenes::write_data_headers()
{
    data_file << "generation;meanp;meant;meanv;meanx;varp;vart;varv;varx;covtp;covtv;covpv;surv_f;surv_m;nf;nm;" << std::endl;
}

// choose surviving male according to its ornament 
void GoodGenes::choose(Individual const &female, double &cost, unsigned chosen_male)
{
    std::vector <double> male_fitness{};
    std::vector <unsigned> male_idxs{};

    // distribution to sample males from
    std::uniform_int_distribution<unsigned> 
        male_sampler(0, males.size() - 1);

    unsigned sampled_male_idx;
    
    double fitness;

    double p = 0.5 * (female.p[0] + female.p[1]);

    double x;

    unsigned mate;

    // set cost to 0
    cost = 0.0;
    
    // develop the prior mean
    double prior_mean = female.prior_mean_x[0] + female.prior_mean_x[1];
    double prior_sigma = female.prior_sigma_x[0] + female.prior_sigma_x[1];
    double prior_sigma_sq_inv = 1.0 / (prior_sigma * prior_sigma);

    double posterior_mean{0.0};
    double posterior_sigma{0.0};
    double posterior_variance{0.0};
    double posterior_var_inv{0.0};

    for (unsigned inspected_male_idx{0}; 
            inspected_male_idx < par.choice_sample_size; 
            ++inspected_male_idx)
    {
        sampled_male_idx = male_sampler(rng_r);

        handicap = males[sampled_male_idx].x;

        // sample perceived ornament distribution from posterior
        estimate = prior_mean + prior_sigma * standard_normal(rng_r);

        if (estimate >= p)
        {
            chosen_male = sampled_male_idx;
        }
        else
        {
            // ok assess mate at a cost c of each assessment the female
            ++n_mates_assessed; 

            posterior_var_inv = 1.0 / posterior_sigma * posterior_sigma;

            // update estimate based on mean and variance
            // this is eq. (11) of Luttbeg (1996)
            posterior_mean = (
                    handicap * prior_sigma_sq_inv + 
                        posterior_mean * 1.0 / posterior_var_inv
                    )/ (posterior_var_inv + prior_sigma_sq_inv);
            
            posterior_variance = 1.0 / prior_sigma_sq_inv;

            mean_ornaments += handicap;
            ss_ornaments += handicap * handicap;
            
            if (n_mates_assessed > 1)
            {
                // this is rho in eq. (12)
                double ornament_variance_sampled_males = ss_ornaments / n_mates_assessed - 
                    mean_ornaments / n_mates_assessed;

                if (ornament_variance_sampled_males <= 0)
                {
                    ornament_variance_sampled_males = 1.0 / prior_sigma_sq_inv;
                }

                // yielding rho_new in eq. (12)
                posterior_variance = 1.0 / (prior_sigma_sq_inv + 
                        1.0/ornament_variance_sampled_males);

            }

            // now sample estimate after assessing male
            estimate = posterior_mean + 
                sqrt(posterior_variance) * standard_normal(rng_r);

            if (estimate >= preference)
            {
                return(sampled_male_idx);
            }
        } // end else (i.e., estimate < preference)
    } //end for inspected male_idx

    // now make distribution of the fitnesses to choose from
    std::discrete_distribution<unsigned> choose_male(male_fitness.begin(),
            male_fitness.end());

    unsigned the_chosen_male = male_idxs[choose_male(rng_r)];

    return(the_chosen_male);
} // choose()



