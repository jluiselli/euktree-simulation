#pragma once

#include <vector>
#include <limits.h>

//#include "rando_xo.hpp"
#include <p2rng/p2rng.hpp>
#include <p2rng/trng/binomial_dist.hpp>


/**
Generate random numbers using the super fast parallel p2rng library (which must be installed, see https://github.com/arminms/p2rng).
Note an issue with this library: the rng is reset after any called to generate_...
This means that unless the seed is changed, two calls to generate_... will produce the same results.
They say it's for reproducibility, but also offer no solution.  So, after each call to generate_..., the seed is updated, 
by generating another random seed.  For that, we use rng_xo.
**/
namespace rando_mp{

	uint64_t seed;	
	std::random_device dev;
	std::mt19937 rngmt;
	std::uniform_int_distribution<std::mt19937::result_type> distro;
	
	void renew_seed(){
		//seed = rando_xo::next_int();
		/*pcg32 rng(seed);
		int max = INT_MAX/2;
		
		auto t = trng::uniform_int_dist(1, max);
		seed = t(rng);*/
		
		

		seed = distro(rngmt);
		
	}
	
	
	//TODO: this uses the basic but fast rng algorithm - which is known to be good enough but does not pass every test
	void generate_random_integers(std::vector<uint32_t> &vec_to_fill, uint32_t min, uint32_t max, size_t size) {
		
		//old slow attempt
		/*std::random_device dev;
		std::mt19937 rngmt(dev());
		std::uniform_int_distribution<std::mt19937::result_type> dist6(min,max-1);

		for (int i = 0; i < vec_to_fill.size(); ++i){
		vec_to_fill[i] = dist6(rngmt);
		}
		return;*/
				
		pcg32 rng(seed);
		p2rng::generate_n
		(   
			std::begin(vec_to_fill), size, p2rng::bind(trng::uniform_int_dist(min, max), rng)
		);
		
		
		renew_seed();
	}
	
	void init(uint64_t initial_seed){
		seed = initial_seed;
		
		rngmt = std::mt19937(initial_seed);
		distro = std::uniform_int_distribution<std::mt19937::result_type>(1,INT_MAX/2);


	}
	
	
	
	
	//TODO: this uses the basic but fast rng algorithm - which is known to be good enough but does not pass every test
	void generate_binomials(std::vector<uint32_t> &vec_to_fill, uint32_t nbtrials, double prob_success, size_t size) {
		pcg32 rng(seed);
		p2rng::generate_n
		(   
			std::begin(vec_to_fill), size, p2rng::bind(trng::binomial_dist(prob_success, nbtrials), rng)
		);			
		renew_seed();
	}
	
	
	//TODO: this uses the basic but fast rng algorithm - which is known to be good enough but does not pass every test
	//template <class It>
	//void generate_binomials(It &iterator, uint32_t nbtrials, double prob_success, size_t size) {
	
	
    
	/*	pcg32 rng(seed);
		p2rng::generate_n
		(   
			iterator, size, p2rng::bind(trng::binomial_dist(prob_success, nbtrials), rng)
		);			
		renew_seed();
	}*/
	/*std::random_device rd;
	std::minstd_rand gen(rd());
	//pcg32 rngx(seed);
	// perform 4 trials, each succeeds 1 in 2 times
	std::binomial_distribution<> d(nbtrials, prob_success);

	for (int i = 0; i < size; ++i){
		*iterator = d(gen);
		iterator++;
	}
	return;*/
		
}





