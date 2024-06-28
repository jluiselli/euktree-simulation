#pragma once

#include <vector>

//#include "xoshiro256ss.h"

namespace rando_xo{

	//std::random_device rng_device;
	//std::mt19937 rng_generator;
	//xoshiro256ss rng_generator(100);



	static inline uint32_t rotl(const uint32_t x, int k) {
		return (x << k) | (x >> (32 - k));
	}
	static uint32_t s[4];


	uint32_t next_int(void) {

		const uint32_t result = rotl(s[1] * 5, 7) * 9;

		const uint32_t t = s[1] << 9;

		s[2] ^= s[0];
		s[3] ^= s[1];
		s[1] ^= s[2];
		s[0] ^= s[3];

		s[2] ^= t;

		s[3] = rotl(s[3], 11);

		return result;
	}
	
	
	
	//TODO: this uses the basic but fast rng algorithm - which is known to be good enough but does not pass every test
	void generate_random_integers(std::vector<uint32_t> &vec_to_fill, uint32_t min, uint32_t max, size_t size) {
	    //uniform_uint32_t_distribution<uint32_t> dist(low, high - 1);
		vec_to_fill.resize(size);
		generate(vec_to_fill.begin(), vec_to_fill.end(), [&]() { return next_int() % (max - min) + min; });
	}
	
	void init(uint32_t intseed1, uint32_t intseed2, uint32_t intseed3, uint32_t intseed4){
		s[0] = intseed1; s[1] = intseed2; s[2] = intseed3; s[3] = intseed4;
		
		//rng_generator.seed(binom_seed);
	}
	
	
	/*
	//TODO: this uses the basic but fast rng algorithm - which is known to be good enough but does not pass every test
	void generate_binomials(std::vector<uint32_t> &vec_to_fill, uint32_t nbtrials, double prob_success, size_t size) {

		std::binomial_distribution<uint32_t> distrib(nbtrials, prob_success);
		vec_to_fill.resize(size);
		for (size_t i = 0; i < vec_to_fill.size(); ++i){
			vec_to_fill[i] = distrib(  rng_generator  );
			
		}	
	}
	*/
	
	
	
}
