#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <array>



/*
#include "xoshiro256ss.h"
#include "pcg_random.hpp"


int main(){
	//xoshiro256ss rng(100);
	pcg_extras::seed_seq_from<std::random_device> seed_source;

	// Make a random number engine
	pcg32_k2_fast rng(seed_source);


	uint32_t nbtrials = 500000;
	double prob_success = 1/(double)nbtrials;
	uint32_t size = 1000000 * 36 * 2;
	uint32_t nb_runs = 10;

	std::binomial_distribution<uint32_t> distrib(nbtrials, prob_success);

	std::vector<uint32_t> arr;
	arr.resize(size);

	for (size_t i = 0; i < nb_runs; ++i){
		std::binomial_distribution<uint32_t> distrib(nbtrials, prob_success);
		
		for (size_t j = 0; j < arr.size(); ++j){
			arr[j] = distrib(  rng  );	
		}
		
		std::cout<<"run "<<i+1<<"/"<<nb_runs<<" done"<<std::endl;
	}
	
	return 0;
}*/


#include <p2rng/p2rng.hpp>
#include <p2rng/trng/binomial_dist.hpp>

int main(){

	uint32_t nbtrials = 500000;
	double prob_success = 1/(double)nbtrials;
	uint32_t size = 200000000;
	uint32_t nb_runs = 10;

	const unsigned long seed{2718281828};
	std::vector<uint32_t> v(size);
	
	pcg32 rng(seed);
	//trng::binomial_dist bd;

	for (size_t i = 0; i < nb_runs; ++i){
		p2rng::generate_n
		(   std::begin(v)
		,   size
		//,   p2rng::bind(trng::uniform_int_dist(10, 100), rng)
		,   p2rng::bind(trng::binomial_dist(prob_success, nbtrials), rng)
		);
		
		std::cout<<"run "<<i+1<<"/"<<nb_runs<<" done"<<std::endl;
		
		/*for (size_t i = 0; i < size; ++i)
		{   if (0 == i % 10)
			 std::cout << '\n';
			std::cout << std::setw(3) << v[i];
		}
		std::cout << '\n' << std::endl;*/
		
		
	}

	/*for (size_t i = 0; i < 100; ++i)
	{   if (0 == i % 10)
	    std::cout << '\n';
	std::cout << std::setw(3) << v[i];
	}
	std::cout << '\n' << std::endl;*/


	return 0;
}
