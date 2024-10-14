#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <algorithm>
#include <array>
#include <set>
#include <map>

#include <stdint.h>
#include <getopt.h>
  
#include "rando_xo.hpp" 
#include "rando_mp.hpp"

#include "pdqsort.h"


#include <omp.h>



#include <future>
#include <string>
#include <mutex>



//#ifdef USETBB
//	#include <execution>
//#endif

#include "ewah/ewah.h"
using namespace ewah;
typedef EWAHBoolArray<uint32_t> bitmap;

using namespace std;


/**
Note: we use uint32_t all over the place - make sure we never overflow.
**/

//used to compile with 
//g++ -O3 main.cpp -march=native -I/usr/include/tbb -lpthread -ltbb -DUSETBB









//TODO: storing all parameters globally is an anti-pattern, no one should ever do that
namespace config{
	uint32_t nbchr = 36;
	uint32_t chrlen = 500000;
	uint32_t nb_gen = 1000;
	uint32_t pop_size = 2000;
	double recomb_rate = 1/((double)chrlen);	
	uint32_t seed = 432498743;
	uint32_t init_nb_indiv = 0;
}




void pause(){
	cout<<"Press any key to continue.";
	cin.get();
	cout<<endl;
}








struct Segment{
        uint32_t a;
        uint32_t b;
        uint32_t individ;
        uint32_t chrno;
        uint32_t chrindex;
        
        
        Segment() : Segment(0, 0, 0, 0, 0){
        
        }
        
        Segment(uint32_t individ, uint32_t chrno, uint32_t chrindex, uint32_t a, uint32_t b) 
			: individ(individ), chrno(chrno), chrindex(chrindex), a(a), b(b){
        	
        }
        
        
        
        inline void setAttr(uint32_t individ, uint32_t chrno, uint32_t chrindex, uint32_t a, uint32_t b){
        	this->individ = individ;
        	this->chrno = chrno;
        	this->chrindex = chrindex;
        	this->a = a;
        	this->b = b;
        }
        
        
        
        
        void print(){
        	cout<<"individ="<<individ<<" chrno="<<chrno<<" chrindex="<<chrindex<<" a="<<a<<" b="<<b<<endl;
        }
        
        
        /* 
        Overloading < operator for sorting purposes
        It sorts by individ, then a, then b.  Sorting by chrno and chrindex is not needed by SegmentList at the moment.
        */
	bool operator<(const Segment& s) const
	{ 
		if (individ < s.individ) return true;
		if (individ > s.individ) return false;
		//if (chrno < s.chrno) return true;
		//if (chrno > s.chrno) return false;
		//if (chrindex < s.chrindex) return true;
		//if (chrindex > s.chrindex) return false;
		if (a < s.a) return true;
		if (a > s.a) return false;
		if (b < s.b) return true;
		return false;
		
		//previous way of sorting, slightly slower
		//return std::tie( individ, chrno, chrindex, a, b ) < std::tie( s.individ, s.chrno, s.chrindex, s.a, s.b ); 
	} 
};





/**
A SegmentList stores the list of segments of a whole population.  The segments are splits into one list per chr (two homologous chr are considered distinct).
That is, for each of the 2*nbchr chromosomes, there is a separate list to store the segments of ALL indivs that belong to a chr.
This accelerates sorting.  It would make more conceptual sense to store two lists per individual, but that would scatter the segments too much and reserve too much memory.
As a result, most functions require the chrno and chrindex, which tell which of the lists to use.
As a note, the memory taken in each list is NEVER cleared, even if the resizing functions are called (resizing a list only updates a size variable, but does not decrease memory).
So, even if segments[c].size() if large, only sizes[c] of its entries are actually used.  This allows reusing memory when adding, as a segment is only created once.
**/
class SegmentList{
public:
	vector<uint32_t> sizes;
	vector<vector<Segment>> segments;

	SegmentList(){
		init();
	}


	void init(){
		sizes.resize(config::nbchr * 2);
		segments.resize(config::nbchr * 2);
		for (uint32_t c = 0;  c < config::nbchr; ++c){
			sizes[2*c] = 0;
			segments[2*c].resize(config::pop_size * 2);	//initial memory, may get larger later
			
			sizes[2*c + 1] = 0;
			segments[2*c + 1].resize(config::pop_size * 2);
		}

	}
	

	/**
	Adds a segment to the list of chrno->chrindex.  Will automatically reserve more memory if needed.
	**/	
	inline void add(uint32_t individ, uint32_t chrno, uint32_t chrindex, uint32_t a, uint32_t b){

		if (sizes[2*chrno+chrindex] + 1 >= segments[2*chrno+chrindex].size()){
			//cout<<"resizing chr "<<chrno<<"."<<chrindex<<" to "<<segments[2*chrno+chrindex].size()*2<<endl;
			segments[2*chrno+chrindex].resize( segments[2*chrno+chrindex].size() * 2 );

		}
		
		set(sizes[2*chrno+chrindex], individ, chrno, chrindex, a, b);
		
		
		sizes[2*chrno+chrindex]++;
	}
	

	/**
	Sets the parameters of the i-th segment of the list of chrno->chrindex.  Make sure that i is less than get_size(chrno, chrindex).
	**/
	void set(size_t i, uint32_t individ, uint32_t chrno, uint32_t chrindex, uint32_t a, uint32_t b){
		segments[2*chrno+chrindex][i].setAttr(individ, chrno, chrindex, a, b);
	}
	
	
	/**
	Returns a REFERENCE to the i-th segment of the list of chrno->chrindex
	**/
	Segment& get(size_t i, uint32_t chrno, uint32_t chrindex){
		return segments[2*chrno+chrindex][i];
	}
	
	
	/**
	Returns the number of segments added to the list of chrno->chrindex
	**/
	size_t get_size(uint32_t chrno, uint32_t chrindex){
		return sizes[2*chrno+chrindex];
	}
	
	
	/**
	Sets the number of segments in the list of chrno->chrindex.  Reserves more memory if needed.
	**/
	void resize(size_t newsize, uint32_t chrno, uint32_t chrindex){
		if (newsize >= segments[2*chrno+chrindex].size()){
			segments[2*chrno+chrindex].resize(newsize);
		}
		sizes[2*chrno+chrindex] = newsize;
	}
	
	
	/**
	Sets the number of segments of EACH list to 0.  Does not reduce memory.
	**/
	void reset_sizes(){
		for (uint32_t c = 0; c < config::nbchr * 2; ++c){
			sizes[c] = 0;
		}
	}
	
	
	
	/**
	Removes the segments of list chrno->chrindex at indices in indices_to_delete.  The latter MUST be sorted.
	**/
	void remove_indices(uint32_t chrno, uint32_t chrindex, vector<size_t> &indices_to_delete){
		
		/*std::set<size_t> indset(indices_to_delete.begin(), indices_to_delete.end());
		vector<Segment> vec;
		
		for (int i = 0; i < sizes[2*chrno+chrindex]; ++i){
			if (!indset.count(i)){
				vec.push_back( segments[2*chrno+chrindex][i] );
			}
			
		}
		segments[2*chrno+chrindex] = vec;
		sizes[2*chrno+chrindex] -= indices_to_delete.size();*/
		
		
		if (indices_to_delete.empty())
			return;

		uint32_t nb_dels = 0;
		
		size_t cur_delindex = 0;
		
		for (size_t i = 0; i < sizes[2*chrno+chrindex] - indices_to_delete.size(); ++i){
			
			while (cur_delindex < indices_to_delete.size() && i + nb_dels == indices_to_delete[cur_delindex]){
				nb_dels++;		
				cur_delindex++;			
			}
			
			if (nb_dels > 0){
				Segment& seg = this->get(i + nb_dels, chrno, chrindex);
				this->set(i, seg.individ, seg.chrno, seg.chrindex, seg.a, seg.b);
			}
		}
		sizes[2*chrno+chrindex] -= indices_to_delete.size();
	}
	
	
	
	/**
	Sorts every lists, using Segment::operator< for the sorting (which sorts by individ, a, b)
	**/
	void sort(){
		
		//this is stuff I tried
		/*vector<uint32_t> v;
		for (uint32_t c = 0; c < config::nbchr*2; ++c){
			v.push_back(c);
		}
		std::for_each(
		    std::execution::par,
		    v.begin(),
		    v.end(),
		    [this](auto&& c)
		    {
			std::sort(segments[c].begin(), segments[c].begin() + sizes[c]);
		});*/
		
		//__gnu_parallel::sort(segments.begin(), segments.begin() + size);
		
		
		//this is quite fast, but unfortunately std::execution::par uses tbb implementation, which has a known memory leak
		/*#ifdef USETBB
			std::sort(std::execution::par, segments.begin(), segments.begin() + size);
		#else
			std::sort(segments.begin(), segments.begin() + size);
		#endif*/
		
		
		#pragma omp parallel for 
		for (uint32_t c = 0; c < config::nbchr * 2; ++c){
			pdqsort(segments[c].begin(), segments[c].begin() + sizes[c]);
			//std::sort(segments[c].begin(), segments[c].begin() + sizes[c]);
		}
	}
	
	
	/**
	Returns the sum of sizes of each segment list.
	**/
	size_t get_total_nb_segments(){
		size_t total = 0;
		for (uint32_t c = 0; c < config::nbchr * 2; ++c){
			total += sizes[c];
		}
		return total;
	}
	
	
	size_t get_total_segment_size(){
		size_t total_size = 0;
		for (uint32_t chrno = 0; chrno < config::nbchr; ++chrno){
			for (uint32_t chrindex = 0; chrindex < 2; ++chrindex){
				for (size_t i = 0; i < get_size(chrno, chrindex); ++i){
					Segment &seg = get(i, chrno, chrindex);
					total_size += seg.b - seg.a + 1;
				}
			}
		}
		return total_size;
	}
	
};









/**
A population simply stores the parameters of its individuals in a set of vectors.  There is no Individual class, as their parameters are stored in 
the population vectors.
**/
class Population{
public:
	//The parents of indiv i are par_ids[2*i], par_ids[2*i+1]
	vector<uint32_t> par_ids; 	
	
	//For indiv i, the chrindex of the parent that gave chrno->chrindex is chr_choices[i*nbchr*2+2*chrno+chrindex].  0 or 1.
	vector<uint32_t> chr_choices;
	
	//members[i] = 1 => is a genealogical ancestor, and 0 otherwise
	vector<uint32_t> members;	
	
	//For indiv i, the nb of breakpoints in chrno->chrindex is chr_nb_breakpoints[i*nbchr*2+2*chrno+chrindex].  0 or 1.
	vector<uint32_t> chr_nb_breakpoints;
	
	//Follow the sets of children in last generation to get time of coalescence and ghosts
	//vector<std::set<uint32_t>> childs_ids_last_gen;
	//vector<bitmap> childs_ids_last_gen;
	//vector<roaring::Roaring> childs_ids_last_gen;
	
		
	Population(){
		par_ids.resize( 2 * config::pop_size );		
		chr_choices.resize( 2 * config::nbchr * config::pop_size ); 	
		members.resize( config::pop_size );		
		chr_nb_breakpoints.resize( 2 * config::nbchr * config::pop_size );
	}
	
	
	pair<uint32_t, uint32_t> get_parents_of(uint32_t individ){
		return make_pair( par_ids[2*individ], par_ids[2*individ+1] );
	}
};







#define MAX_ARG_GEN 101

class ARGInfo{
public:
	
	/**
	IMPORTANT: because children are stored in a bitmap, the set of children of an indiv MUST be added in increasing order.
	This is because bitmap wants its bit sets in increasing order, otherwise it breaks.
	Learned this the hard way when I had to debug - consider changing to another data structure at some point.
	Also see comment in set_parents
	**/
	vector< map<uint32_t, bitmap> > children_per_gen;	//children_per_gen[g][i] = set of children of indiv i at gen g
	uint32_t chunk_size;
	
	//don't forget to assign bits in order!
	vector<bitmap> genetic_anc_per_gen;					//genetic_anc_per_gen[g] = set of genetic ancestors per generation

	//vector< map<uint32_t, bitmap> > all_anc_per_chunk;	//all_anc_per_chunk[g][c] = indivs at gen g that are ancestors of everyone in chunk c
	
	
	ARGInfo(){
		set_chunk_size(1000);
		children_per_gen.resize(MAX_ARG_GEN);
		genetic_anc_per_gen.resize(MAX_ARG_GEN);
		//all_anc_per_chunk.resize(MAX_ARG_GEN);
	}
	
	
	
	void set_chunk_size(uint32_t chunk_size){
		this->chunk_size = chunk_size;
	
	}
	
	
	uint32_t get_number_of_chunks(){
		return ceil( (float)config::init_nb_indiv / (float)chunk_size );
	}
	
	
	
	/**
	This tells parent1 and parent2 that indiv is one of their children.
	IMPORTANT: because children are stored in a bitmap, set_parents MUST be called with indiv in increasing order.
	This is because bitmap wants its bit sets in increasing order, otherwise it breaks.
	Learned this the hard way when I had to debug - consider changing to another data structure at some point.
	**/
	void set_parents(uint32_t indiv, uint32_t gen_of_indiv, uint32_t parent1, uint32_t parent2){
		if (gen_of_indiv + 1 >= children_per_gen.size()){
			cout<<"Trying to set parents at gen="<<gen_of_indiv<<" but MAX_ARG_GEN="<<MAX_ARG_GEN<<", will ignore"<<endl;
			return;
		}
		
		children_per_gen[ gen_of_indiv + 1 ][parent1].set(indiv);
		children_per_gen[ gen_of_indiv + 1 ][parent2].set(indiv);
	}
	
	
	
	
	//returnvalue[g] = set of indivs at gen g that are ancestors of everyone in chunk, g up to maxgen
	vector<bitmap> compute_all_anc_chunk(uint32_t chunk_no, uint32_t maxgen){
		
		vector<bitmap> retvec(maxgen + 1);	//don't forget to set bits in order

		//[i1, i2) = range of indivs in chunk
		uint32_t i1 = chunk_no * this->chunk_size;
		uint32_t i2 = min( (chunk_no + 1) * this->chunk_size, config::pop_size );
		uint32_t nb_followed_ind_in_chunk = 0;
		
		
		//cur_desc[i] = set of descendants of indiv i in the range [i1, i2 - 1] in current generation (init gen = 0)
		//updated each iteration
		vector<bitmap> cur_desc;
		cur_desc.resize(config::pop_size);
		
		//at gen = 0, descendants of i is only i itself
		for (uint32_t i = i1; i < i2; ++i){
			if (i < config::init_nb_indiv){
				cur_desc[i].set(i);	//sets bits in order
				nb_followed_ind_in_chunk++;
			}
			// else : not an offspring we want to follow
		}
		
		
		vector<bitmap> next_desc;
		next_desc.resize(config::pop_size);
		
		uint32_t nbanc = 0;
		
		//then update cur_desc for gen=1, 2, ..., gen
		for (uint32_t g = 1; g <= maxgen; ++g){
			
			//the loop below assumes that cur_desc contains the chunk descendants of every indiv for 
			//gen g-1, and puts into next_desc the descendants of gen g for each indiv, by taking the union of desc of their children.
			
			#pragma omp parallel for num_threads(4)
			for (uint32_t i = 0; i < config::pop_size; ++i){
				next_desc[i].reset();
				
				if (children_per_gen[g].count(i)){
					
					bitmap& children = children_per_gen[g][i];
					
					for (uint32_t ch : children){
						next_desc[i] = next_desc[i] | cur_desc[ch];
					}
				}
			}
			
			cur_desc.swap(next_desc);	//omg I use swap so clever, previously was "cur_desc = next_desc;"
			
			
			//if any indiv has the whole chunk as descendants, we record it
			for (uint32_t i = 0; i < config::pop_size; ++i){
				if (cur_desc[i].numberOfOnes() == nb_followed_ind_in_chunk){
					retvec[g].set(i);	//sets bits in order
				}
			}
		}

		return retvec;
	}
	
	
	
	
	
	vector<bitmap> get_all_anc(uint32_t maxgen){
		
		uint32_t nb_chunks = ceil( (float)config::pop_size / (float)chunk_size );
		
		//all_anc[g] = set of indivs that are ancestors of everyone in generation g
		vector<bitmap> all_anc(maxgen + 1);

		
		for (uint32_t c = 0; c < nb_chunks; ++c){

			vector<bitmap> all_anc_in_chunk = compute_all_anc_chunk(c, maxgen);

			for (uint32_t g = 1; g <= maxgen; ++g) {

				//if c == 0, all_anc starts with the all-ancestors of chunk 0 (since it is not initialized)
				//then for c > 0, we retain only those that are all-ancestors so far, and ancestor of the next chunk
				if (c == 0)
					all_anc[g] = all_anc[g] | all_anc_in_chunk[g];
				else 
					all_anc[g] = all_anc[g] & all_anc_in_chunk[g];
			}		
		}

		
		
		return all_anc;
		
	}

};






class Lineage{
public:
	//we maintain two Population objects: the current and a temporary one, for memory purposes.  Their memory can be reused.
	//we make them pointers so that we can easily swap them without accidental copies
	Population* cur_pop;
	Population* temp_pop;
	
	//we also maintain two segment lists, again to reuse their meomory
	SegmentList* cur_seglist;
	SegmentList* temp_seglist;

	//data we want to store and export in the end
	vector<uint32_t> nb_ind_genealogical_ancestors;
	vector<uint32_t> nb_ind_genetic_ancestors;
	vector<uint32_t> nb_chr_genetic_ancestors;
	vector<uint64_t> nb_bases;
	vector<uint32_t> nb_segments;
	vector<uint32_t> nb_fusions;
	vector<uint32_t> nb_seg_coalescences;
	vector<uint32_t> nb_splits;
	vector<uint32_t> super_ghosts;

	// to use at each generation
	uint32_t tmp_nb_segments;
	uint32_t tmp_nb_ind_genetic_anc;
	uint32_t tmp_nb_chr_genetic_anc;
	uint32_t tmp_nb_fusions;
	uint32_t tmp_nb_seg_coalescences;
	uint32_t tmp_nb_splits;
	uint64_t tmp_nb_bases;

	// Genealogical data (mainly for ghosts)
	uint32_t first_common_anc;
	uint32_t all_common_anc;
	ARGInfo arg_info;

	// Follow time
	uint32_t back_time = 0;

	Lineage(){
		back_time = 0;
		arg_info.set_chunk_size(1000);
		cur_pop = new Population();
		std::fill(cur_pop->members.begin(), cur_pop->members.end(), 1);

		temp_pop = new Population();
		
		
		cur_seglist = new SegmentList();
		temp_seglist = new SegmentList();

		for (size_t i = 0; i < config::init_nb_indiv; ++i){
			for (size_t c = 0; c < config::nbchr; ++c){
				cur_seglist->add( i, c, 0, 0, config::chrlen -1 );    
				cur_seglist->add( i, c, 1, 0, config::chrlen -1 );
			}
		}

		// data storage TODO : if we have the information of number of generations, 
		// we could set sizes and prevent further resizing of vectors ? 
		nb_ind_genealogical_ancestors.push_back(config::pop_size);
		nb_ind_genetic_ancestors.push_back(config::pop_size);
		nb_chr_genetic_ancestors.push_back(config::pop_size * config::nbchr * 2);
		super_ghosts.push_back(0);
		nb_segments.push_back(config::pop_size * config::nbchr * 2);
		// ^ hard-coded that we take the whole population to start with
		nb_fusions.push_back(0);
		nb_seg_coalescences.push_back(0);
		nb_splits.push_back(0);
		nb_bases.push_back(config::pop_size * config::nbchr * 2 * config::chrlen);

		// To use at each generation
		tmp_nb_segments = 0;
		tmp_nb_ind_genetic_anc = 0;
		tmp_nb_chr_genetic_anc = 0;
		tmp_nb_fusions = 0;
		tmp_nb_seg_coalescences = 0;
		tmp_nb_splits = 0;
		tmp_nb_bases = 0;

		// init genealogical data (mainly for ghosts): will be replaced with actual value
		first_common_anc = 0;
		all_common_anc = 0;
	}
	
	
	~Lineage(){
		delete cur_pop;
		delete temp_pop;
		delete cur_seglist;
		delete temp_seglist;
	}
	

	void write_data(string filename){
		std::ofstream data(filename);
		data<<"backtime,nb_ind_genealogical_ancestors,nb_ind_genetic_ancestors,nb_chr_genetic_ancestors";
		data<<",nb_segments,nb_fusions,nb_seg_coal,nb_splits,nb_bases,nb_super_ghosts"<<std::endl;
		for (uint32_t t = 0; t < back_time; t++){
			data<<t<<","<<nb_ind_genealogical_ancestors[t]<<","<<nb_ind_genetic_ancestors[t]<<","<<nb_chr_genetic_ancestors[t]<<",";
			data<<nb_segments[t]<<","<<nb_fusions[t]<<","<<nb_seg_coalescences[t]<<","<<nb_splits[t]<<","<<nb_bases[t]<<","<<super_ghosts[t]<<'\n';
		}
		data.close();
	}

	void write_coal_data(string filename){
		std::ofstream coal("coal_"+filename);
		coal<<"first_common_anc_time,all_common_anc_time\n";
		coal<<first_common_anc<<","<<all_common_anc<<std::endl;
		coal.close();
	}
	

	void backstep(){

		//fill the temp_pop with the next_pop.  
		Lineage::generate_population(*cur_pop, *temp_pop);

		tmp_nb_splits = 0;
		//fill the temp_seglist with the new segments
		Lineage::do_recombinations(*cur_pop, *temp_pop, *cur_seglist, *temp_seglist);

		//now that next pop is filled, swap it for the current
		Population* even_more_temporary_pop = cur_pop;
		cur_pop = temp_pop;
		temp_pop = even_more_temporary_pop;
		
		//Also swap segments
		SegmentList* even_more_temporary_seglist = cur_seglist;
		cur_seglist = temp_seglist;
		temp_seglist = even_more_temporary_seglist;
		
		tmp_nb_fusions = 0;
		tmp_nb_seg_coalescences = 0;
		Lineage::check_fused_segments(*cur_seglist);

		Lineage::record_stats_for_this_generation(*cur_pop, *cur_seglist);
		if (all_common_anc == 0){
			// Deal with genealogical data to check coalescence
			Lineage::check_coalescence();
		}
		
		back_time++;
	}


	
	/*
	* Assuming that we never reached a state of only superancestors, 
	* this function:
	* 1) finds the first generation with a superancestor (first_common_anc)
	* 2) finds the first generation with ONLY superancestors (all_common_anc)
	* 3) fills the superghosts array
	* The function replays evolution from gen 0 to back_time, so we only run it at certain generations.
	* If run twice, data from previous calculation is reset.
	*/
	void check_coalescence(){
	
		if (all_common_anc != 0)
			return;	//nothing more to do here, only super-ancestors from now on
	
		if (back_time % 70 == 0 && back_time > 0) {	//70 generations should be enough to reach all_common_anc state, even with 1M pop

			super_ghosts.clear();
			first_common_anc = 0;
			all_common_anc = 0;
			
			vector<bitmap> all_anc = arg_info.get_all_anc(back_time);

			for (uint32_t g = 0; g <= back_time; g++) {


				//debug stuff 
				//cout << "all_anc[g]=" << g << " : " << all_anc[g].numberOfOnes() << 
				//		"   nb_ind_genealogical_ancestors[g]="<< nb_ind_genealogical_ancestors[g]<<endl;


				//check if g is first generation with a super-ancestor (genealogical)
				if (first_common_anc == 0 && !all_anc[g].empty()) {
					first_common_anc = g;

					//cout << "first_common_anc=" << g << endl;
				}


				//check if g is first generation with ONLY super-ancestors
				if (first_common_anc != 0 && all_common_anc == 0) {
					
					//nb gen ancestors should have been computed in record_stats_for_this_generation
					uint32_t nb_gen_ancestors = nb_ind_genealogical_ancestors[g];
					if (all_anc[g].numberOfOnes() == nb_gen_ancestors) {
						all_common_anc = g;

						//cout << "all_common_anc=" << g << endl;
					}
				}


				//then record number of superghosts
				if (first_common_anc == 0) {
					super_ghosts.emplace_back(0);
				}
				else{
					//super ghosts = superancestors but NOT genetic ancestor
					bitmap superghosts = (all_anc[g] - arg_info.genetic_anc_per_gen[g]);
					super_ghosts.emplace_back(superghosts.numberOfOnes());
				}
			}


		}
	
		/*uint32_t max_nb_childs = 0;
		uint32_t min_nb_childs = config::pop_size;
		for (auto s : pop.childs_ids_last_gen){
			uint32_t size = s.size();
			if (size > max_nb_childs) max_nb_childs = size;
			if (size < min_nb_childs) min_nb_childs = size;
		}
		if (first_common_anc == 0 && max_nb_childs == config::pop_size){
			first_common_anc = back_time;
		}
		if (min_nb_childs == config::pop_size){
			all_common_anc = back_time;
		}*/
	}


	void record_stats_for_this_generation(Population &pop, SegmentList &seglist) {
	
		tmp_nb_bases = 0;
		tmp_nb_segments = 0;

		std::set<int> chrsm_ids;
		std::set<int> ind_ids;

		for (uint32_t i = 0; i < 2*config::nbchr; i++){
			tmp_nb_segments += seglist.sizes[i];
			for (uint32_t size = 0; size < seglist.sizes[i]; size++){
				Segment seg = seglist.segments[i][size];
				ind_ids.insert(seg.individ);
				chrsm_ids.insert(2*config::nbchr*seg.individ + 2*seg.chrno + seg.chrindex);
				tmp_nb_bases += seg.b - seg.a + 1;
			}
		}
		
		
		//we need to remember genetic ancestors for when we compute superghosts later on
		if (back_time < MAX_ARG_GEN - 1) {
			
			//must build bitmap for arginfo.  The nice thing about std::set is that enumeration is in sorted order.
			bitmap ind_bitmap;
			for (int id : ind_ids){
				ind_bitmap.set(id);
			}
			
			//why back_time + 1?  because ind_ids are the parents of the upcoming generation, but back_time was not incremented yet
			arg_info.genetic_anc_per_gen[back_time + 1] = ind_bitmap;	
		}
		

		tmp_nb_ind_genetic_anc = ind_ids.size();
		tmp_nb_chr_genetic_anc = chrsm_ids.size();

		nb_fusions.emplace_back(tmp_nb_fusions);
		nb_seg_coalescences.emplace_back(tmp_nb_seg_coalescences);
		nb_splits.emplace_back(tmp_nb_splits);
		nb_segments.emplace_back(tmp_nb_segments);
		nb_ind_genetic_ancestors.emplace_back(tmp_nb_ind_genetic_anc);
		nb_chr_genetic_ancestors.emplace_back(tmp_nb_chr_genetic_anc);
		nb_bases.emplace_back(tmp_nb_bases);
		nb_ind_genealogical_ancestors.emplace_back(std::accumulate(pop.members.begin(), pop.members.end(), 0));

		if (all_common_anc != 0){
			// all genealogical ancestors are the ancestors of everybody
			super_ghosts.emplace_back(nb_ind_genealogical_ancestors.back() - nb_ind_genetic_ancestors.back());
		}
		//below is now handled by check coalescence
		/*else if (first_common_anc == 0){
			// nobody is the ancestor of everybody
			super_ghosts.emplace_back(0);
		}
		else {
			// tricky case
			uint32_t acc_super_ghosts = 0;
			for (uint32_t individ = 0; individ < config::pop_size; individ++){
				//for each individual check if its a superancestor and if its a genealogical ancestor
				if (pop.childs_ids_last_gen[individ].size() == config::pop_size
					&& (ind_ids.find(individ) == ind_ids.end())	){
						acc_super_ghosts++;
					}
			}
			super_ghosts.emplace_back(acc_super_ghosts);
		}*/
	}
	
	//TODO: the next functions are static only because they did not belong to a class before, and this was the easiest way to add them to the Lineage class.
	void generate_population(Population &prev_pop, Population &pop_to_fill) {

		//TODO: could be optimized by computing rand parents only for genealogical indivs
		//chooses two parents for each member of prev_pop

		//args are vector_to_fill, min, max (exclusively), number of entries
		rando_mp::generate_random_integers(prev_pop.par_ids, 0, config::pop_size, 2 * config::pop_size);



		std::fill(pop_to_fill.members.begin(), pop_to_fill.members.end(), 0);
		


		for (uint32_t i = 0; i < config::pop_size; ++i) {
			if (prev_pop.members[i] != 0) {
				pop_to_fill.members[prev_pop.par_ids[2 * i]] = 1;
				pop_to_fill.members[prev_pop.par_ids[2 * i + 1]] = 1;
			}
			if (all_common_anc == 0){
				// We have not reached coalescence yet and need to follow genealogical ancestry
				
				if (back_time < MAX_ARG_GEN - 1){
				
					if (prev_pop.members[i] == 1){
						uint32_t parent1 = prev_pop.par_ids[2 * i];
						uint32_t parent2 = prev_pop.par_ids[2 * i + 1];
						
						//NOTE: SET_PARENT ONLY WORKS BECAUSE WE ARE SENDING i's IN INCREASING ORDER
						//otherwise, bitset data structure would break
						arg_info.set_parents(i, back_time, parent1, parent2);
					}
				}			
			}
		}

	}




	


	
	
	void do_recombinations(Population &prev_pop, Population &next_pop, SegmentList &prev_seglist, SegmentList &next_seglist){
	
		    
		//chooses the recombined chromosome obtained from each parent (0 = first chr, 1 = second chr)
		rando_mp::generate_random_integers(prev_pop.chr_choices, 0, 2, 2 * config::nbchr * config::pop_size);

		
		//TODO: could be optimized by choosing nb breakpts only for genetic chrsms
		//chooses the number of breakpoints for each chromosome
		rando_mp::generate_binomials( prev_pop.chr_nb_breakpoints, config::chrlen, config::recomb_rate, 2 * config::nbchr * config::pop_size );


		/*
		Next we pre-calculate the set of breakpoints on each chr.  nb_breakpoints is the sum of nb break across *all* chr, then we generate all break positions one shot 
		into all_breakpoint_positions.  Each specific chr has its breaks in a sublist of that vector.
		chr_breakpoint_indices[ index_of_chr ] tells where that sublist starts, that is, the set of break of the chr are 
		 all_breakpoint_positions[ chr_breakpoint_indices[ index_of_chr ]  :  chr_breakpoint_indices[ index_of_chr ] + prev_pop.chr_nb_breakpoints[ index_of_chr ] ]
		This is a bit complicated, but this is the price to pay for optimization... 
		*/
		uint64_t nb_breakpoints = accumulate(prev_pop.chr_nb_breakpoints.begin(), prev_pop.chr_nb_breakpoints.begin() + 2 * config::nbchr * config::pop_size, 0);
		vector<uint32_t> all_breakpoint_positions(nb_breakpoints);
		rando_mp::generate_random_integers(all_breakpoint_positions, 1, config::chrlen, nb_breakpoints);
		
		
		size_t cur_index = 0;
		vector<int> chr_breakpoint_indices( 2 * config::nbchr * config::pop_size );
		for (uint32_t chrno = 0; chrno < config::nbchr; ++chrno){
			for (int chrindex = 0; chrindex < 2; ++chrindex){
				for (uint32_t individ = 0; individ < config::pop_size; ++individ){
					size_t index_of_chr = individ * 2 * config::nbchr + 2 * chrno + chrindex;
					
					uint32_t nbrecombs = prev_pop.chr_nb_breakpoints[index_of_chr];
					if (nbrecombs == 0){
						chr_breakpoint_indices[index_of_chr] = -1;	//indicates no recomb
					}
					else{
						chr_breakpoint_indices[index_of_chr] = cur_index;
						cur_index += nbrecombs;
					}				
				}
			}	
		}
		


		
		prev_seglist.sort();
		next_seglist.reset_sizes();
		
		//Since prev_seglist is split into chromosomes, we iterate over each chromosome and handle the segments on it.  Recall that 
		//prev_seglist stores segments on that chromosome for all indivs - we sort the lists by individual, and handle all segments on that indiv in a single loop.
		//When we encouter a new indiv, a new loop starts.
		//Note: same warning on parallelization as in check_fused_segment
		#pragma omp parallel for 
		for (uint32_t chrno = 0; chrno < config::nbchr; ++chrno){
			for (int chrindex = 0; chrindex < 2; ++chrindex){
			
				int next_seg_cpt = 0;
				uint32_t s = 0; 	//current segment on chromosome (c, chrindex)

				
				
				while (s < prev_seglist.get_size(chrno, chrindex)){
					
					//when we get here, we meet the first segment of chrno->chrindex from individ = seg.indiv_id
					//we traverse of segments with that indiv, until we get to another indiv
					Segment &seg = prev_seglist.get(s, chrno, chrindex);
					
					//this is to know when we have reached the segments of a new indiv
					int cur_segindivid = seg.individ;
					int last_segindivid = cur_segindivid;
					
					uint32_t parent1_id, parent2_id;
					std::tie(parent1_id, parent2_id) = prev_pop.get_parents_of(seg.individ);	//sets two variables at once

					size_t index_of_chr = seg.individ * 2 * config::nbchr + 2 * chrno + chrindex;
					
					
					uint32_t nbrecombs = prev_pop.chr_nb_breakpoints[ index_of_chr ];
					int bpindex = chr_breakpoint_indices[ index_of_chr ];


					std::vector<uint32_t> recomb_pos_list;
					if (nbrecombs > 0){
						recomb_pos_list = vector<uint32_t>(all_breakpoint_positions.begin() + bpindex, all_breakpoint_positions.begin() + bpindex + nbrecombs);
						std::sort(recomb_pos_list.begin(), recomb_pos_list.end());
						
					}
					
					/*
					uint32_t nbrecombs = prev_pop.chr_nb_breakpoints[ index_of_chr ];	//hmm indexing has become complicated
					
					
					std::vector<uint32_t> recomb_pos_list(nbrecombs);

					if (nbrecombs > 0){
						rando_mp::generate_random_integers(recomb_pos_list, 1, config::chrlen, nbrecombs);
						std::sort(recomb_pos_list.begin(), recomb_pos_list.end());

					}
					*/
					
					
					
					
					int par_gave_chr, par_id;
					if (chrindex == 0) { // chrsm A : from parent 1
					    par_gave_chr = prev_pop.chr_choices[index_of_chr];
					    par_id = parent1_id;
					} else { // chrsm B : from parent 2
					    par_gave_chr = prev_pop.chr_choices[index_of_chr];
					    par_id = parent2_id;
					}
					

					uint32_t r = 0;	//position in cur_recomb_pos_list        			
		    			
		    			//This is the loop for the current individual.
		    			//We iterate the segment and recomb_pos_list together.  Whenever appropriate, a segment portion to the parent.  After we have passed all of a segment,  
		    			//we increase s to move to the next segment (until we get into another indiv).  When we handle a breakpoint, we increment r.
		    			while (last_segindivid == cur_segindivid){
						bool increment_s = false;
						if (r >= recomb_pos_list.size()){
						    //there are no more recombinations -> unsent part of current segment is given to parent
						    int seg_begin = seg.a;
						    if (recomb_pos_list.size() > 0 and recomb_pos_list[r - 1] > seg.a and recomb_pos_list[r - 1] <= seg.b)
							seg_begin = recomb_pos_list[r - 1];
	    
						    next_seglist.add(par_id, seg.chrno, par_gave_chr, seg_begin, seg.b);
						    increment_s = true;
						}
						else if (recomb_pos_list[r] <= seg.a){
						    //current recombination point is before a -> increment r and swap parent
						    r += 1;
						    par_gave_chr = abs(par_gave_chr - 1);
						}
						else if (recomb_pos_list[r] > seg.b){
						    //current recombination is after b -> unsent part of current segment is given to parent, increment s
						    int seg_begin = seg.a;
						    if (r > 0 and recomb_pos_list[r - 1] > seg.a and recomb_pos_list[r - 1] <= seg.b)
							seg_begin = recomb_pos_list[r - 1];
						    
						    next_seglist.add(par_id, seg.chrno, par_gave_chr, seg_begin, seg.b);
						    increment_s = true;
						}
						else{   //recomb_pos_list[r] > seg.a and <= seg.b
						    //current recombination is right inside seg -> send part to parent, increment r
							if (r > 0 and recomb_pos_list[r - 1] == recomb_pos_list[r]){
								r += 1; //nothing happens
							}
							else {
								tmp_nb_splits++;
								int seg_begin = seg.a;
								if (r > 0 and recomb_pos_list[r - 1] > seg.a and recomb_pos_list[r - 1] <= seg.b)
								seg_begin = recomb_pos_list[r - 1];

								next_seglist.add(par_id, seg.chrno, par_gave_chr, seg_begin, recomb_pos_list[r] - 1);

								r += 1;
								par_gave_chr = abs(par_gave_chr - 1);
							}
						}
						
						if (increment_s){
						    s++;
						    					    uint32_t preva = seg.a, prevb = seg.b;
						    if (s < prev_seglist.get_size(chrno, chrindex)){
							seg = prev_seglist.get(s, chrno, chrindex);
							cur_segindivid = seg.individ;
						    }
						    else
							cur_segindivid = -999;	//forces while loop to stop
						}
					}
				}
			}

		}	//number of brackets indicates that function could be split
			//cout<<"is="<<indivs.size()<<" c0="<<cpt0<<" c1="<<cpt1<<endl;
			//pause();
	}
	
	
	
	
	void check_fused_segments(SegmentList &cur_seglist){
	
		cur_seglist.sort();

		
		//Note: parallelizing this for loop works here because each iteration affects a different vector in cur_seglist.
		//	If the storing in cur_seglist ever changes, running this loop in parallel may become problematic, e.g. if distinct iterations
		//      affect the same vector in cur_seglist.
		#pragma omp parallel for 
		for (uint32_t chrno = 0; chrno < config::nbchr; ++chrno){
			for (uint32_t chrindex = 0; chrindex < 2; ++chrindex){
				vector<size_t> indices_to_delete;        	
			
				for (size_t i = 1; i < cur_seglist.get_size(chrno, chrindex); ++i){
					//if seg[i] is right next to seg[i-1], on the same indiv and chrsm, fuse it --> index i - 1 will get deleted later on
					Segment& seg = cur_seglist.get(i, chrno, chrindex);

					//cout<<"segbef.a="<<seg.a<<" segbef.b="<<seg.b<<endl;
					Segment& segprev = cur_seglist.get(i - 1, chrno, chrindex);
			    		if (seg.individ == segprev.individ && seg.a <= segprev.b+1){
							if (segprev.b+1 == seg.a){ tmp_nb_seg_coalescences++; }
							else { tmp_nb_fusions++; }
						
						//seg.a = segprev.a;
						//seg.b = max(seg.b,segprev.b);
						
						cur_seglist.set(i, seg.individ, seg.chrno, seg.chrindex, segprev.a, max(seg.b, segprev.b));
						
						indices_to_delete.push_back(i - 1);
			    		}
				}

				
				cur_seglist.remove_indices(chrno, chrindex, indices_to_delete);	
					
			}
		}
		//cout<<"nbdels="<<nbdels<<endl;
	}
	
	
	

	
};


void print_help(){
	std::cout<<"Command line usage:\n"
	<<"./simchr -c nb_of_chrsm -l chrsm_len -g nb_generations -p population_size -i init_nb_indiv -r recombination_rate -s seed_for_prng\n"
	<<"or \n"
	<<"./simchr --nbchr nb_of_chrsm --chrlen chrsm_len --nb_gen nb_generations --pop_size population_size --init_nb_indiv nb_of_indiv_to_init_segments --recomb_rate recombination_rate --seed seed_for_prng\n"
	<<"or any combination of short/long name ! No parameter is mandatory. In absence of specifications, default values are:\n"
	<<" nbchr : "<<config::nbchr<<"\n"
	<<" chrlen : "<<config::chrlen<<"\n"
	<<" nb_gen : "<<config::nb_gen<<"\n"
	<<" pop_size : "<<config::pop_size<<"\n"
	<<" init_nb_indiv : whole population.\n"
	<<" recomb_rate : "<<config::recomb_rate<<"\n"
	<<" seed : "<<config::seed<<"\n"
	<<"\nIf nb_gen provided is 0, simulation will run until it approaches the equilibrium.\n";
}


void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char * options_list = "hc:l:g:p:r:s:i:";
  static struct option long_options_list[] = {
      {"help",      no_argument,        nullptr, 'h'},
      {"nbchr",     required_argument,  nullptr, 'c'},
      {"chrlen",    required_argument,  nullptr, 'l'},
      {"nb_gen",    required_argument,  nullptr, 'g'},
      {"pop_size",  required_argument,  nullptr, 'p'},
	  {"recomb_rate",  required_argument,  nullptr, 'r'},
	  {"seed",      required_argument,  nullptr, 's'},
	  {"init_nb_indiv", required_argument, nullptr, 'i'}
	//   {"recomb_nb",  no_argument,  nullptr, 'R'},
  };

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list,
                               long_options_list, nullptr)) != -1) {
    switch (option) {
      case 'h' : {
        print_help();
        exit(EXIT_SUCCESS);
      }
      case 'c' : {
        config::nbchr = atol(optarg);
		break;
      }
      case 'l' : {
        config::chrlen = atol(optarg);
        break;
      }
      case 'g' : {
        config::nb_gen = atol(optarg);
        break;
      }
	  case 'p' : {
        config::pop_size = atol(optarg);
        break;
      }
      case 'r' : {
        config::recomb_rate = atof(optarg);
        break;
      }
	  case 's' : {
        config::seed = atol(optarg);
        break;
      }
	  case 'i' : {
		config::init_nb_indiv = atol(optarg);
	  }
    }
  }
}



int main(int argc, char* argv[]) {

	interpret_cmd_line_options(argc, argv);
	if (config::init_nb_indiv == 0){
		std::cout<<"init_nb_indiv not specified, follow whole population"<<std::endl;
		// If not specified, we use the whole population as starting point
		config::init_nb_indiv = config::pop_size;
	}
	
	rando_mp::init(config::seed);



	uint32_t step_print = 10;

	Lineage lineage;
	cout<<"Starting simulation   s="<<lineage.cur_seglist->get_total_nb_segments()<<"   nbases="<<lineage.cur_seglist->get_total_segment_size()<<endl;

	bool should_continue = true;
	float threshold_factor = 1;
	uint32_t target_bases = config::nbchr * 2 * config::chrlen * threshold_factor;
	if (config::nb_gen == 0){
		std::cout << "will end when nb_bases reach " << target_bases<<std::endl;
	}

	while ( should_continue ) {

		lineage.backstep();

		if (lineage.back_time % step_print == 0) {
			size_t nbsegments = lineage.nb_segments.back();
			size_t segsizes = lineage.nb_bases.back();
			size_t npop = lineage.nb_ind_genealogical_ancestors.back();
			cout << "g=" << lineage.back_time << "   npop=" << npop << "   nbseg="<<nbsegments<<"   nbbases="<<segsizes<<endl;
		}
		if ((config::nb_gen != 0) && (lineage.back_time >= config::nb_gen)){
			should_continue = false;
		}
		if ((config::nb_gen == 0) && (lineage.nb_bases.back() < target_bases)){
			should_continue = false;
		}
	}
	stringstream filename;
	filename << "nbchr-"<<config::nbchr<<"-chrlen-"<<config::chrlen<<"-nb_gen-"<<config::nb_gen
	<<"-pop_size-"<<config::pop_size<<"-init_nb_indiv-"<<config::init_nb_indiv<<"-recomb_rate-"
	<<config::recomb_rate<<"-seed-"<<config::seed<<".csv";
	lineage.write_data(filename.str());
	lineage.write_coal_data(filename.str());

	return 0;
}






