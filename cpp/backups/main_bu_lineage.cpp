#include <iostream>
#include <vector>
#include <random>
//backup with lineage class

#include <algorithm>
#include <array>
#include <set>

#include <stdint.h>

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
    uint32_t nb_gen = 200;
    uint32_t pop_size = 100;
    double recomb_rate = 1/((double)chrlen);	
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
        
        Segment(uint32_t individ, uint32_t chrno, uint32_t chrindex, uint32_t a, uint32_t b) : individ(individ), chrno(chrno), chrindex(chrindex), a(a), b(b){
        	
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











class Lineage{
public:
	//we maintain two Population objects: the current and a temporary one, for memory purposes.  Their memory can be reused.
	//we make them pointers so that we can easily swap them without accidental copies
	Population* cur_pop;
	Population* temp_pop;
	
	//we also maintain two segment lists, again to reuse their meomory
	SegmentList* cur_seglist;
	SegmentList* temp_seglist;


	Lineage(){
		
		cur_pop = new Population();
		std::fill(cur_pop->members.begin(), cur_pop->members.end(), 1);

		temp_pop = new Population();
		
		
		cur_seglist = new SegmentList();
		temp_seglist = new SegmentList();

		for (size_t i = 0; i < config::pop_size; ++i){
			for (size_t c = 0; c < config::nbchr; ++c){
				cur_seglist->add( i, c, 0, 0, config::chrlen -1 );    
				cur_seglist->add( i, c, 1, 0, config::chrlen -1 );
			}
		}
	}
	
	
	~Lineage(){
		delete cur_pop;
		delete temp_pop;
		delete cur_seglist;
		delete temp_seglist;
	}
	
	
	
	

	void backstep(){

		//fill the temp_pop with the next_pop.  
		Lineage::generate_population(*cur_pop, *temp_pop);

		//fill the temp_seglist with the new segments
		//Lineage::do_recombinations(*cur_pop, *temp_pop, *cur_seglist, *temp_seglist);

		//now that next pop is filled, swap it for the current
		Population* even_more_temporary_pop = cur_pop;
		cur_pop = temp_pop;
		temp_pop = even_more_temporary_pop;
		
		//Also swap segments
		SegmentList* even_more_temporary_seglist = cur_seglist;
		cur_seglist = temp_seglist;
		temp_seglist = even_more_temporary_seglist;
		
		
		//Lineage::check_fused_segments(*cur_seglist);
		
	}
	
	
	
	static void generate_population(Population &prev_pop, Population &pop_to_fill) {

		//TODO: could be optimized by computing rand parents only for genealogical indivs
		//chooses two parents for each member of prev_pop

		//args are vector_to_fill, min, max (exclusively), number of entries
		rando_mp::generate_random_integers(prev_pop.par_ids, 0, config::pop_size, 2 * config::pop_size);


		for (int i = 0; i < 2 * config::pop_size; ++i){
			cout<<prev_pop.par_ids[i]<<" ";
		}
		cout<<endl;


		std::fill(pop_to_fill.members.begin(), pop_to_fill.members.end(), 0);

				uint32_t npop = accumulate(pop_to_fill.members.begin(), pop_to_fill.members.end(), 0);

		for (uint32_t i = 0; i < config::pop_size; ++i) {
			if (prev_pop.members[i] != 0) {
				pop_to_fill.members[prev_pop.par_ids[2 * i]] = 1;
				pop_to_fill.members[prev_pop.par_ids[2 * i + 1]] = 1;
			}
		}
				uint32_t npop2 = accumulate(pop_to_fill.members.begin(), pop_to_fill.members.end(), 0);
				cout<<"n="<<npop<<" n2="<<npop2<<endl;
	}




	


	
	
	static void do_recombinations(Population &prev_pop, Population &next_pop, SegmentList &prev_seglist, SegmentList &next_seglist){
	
		    
		//chooses the recombined chromosome obtained from each parent (0 = first chr, 1 = second chr)
		rando_mp::generate_random_integers(prev_pop.chr_choices, 0, 2, 2 * config::nbchr * config::pop_size);

		
		//TODO: could be optimized by choosing nb breakpts only for genetic chrsms
		//chooses the number of breakpoints for each chromosome
		rando_mp::generate_binomials( prev_pop.chr_nb_breakpoints, config::chrlen, config::recomb_rate, 2 * config::nbchr * config::pop_size );


		
		prev_seglist.sort();
		next_seglist.reset_sizes();
		
		set<uint32_t> indivs; 
				
		
		//Since prev_seglist is split into chromosomes, we iterate over each chromosome and handle the segments on it.  Recall that 
		//prev_seglist stores segments on that chromosome for all indivs - we sort the lists by individual, and handle all segments on that indiv in a single loop.
		//When we encouter a new indiv, a new loop starts.
		//Note: same warning on parallelization as in check_fused_segment
		//#pragma omp parallel for 
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


					//here we choose and sort the number of breakpoints in chromosome from seg.individ, chrno, chrindex
					//could be optimized by doing batch random choices...
					uint32_t nbrecombs = prev_pop.chr_nb_breakpoints[ seg.individ * 2 * chrno + 2 * chrno + chrindex ];	//hmm indexing has become complicated
					//nbrecombs++;	//to do like the python code
					std::vector<uint32_t> recomb_pos_list(nbrecombs);

					if (nbrecombs > 0){
						rando_mp::generate_random_integers(recomb_pos_list, 1, config::chrlen, nbrecombs);
						std::sort(recomb_pos_list.begin(), recomb_pos_list.end());
					}
					
					
					int par_gave_chr, par_id;
					if (chrindex == 0) { // chrsm A : from parent 1
					    par_gave_chr = prev_pop.chr_choices[seg.individ * 2 * chrno + 2 * chrno + chrindex];	//not sure about indices
					    par_id = parent1_id;
					} else { // chrsm B : from parent 2
					    par_gave_chr = prev_pop.chr_choices[seg.individ * 2 * chrno + 2 * chrno + chrindex];
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
						    int seg_begin = seg.a;
						    if (r > 0 and recomb_pos_list[r - 1] > seg.a and recomb_pos_list[r - 1] <= seg.b)
							seg_begin = recomb_pos_list[r - 1];
						    
						    next_seglist.add(par_id, seg.chrno, par_gave_chr, seg_begin, recomb_pos_list[r] - 1);

						    r += 1;
						    par_gave_chr = abs(par_gave_chr - 1);
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
	
	
	
	
	static void check_fused_segments(SegmentList &cur_seglist){
	
		cur_seglist.sort();

		int nbdels = 0;
		
		//Note: parallelizing this for loop works here because each iteration affects a different vector in cur_seglist.
		//	If the storing in cur_seglist ever changes, running this loop in parallel may become problematic, e.g. if distinct iterations
		//      affect the same vector in cur_seglist.
		//#pragma omp parallel for 
		for (uint32_t chrno = 0; chrno < config::nbchr; ++chrno){
			for (uint32_t chrindex = 0; chrindex < 2; ++chrindex){
				vector<size_t> indices_to_delete;        	
			
				for (size_t i = 1; i < cur_seglist.get_size(chrno, chrindex); ++i){
					//if seg[i] is right next to seg[i-1], on the same indiv and chrsm, fuse it --> index i - 1 will get deleted later on
					Segment& seg = cur_seglist.get(i, chrno, chrindex);

					//cout<<"segbef.a="<<seg.a<<" segbef.b="<<seg.b<<endl;
					Segment& segprev = cur_seglist.get(i - 1, chrno, chrindex);
			    		if (seg.individ == segprev.individ && seg.a <= segprev.b+1){
						
						//seg.a = segprev.a;
						//seg.b = max(seg.b,segprev.b);
						
						cur_seglist.set(i, seg.individ, seg.chrno, seg.chrindex, segprev.a, max(seg.b, segprev.b));
						
						indices_to_delete.push_back(i - 1);
						nbdels++;
			    		}
				}

				
				cur_seglist.remove_indices(chrno, chrindex, indices_to_delete);	
					
			}
		}
		//cout<<"nbdels="<<nbdels<<endl;
	}
	
	
	

	
};









int main() {

	/*
	these set the seed for the randomizers
	rando_mp needs and int, plus four uint32_t
	rando_xo needs four uint32_t for seeds, but it is set by rando_mp
	*/
	rando_xo::init(11,22,33,44);
	rando_mp::init(12222223, 11, 22, 33, 44);

	Lineage lineage;
	cout<<"Starting simulation   s="<<lineage.cur_seglist->get_total_nb_segments()<<endl;


	for (uint32_t g = 0; g < config::nb_gen; ++g) {

		lineage.backstep();

		if (g % 1 == 0) {
		    size_t nbsegments = lineage.cur_seglist->get_total_nb_segments();
		    uint32_t npop = accumulate(lineage.cur_pop->members.begin(), lineage.cur_pop->members.end(), 0);
		    cout << "g=" << g << "   npop=" << npop << "   s="<<nbsegments<<endl;
		}
	}


	return 0;
}





