import numpy as np
import random as random
#import json
import os
import sys

#from numba import jit
#from numba.experimental import jitclass
#from numba import int64, int32 

#from line_profiler import profile
#from execo_engine import slugify


'''
cython3 -3 --embed recomb_multichrsm.pyx
gcc -O3 -I/usr/include/python3.10 -I/usr/include/python3.10 recomb_multichrsm.c -lpython3.10 -o test
'''


#TODO: il y a toujours au moins une recomb par chr, même on diminue le rate.  Est-ce ok?





pop_size = 20000       # nb of individuals (of 2*nb_chrsm chromosomes each)
final_time = 1000       # Simulation length in generation 
chrsm_len = 500000       # "blocks" in *each* chromosome
recomb_rate = 1/chrsm_len   # probability to recombine between any 2 blocks
nb_chrsm = 36          # nb of *pairs* of chromosomes


verbose = True
write_data = False
filename = "tmp.json"
seed = 47098438



## prepare the random draws
np.random.seed(seed)
chrsm_choice_rd = np.random.randint(0, 2, nb_chrsm*pop_size)
expected_nb_recomb = int(pop_size*chrsm_len*recomb_rate*2*nb_chrsm)
recomb_counter = 0
recomb_pos_random = np.random.randint(1,
                chrsm_len,
                2*expected_nb_recomb)




INDIVID = 0
CHRNO = 1
CHRINDEX = 2
SEG_A = 3
SEG_B = 4


'''
class Simulation:
    def __init__(self, dictionary):
        ## default values
        self.final_time = 1000       # Simulation length in generation 
        self.chrsm_len = 10000       # "blocks" in *each* chromosome
        self.recomb_rate = 1/10000   # probability to recombine between any 2 blocks
        self.nb_chrsm = 1            # nb of *pairs* of chromosomes
        self.pop_size = 1000         # nb of individuals (of 2*nb_chrsm chromosomes each)

        self.seg_list = []           # additionnal parameter to get the segments to follow
                                     # if empty, follow the whole population
        self.verbose = False
        self.write_data = True
        self.filename = "tmp.json"
        self.seed = 34968374

        ## override with given values
        for k,v in dictionary.items():
            setattr(self, k, v)
            print(k,v)

        ## prepare the random draws
        np.random.seed(self.seed)
        self.chrsm_choice_rd = np.random.randint(0, 2, self.nb_chrsm*self.pop_size)
        self.expected_nb_recomb = int(self.pop_size*self.chrsm_len*self.recomb_rate*2*self.nb_chrsm)
        self.recomb_counter = 0
        self.recomb_pos_random = np.random.randint(1,
                        self.chrsm_len,
                        2*self.expected_nb_recomb)


    def run(self):
        pop = Population(self)
        for i in range(self.pop_size):
            pop.individuals[i] = Individual(self, i, None)
            # We initialize the whole population anyway to follow the genealogy and get the
            # genealogical ancestors
        lin = Lineage(pop, self.seg_list)

        while (lin.back_time < self.final_time-1):
            if (self.verbose) & (lin.back_time%2==0):
               print(f"T={lin.back_time}  S={sum(len(lc) for lc in lin.cur_segments)} "\
                     f" Ngeneal={lin.nb_ind_genealogical_ancestors[lin.back_time]}  Ngeneti={lin.nb_ind_genetic_ancestors[lin.back_time]}"\
                     f" Nb_bases={lin.genetic_mat[lin.back_time]}")

            lin.backward_step()
            lin.check_fused_segments()

        if self.write_data:
            lin.write_data(self.filename)

'''


class Individual:
    
    
    
    def __init__(self, id):
        

        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random

        # In case of multi chromosomes, we assume they all have the same length and the same average
        # recomb_rate per chromosomes
#        nb_recombs = [1] * nb_chrsm  #np.random.binomial(chrsm_len, recomb_rate, nb_chrsm)
#        nb_recombs2 = [1] * nb_chrsm   #np.random.binomial(chrsm_len, recomb_rate, nb_chrsm)

        nb_recombs = np.random.binomial(chrsm_len, recomb_rate, nb_chrsm)
        nb_recombs2 = np.random.binomial(chrsm_len, recomb_rate, nb_chrsm)

        # Positions de recombinaisons
        self.recomb_par1, self.recomb_par2 = [], []
        for i in range(nb_chrsm):
            self.recomb_par1 += [np.random.randint(1, chrsm_len, nb_recombs[i])]    #+1 ?
            self.recomb_par2 += [np.random.randint(1, chrsm_len, nb_recombs2[i])]   #+1 ?
       
        self.par1_gave_chrsm = chrsm_choice_rd[2*nb_chrsm*id : 2*nb_chrsm*id+nb_chrsm]
        self.par2_gave_chrsm = chrsm_choice_rd[2*nb_chrsm*id+nb_chrsm : (2*nb_chrsm*id+nb_chrsm)+nb_chrsm]
        

        # Les chromosomes, après recombinaison, sont à "droite" ou à "gauche"
        #                  (gauche) 0  1 (droite)
        # |  (                      |  (
        # |  (                      |  (
        # | \(  recombinaision      (  |
        # |  (                      (  |
        # |  (                      (  |
        # Donc si on donne le gauche (0), c’est anciennement ce qui était à gauche 
        # au-dessus du point de coupure, et ce qui était à droite en-dessous du point
        # de coupure

        self.id = id # position dans la population


class Population:
    def __init__(self):

        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        self.individuals = {}
        chrsm_choice_rd = np.random.randint(0,2, 2 * nb_chrsm * pop_size)




    #this is the previous forward_step, now it generates the parent indivs from the child indivs
    #side effect: child_pop indivs get parent1 and parent2 assigned
    def generate_from_child_population(self, child_pop, to_check):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        chrsm_choice_rd = np.random.randint(0,2, 2 * nb_chrsm * pop_size)
        par_ids = np.random.randint(0, pop_size, 2*pop_size)
        recomb_pos_random = np.random.randint(1,
                        chrsm_len,
                        2*expected_nb_recomb)
        recomb_counter = 0

        for id, indiv in child_pop.individuals.items():
            par1_id =  par_ids[2*id]
            par2_id = par_ids[2*id+1]


            if par1_id not in self.individuals:
                self.individuals[par1_id] = Individual(par1_id)

            if par2_id not in self.individuals:
                self.individuals[par2_id] = Individual(par2_id)


            indiv.parent1_id = par1_id
            indiv.parent2_id = par2_id


class Segment:
    __slots__ = 'indiv_id', 'chrsm', 'a', 'b'
    def __init__(self,id,chrsm,a,b):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        self.indiv_id = id
        self.chrsm = chrsm
        self.a = a
        self.b = b

    #def print(self):
    #    print(self.indiv_id, self.chrsm, self.a, self.b)


class Lineage:
    def __init__(self, initial_pop):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        self.pop = initial_pop
        self.back_time = 0

        #if seg_list != []:
        #    self.cur_segments = seg_list # [ Segment(indiv_id, chrsm, a, b)]
        #else: # We follow the whole population
        
        
        if True:
            self.test_segments = np.empty( (pop_size * nb_chrsm * 2, 5), np.int32 )
            seg_index = 0
            for indiv_id in range(pop_size):
                for nc in range(nb_chrsm):
                    for chrindex in [0, 1]:
                        self.test_segments[seg_index] = [indiv_id, nc, chrindex, 0, chrsm_len-1]
                        seg_index += 1
                        
            print(self.test_segments[0])   
            print(self.test_segments[1])                        
            print(self.test_segments[2454])
            print(self.test_segments[2455])
            print(self.test_segments[2456])
            
            
            self.cur_segments = []
            for nc in range(nb_chrsm):
                self.cur_segments.append([])
                for indiv_id in range(pop_size):
                    self.cur_segments[nc].append(Segment(indiv_id, (nc,0), 0, chrsm_len-1))
                    self.cur_segments[nc].append(Segment(indiv_id, (nc,1), 0, chrsm_len-1))
            
            
            
        '''    
        self.nb_segments = np.zeros(final_time)
        self.nb_segments[0] = sum(len(seg_nc) for seg_nc in self.cur_segments)
        anc_ind_list = set()
        anc_chrsm_list = set()
        nb_bases = 0
        for nc in range(nb_chrsm):
            for seg in self.cur_segments[nc]:
                nb_bases += seg.b - seg.a +1
                anc_ind_list.add(seg.indiv_id)
                anc_chrsm_list.add((seg.indiv_id, seg.chrsm))
        self.nb_ind_genetic_ancestors =  np.zeros(final_time)
        self.nb_ind_genetic_ancestors[0] = len(anc_ind_list)
        self.nb_ind_genealogical_ancestors = np.zeros(final_time)
        self.nb_ind_genealogical_ancestors[0] = len(anc_ind_list)
        self.nb_chr_genetic_ancestors = np.zeros(final_time)
        self.nb_chr_genetic_ancestors[0] = len(anc_chrsm_list)

        self.genetic_mat = np.zeros(final_time)
        self.genetic_mat[0] = nb_bases

        self.has_separated = False
        self.first_common_anc = 0
        self.all_common_anc = 0
        '''


    def legacy_backstep(self):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        next_segments = []                
                
        #at this point, each indiv in self.pop must has parent1_id and parent2_id set
        for nc in range(nb_chrsm):
            next_segments.append([])
            for segment in self.cur_segments[nc]:
            
                #ML: NOTE for whatever reason, without the .copy() sorting could affect the results
                if segment.chrsm[1] == 0: # chrsm A : from parent 1
                    par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par1_gave_chrsm[nc]
                    par_id =  self.pop.individuals[segment.indiv_id].parent1_id
                    recomb_pos_list =  self.pop.individuals[segment.indiv_id].recomb_par1[nc].copy()
                else: # chrsm B : from parent 2
                    par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par2_gave_chrsm[nc]
                    par_id =  self.pop.individuals[segment.indiv_id].parent2_id
                    recomb_pos_list =  self.pop.individuals[segment.indiv_id].recomb_par2[nc].copy()

                
                recomb_pos_list.sort()
                
                
                if self.back_time == 1 and par_id == 0 and nc == 2:
                    print(f"sid={segment.indiv_id} schr={segment.chrsm}  pos={recomb_pos_list}")
                
                nb_recomb_before_a, nb_recomb_in_seg = 0,0
                for pos in recomb_pos_list:
                    if pos <= segment.a:
                        nb_recomb_before_a += 1
                    elif pos <= segment.b:
                        nb_recomb_in_seg += 1
                    else:
                        break


                if par_gave_chrsm == 0:
                    # Le parent a donné le chrsm de gauche après la recombinaison
                    start_chrsm = (nb_recomb_before_a % 2) == 1
                else:
                    # Le parent a donné le chrsm de droite après la recombinaison
                    start_chrsm = (nb_recomb_before_a % 2) == 0

                if start_chrsm == True:
                    start_chrsm = 1
                if start_chrsm == False:
                    start_chrsm = 0

                nb_recomb_in_seg = len([pos for pos in recomb_pos_list if (pos > segment.a and pos <= segment.b)])
                if nb_recomb_in_seg == 0:
                    # Pas de coupure du segment suivi
                    # le chromosome où on l’envoie dépend seulement du nombre de coupures *avant*
                    next_segments[nc].append(Segment(par_id, (nc, start_chrsm), segment.a, segment.b))
                else: # le segment est coup en x morceaux
                    self.has_separated = True
                    positions = [segment.a] + list(recomb_pos_list[nb_recomb_before_a:nb_recomb_before_a+nb_recomb_in_seg]) + [segment.b+1]

                    for i in range(nb_recomb_in_seg + 1):
                        if positions[i+1] < positions[i]:
                            print(par_id, nb_recomb_in_seg ,start_chrsm, positions, int((len([pos for pos in recomb_pos_list if pos <= segment.a]) % 2 ) == 1) )
                        next_segments[nc].append(Segment(par_id, (nc, start_chrsm), positions[i], positions[i+1]-1))
                        start_chrsm = abs(start_chrsm - 1)

        
        self.cur_segments = next_segments
    
    
    
    
    
    def backward_step(self):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random

        INDIVID, CHRNO, CHRINDEX, A, B = 0, 1, 2, 3, 4

        next_pop = Population()
        next_pop.generate_from_child_population(self.pop, False) #self.all_common_anc == 0)


        
        next_seg_cpt = 0
        test_next_segments = np.empty( (len(self.test_segments) * 100, 5), np.int32 )  #next going backwards in time
        
        self.test_segments.view('i4,i4,i4,i4,i4').sort(order=['f0'], axis=0)
        #self.test_segments.sort()
        

        
        
        
        last_unit = (-1, -1, -1)	#this tuple represents last (indiv_id, chrno, chrindex) seen in the loop, so we know when it changes
        
        cur_par_gave_chr = -1
        cur_par_id = -1
        cur_recomb_pos_list = None
        
        s = 0
        

	#recall that a segment object is a tuple (individ, chrno, chrindex, a, b)
        while s < len(self.test_segments):
            seg = self.test_segments[s]
            last_unit = (seg[INDIVID], seg[CHRNO], seg[CHRINDEX])
            cur_unit = last_unit
            

            
            if seg[CHRINDEX] == 0: # chrsm A : from parent 1
                cur_par_gave_chr = self.pop.individuals[seg[INDIVID]].par1_gave_chrsm[seg[CHRNO]]
                cur_par_id =  self.pop.individuals[seg[INDIVID]].parent1_id
                cur_recomb_pos_list = self.pop.individuals[seg[INDIVID]].recomb_par1[seg[CHRNO]]
            else: # chrsm B : from parent 2
                cur_par_gave_chr = self.pop.individuals[seg[INDIVID]].par2_gave_chrsm[seg[CHRNO]]
                cur_par_id =  self.pop.individuals[seg[INDIVID]].parent2_id
                cur_recomb_pos_list = self.pop.individuals[seg[INDIVID]].recomb_par2[seg[CHRNO]]
            

            cur_recomb_pos_list.sort()
            
            r = 0	#position in cur_recomb_pos_list
            
            #TODO: what if cur_recomb_pos_list has len 0?
            
            while last_unit == cur_unit:
                news = s
                if r >= len(cur_recomb_pos_list):
                    #there are no more recombinations -> unsent part of current segment is given to parent
                    seg_begin = seg[A]
                    if len(cur_recomb_pos_list) > 0 and cur_recomb_pos_list[r - 1] > seg[A] and cur_recomb_pos_list[r - 1] <= seg[B]:
                        seg_begin = cur_recomb_pos_list[r - 1]
                    
                    test_next_segments[next_seg_cpt] = [cur_par_id, seg[CHRNO], cur_par_gave_chr, seg_begin, seg[B]]
                    next_seg_cpt += 1
                    
                    news = s + 1

                elif cur_recomb_pos_list[r] <= seg[A]:
                    #current recombination point is before a -> increment r and swap parent
                    r += 1
                    cur_par_gave_chr = abs(cur_par_gave_chr - 1)
                elif cur_recomb_pos_list[r] > seg[B]:
                    #current recombination is after b -> unsent part of current segment is given to parent, increment s
                    seg_begin = seg[A]
                    if r > 0 and cur_recomb_pos_list[r - 1] > seg[A] and cur_recomb_pos_list[r - 1] <= seg[B]:
                        seg_begin = cur_recomb_pos_list[r - 1]
                    
                    test_next_segments[next_seg_cpt] = [cur_par_id, seg[CHRNO], cur_par_gave_chr, seg_begin, seg[B]]
                    next_seg_cpt += 1
                    
                    news = s + 1

                else:   # > seg[A] and <= seg[B]
                    #current recombination is right inside seg -> sent part to parent, increment r
                    seg_begin = seg[A]
                    if r > 0 and cur_recomb_pos_list[r - 1] > seg[A] and cur_recomb_pos_list[r - 1] <= seg[B]:
                        seg_begin = cur_recomb_pos_list[r - 1]
                    
                    test_next_segments[next_seg_cpt] = [cur_par_id, seg[CHRNO], cur_par_gave_chr, seg_begin, cur_recomb_pos_list[r] - 1]
                    next_seg_cpt += 1
                    
                    
                    r += 1
                    cur_par_gave_chr = abs(cur_par_gave_chr - 1)
                
                if news != s:
                    s = news
                    if s < len(self.test_segments):
                        seg = self.test_segments[s]
                        cur_unit = (seg[INDIVID], seg[CHRNO], seg[CHRINDEX])                    
                    else:
                        cur_unit = None
                
        #test_next_segments.resize(next_seg_cpt)
        self.test_segments = np.resize(test_next_segments, (next_seg_cpt, 5))

        
        #self.test_segments = np.sort(self.test_segments.view('i4,i4,i4,i4,i4'), order=['f0'], axis=0).view(np.int32)
        #self.test_segments.view('i4,i4,i4,i4,i4').sort(order=['f0'], axis=0)

        #for i in range(len(self.test_segments)):
        #    print(self.test_segments[i])
        #print(self.test_segments)
        #sys.exit()             
        
        
        
        #FOR DEBUGGING
        #self.legacy_backstep()
        global DEBUG
        if DEBUG:
            self.legacy_backstep()
            tmp_segments = []
            for nc in range(nb_chrsm):
                for seg in self.cur_segments[nc]:
                    tmp_segments.append( (seg.indiv_id, seg.chrsm[0], int(seg.chrsm[1]), seg.a, seg.b) )
            tmp_segments.sort()

            self.test_segments.view('i4,i4,i4,i4,i4').sort(order=['f0'], axis=0)
        
            if len(self.test_segments) != len(tmp_segments):
                print(f"difflen: test_seg={len(self.test_segments)} vs tmp={len(tmp_segments)}")
                for i in range(100): #range(len(self.test_segments)):
                     issame = True
                     for j in range(5):
                         if self.test_segments[i,j] != tmp_segments[i][j]:
                             issame = False
                     if not issame or True:
                         print(f"test[{i}] = {self.test_segments[i]}   \nseg[{i}] =  {tmp_segments[i]}")
            else:
                for i in range(min(200, len(self.test_segments))):
                     issame = True
                     for j in range(5):
                         if self.test_segments[i,j] != tmp_segments[i][j]:
                              print(f"test[i] = {self.test_segments[i]}   \nseg[i] =  {tmp_segments[i]}")


        self.back_time += 1
        self.pop.individuals = next_pop.individuals



        '''
        if self.all_common_anc == 0:
            min_nb_desc, max_nb_desc = pop_size, 0
            for indiv in self.pop.individuals.values():
                indiv.final_desc = list(set(indiv.final_desc))
                nb_desc = len(indiv.final_desc)
                if nb_desc < min_nb_desc:
                    min_nb_desc = nb_desc
                if nb_desc > max_nb_desc:
                    max_nb_desc = nb_desc
            if self.first_common_anc == 0 and max_nb_desc == pop_size:
#                print("Most recent common ancestor of the population at time ", self.back_time)
                self.first_common_anc = self.back_time
            if min_nb_desc == pop_size:
#                print("All common ancestor of the population at time ", self.back_time)
                self.all_common_anc = self.back_time
#            print("min", min_nb_desc, "max", max_nb_desc)
        '''


    def check_fused_segments(self):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        INDIVID, CHRNO, CHRINDEX, A, B = 0, 1, 2, 3, 4
        
        #self.test_segments.sort()
        self.test_segments.view('i4,i4,i4,i4,i4').sort(order=['f0'], axis=0)

        #for i in range(len(self.test_segments)):
        #    print(self.test_segments[i])
        #print(self.test_segments)
        #sys.exit()             


        indices_to_delete = set()
        for i in range(1, len(self.test_segments)):
            
            #if seg[i] is right next to seg[i-1], on the same indiv and chrsm, fuse it --> index i - 1 will get deleted later on
            if (self.test_segments[i, INDIVID] == self.test_segments[i - 1, INDIVID] and 
                self.test_segments[i, CHRNO] == self.test_segments[i - 1, CHRNO] and 
                self.test_segments[i, CHRINDEX] == self.test_segments[i - 1, CHRINDEX] and 
                self.test_segments[i, A] <= self.test_segments[i - 1, B]+1):
                    self.test_segments[i, A] = self.test_segments[i - 1, A]
                    self.test_segments[i, B] = max(self.test_segments[i - 1, B],self.test_segments[i, B])
                    indices_to_delete.add(i - 1)
        
        
        #print(f"nb_segments={len(self.test_segments)}    nb_fusions={len(indices_to_delete)}")
        
        tmp_segments = np.empty( (len(self.test_segments), 5), np.int32 )
        tmp_cpt = 0
        #it is apparently faster to delete a list of indices by recreating the array but omitting deleted indices
        for i in range(len(self.test_segments)):
            if i not in indices_to_delete:
                tmp_segments[tmp_cpt] = self.test_segments[i]
                tmp_cpt += 1

        
        self.test_segments = np.resize(tmp_segments, (tmp_cpt, 5))

        self.test_segments.view('i4,i4,i4,i4,i4').sort(order=['f0'], axis=0)

        #for i in range(len(self.test_segments)):
        #    print(self.test_segments[i])
        #print(self.test_segments)
        #sys.exit()

        #self.check_fused_segments_legacy()


        global DEBUG
        if DEBUG:
            self.check_fused_segments_legacy()
            tmp_segments = []
            for nc in range(nb_chrsm):
                for seg in self.cur_segments[nc]:
                    tmp_segments.append( (seg.indiv_id, seg.chrsm[0], int(seg.chrsm[1]), seg.a, seg.b) )
            tmp_segments.sort()

            self.test_segments.view('i4,i4,i4,i4,i4').sort(order=['f0'], axis=0)
        
            if len(self.test_segments) != len(tmp_segments):
                print(f"difflen: test_seg={len(self.test_segments)} vs tmp={len(tmp_segments)}")
                for i in range(100): #range(len(self.test_segments)):
                     issame = True
                     for j in range(5):
                         if self.test_segments[i][j] != tmp_segments[i][j]:
                             issame = False
                     if not issame or True:
                         print(f"test[{i}] = {self.test_segments[i]}   \nseg[{i}] =  {tmp_segments[i]}")
            else:
                for i in range(len(self.test_segments)):
                     issame = True
                     for j in range(5):
                         if self.test_segments[i][j] != tmp_segments[i][j]:
                              print(f"test[i] = {self.test_segments[i]}   \nseg[i] =  {tmp_segments[i]}")




    #side effect: sorts self.cur_segments
    def check_fused_segments_legacy(self):
        global pop_size, final_time, chrsm_len, recomb_rate, nb_chrsm, verbose, write_data, seed, chrsm_choice_rd, expected_nb_recomb, recomb_counter, recomb_pos_random
        tmp_segments = []
        genetic_mat = 0
        anc_ind_list = set()
        anc_chrsm_list = set()

        for nc in range(nb_chrsm):
            self.cur_segments[nc].sort(key=lambda s: (s.indiv_id, s.chrsm, s.a, s.b))  # sorts in place
            tmp_segments.append([])

            indices_to_delete = set()

            for i in range(1, len(self.cur_segments[nc])):

                #if seg[i] is right next to seg[i-1], on the same indiv and chrsm, fuse it --> index i - 1 will get deleted later on
                if (self.cur_segments[nc][i].indiv_id == self.cur_segments[nc][i - 1].indiv_id and 
                self.cur_segments[nc][i].chrsm == self.cur_segments[nc][i - 1].chrsm and 
                self.cur_segments[nc][i].a <= self.cur_segments[nc][i - 1].b+1):

                    self.cur_segments[nc][i].a = self.cur_segments[nc][i - 1].a
                    self.cur_segments[nc][i].b = max(self.cur_segments[nc][i - 1].b,self.cur_segments[nc][i].b)
                    indices_to_delete.add(i - 1)

            #it is apparently faster to delete a list of indices by recreating the array but omitting deleted indices
            for i in range(len(self.cur_segments[nc])):
                if i not in indices_to_delete:
                    tmp_segments[nc].append(self.cur_segments[nc][i])

        


        self.cur_segments = tmp_segments


        
    

'''
    def write_data(self,filename):
        d = {"nb_ind_genealogical_ancestors" : list(self.nb_ind_genealogical_ancestors),
             "nb_ind_genetic_ancestors" : list(self.nb_ind_genetic_ancestors),
             "nb_chr_genetic_ancestors" : list(self.nb_chr_genetic_ancestors),
             "nb_segments" : list(self.nb_segments),
             "genetic_mat" : list(self.genetic_mat),
             "first_commom_anc" : self.first_common_anc,
             "all_common_anc" : self.all_common_anc 
             }
        json_object = json.dumps(d)
        with open(filename, "w") as outfile:
            outfile.write(json_object)

    def write_all_distrib(self, filename):

        d = {"indiv_id" : [seg.indiv_id for seg in self.cur_segments],
             "chrsm" :    [seg.chrsm for seg in self.cur_segments],
             "a" :        [seg.a for seg in self.cur_segments],
             "b" :        [seg.b for seg in self.cur_segments],
             "time" : self.back_time}
        json_object = json.dumps(d)
        with open(filename, "w") as outfile:
            outfile.write(json_object)

'''



DEBUG = False

pop = Population()
for i in range(pop_size):
    pop.individuals[i] = Individual(i)
    # We initialize the whole population anyway to follow the genealogy and get the
    # genealogical ancestors
lin = Lineage(pop)

while (lin.back_time < final_time-1):
    if (verbose) & (lin.back_time%1==0):
       print(f"T={lin.back_time}  S={len(lin.test_segments)}")
       #print(f"T={lin.back_time}  S={sum(len(lc) for lc in lin.cur_segments)} "\
       #      f" Ngeneal={lin.nb_ind_genealogical_ancestors[lin.back_time]}  Ngeneti={lin.nb_ind_genetic_ancestors[lin.back_time]}"\
       #      f" Nb_bases={lin.genetic_mat[lin.back_time]}")

    lin.backward_step()
    lin.check_fused_segments()







