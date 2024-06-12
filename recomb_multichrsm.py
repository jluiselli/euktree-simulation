import numpy as np
import random as random
import json
import os


#apparently the fastest way to rand (faster would be to do them in batch)
def my_randint(m1, m2):
    return int((m2 + 1 - m1) * random.random() + m1)




class Individual:
    def __init__(self, id, parent1, parent2, chrsm_len, recomb_rate, child=None, nb_chrsm=1):
        if child != None:
            self.final_desc = child.final_desc
            self.nb_chrsm = child.nb_chrsm
        else:
            self.final_desc = [id]
            self.nb_chrsm = nb_chrsm

        # In case of multi chromosomes, we assume they all have the same length and the same average
        # recomb_rate per chromosomes
        nb_recombs = np.random.binomial(chrsm_len, recomb_rate, self.nb_chrsm)
        # Positions de recombinaisons
        self.recomb_par1 = [np.random.choice(chrsm_len-1, nb_recombs[i], replace=False)+1
                            for i in range(self.nb_chrsm)]

        nb_recombs = np.random.binomial(chrsm_len, recomb_rate, self.nb_chrsm)
        self.recomb_par2 = [np.random.choice(chrsm_len-1, nb_recombs[i], replace=False)+1
                            for i in range(self.nb_chrsm)]

        self.par1_gave_chrsm = np.random.randint(0,2, self.nb_chrsm)
        self.par2_gave_chrsm = np.random.randint(0,2, self.nb_chrsm)

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

        if parent1 != None:
            self.parent1_id = parent1.id
            self.parent2_id = parent2.id

        self.id = id # position dans la population


class Population:
    def __init__(self, size, chrsm_len, recomb_rate, nb_chrsm):
        self.size = size
        self.individuals = {}
        for i in range(size):
            self.individuals[i] = Individual(i, None, None, chrsm_len, recomb_rate, None, nb_chrsm)
            # We initialize the whole population anyway to follow the genealogy and get the
            # genealogical ancestors
        self.chrsm_len = chrsm_len
        self.r = recomb_rate
        self.nb_chrsm = nb_chrsm

    #this is the previous forward_step, now it generates the parent indivs from the child indivs
    #side effect: child_pop indivs get parent1 and parent2 assigned
    def generate_from_child_population(self, child_pop, to_check):

        for indiv in child_pop.individuals.values():
            par1_id = my_randint(0,self.size-1)
            par2_id = my_randint(0,self.size-1)
            while par1_id == par2_id:
                par2_id = my_randint(0, self.size-1)

            if par1_id not in self.individuals:
                self.individuals[par1_id] = Individual(par1_id, None, None, self.chrsm_len, self.r, indiv)
            elif to_check:
                self.individuals[par1_id].final_desc[0:0] = indiv.final_desc

            if par2_id not in self.individuals:
                self.individuals[par2_id] = Individual(par2_id, None, None, self.chrsm_len, self.r, indiv)
            elif to_check:
                self.individuals[par2_id].final_desc[0:0] = indiv.final_desc

            indiv.parent1_id = par1_id
            indiv.parent2_id = par2_id


class Segment:
    def __init__(self,id,chrsm,a,b):
        if chrsm == 0:
            chrsm = (0,0)
        elif chrsm == 1:
            chrsm = (0,1)
        elif type(chrsm) != tuple:
            print("error in chrsm type, should be (chrsm number, 0 or 1) as we hav multiple chrsms")
            exit(1)
        self.indiv_id = id
        self.chrsm = chrsm
        self.a = a
        self.b = b

    def print(self):
        print(self.indiv_id, self.chrsm, self.a, self.b)


class Lineage:
    def __init__(self, initial_pop, seg_list=None):
        self.pop = initial_pop
        self.back_time = 0

        if seg_list != None:
            self.cur_segments = seg_list # [ Segment(indiv_id, chrsm, a, b)]
        else: # We follow the whole population
            self.cur_segments = []
            for nc in range(self.pop.nb_chrsm):
                self.cur_segments.append([])
                for indiv_id in range(self.pop.size):
                    self.cur_segments[nc].append(Segment(indiv_id, (nc,0), 0, self.pop.chrsm_len-1))
                    self.cur_segments[nc].append(Segment(indiv_id, (nc,1), 0, self.pop.chrsm_len-1))

        self.nb_segments = [sum(len(seg_nc) for seg_nc in self.cur_segments)]
        anc_ind_list = []
        anc_chrsm_list = []
        nb_bases = 0
        for nc in range(self.pop.nb_chrsm):
            for seg in self.cur_segments[nc]:
                nb_bases += seg.b - seg.a +1
                if not seg.indiv_id in anc_ind_list:
                    anc_ind_list.append(seg.indiv_id)
                if not (seg.indiv_id, seg.chrsm) in anc_chrsm_list:
                    anc_chrsm_list.append((seg.indiv_id, seg.chrsm))
        self.nb_ind_genetic_ancestors = [len(anc_ind_list)]
        self.nb_ind_genealogical_ancestors = [len(anc_ind_list)]
        self.nb_chr_genetic_ancestors = [len(anc_chrsm_list)]

        self.genetic_mat = [nb_bases]

        self.has_separated = False
        self.first_common_anc = 0
        self.all_common_anc = 0


    def backward_step(self):
        next_segments = []  #next going backwards in time
        next_pop = Population(self.pop.size, self.pop.chrsm_len, self.pop.r, self.pop.nb_chrsm)
        next_pop.generate_from_child_population(self.pop, self.all_common_anc == 0)
        # next pop has only indivs with a descendant

        #at this point, each indiv in self.pop must has parent1_id and parent2_id set
        for nc in range(self.pop.nb_chrsm):
            next_segments.append([])
            for segment in self.cur_segments[nc]:
                if segment.chrsm[1] == 0: # chrsm A : from parent 1
                    par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par1_gave_chrsm[nc]
                    par_id =  self.pop.individuals[segment.indiv_id].parent1_id
                    recomb_pos_list =  self.pop.individuals[segment.indiv_id].recomb_par1[nc]
                else: # chrsm B : from parent 2
                    par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par2_gave_chrsm[nc]
                    par_id =  self.pop.individuals[segment.indiv_id].parent2_id
                    recomb_pos_list =  self.pop.individuals[segment.indiv_id].recomb_par2[nc]

                if par_gave_chrsm == 0:
                    # Le parent a donné le chrsm de gauche après la recombinaison
                    # start_chrsm = 0
                    start_chrsm = int((len([pos for pos in recomb_pos_list if pos <= segment.a]) % 2 ) == 1)
                else:
                    # Le parent a donné le chrsm de droite après la recombinaison
                    # start_chrsm = 1
                    start_chrsm = int((len([pos for pos in recomb_pos_list if pos <= segment.a]) % 2 ) == 0)

                nb_recomb_in_seg = len([pos for pos in recomb_pos_list if (pos > segment.a and pos <= segment.b)])
                if nb_recomb_in_seg == 0:
                    # Pas de coupure du segment suivi
                    # le chromosome où on l’envoie dépend seulement du nombre de coupures *avant*
                    next_segments[nc].append(Segment(par_id, start_chrsm, segment.a, segment.b))
                else: # le segment est coup en x morceaux
                    self.has_separated = True
                    positions = [segment.a]
                    for pos in [pos for pos in recomb_pos_list if (pos > segment.a and pos <= segment.b)]:
                        positions.append(pos)
                    positions.append(segment.b+1)
                    positions.sort()
                    for i in range(nb_recomb_in_seg + 1):
                        if positions[i+1] < positions[i]:
                            print(par_id, nb_recomb_in_seg ,start_chrsm, positions, int((len([pos for pos in recomb_pos_list if pos <= segment.a]) % 2 ) == 1) )
                        next_segments[nc].append(Segment(par_id, start_chrsm, positions[i], positions[i+1]-1))
                        start_chrsm = abs(start_chrsm - 1)

        self.back_time += 1
        self.cur_segments = next_segments
        self.pop.individuals = next_pop.individuals

        if self.all_common_anc == 0:
            min_nb_desc, max_nb_desc = self.pop.size, 0
            for indiv in self.pop.individuals.values():
                indiv.final_desc = list(set(indiv.final_desc))
                nb_desc = len(indiv.final_desc)
                if nb_desc < min_nb_desc:
                    min_nb_desc = nb_desc
                if nb_desc > max_nb_desc:
                    max_nb_desc = nb_desc
            if self.first_common_anc == 0 and max_nb_desc == self.pop.size:
#                print("Most recent common ancestor of the population at time ", self.back_time)
                self.first_common_anc = self.back_time
            if min_nb_desc == self.pop.size:
#                print("All common ancestor of the population at time ", self.back_time)
                self.all_common_anc = self.back_time
#            print("min", min_nb_desc, "max", max_nb_desc)


    #side effect: sorts self.cur_segments
    def check_fused_segments(self):
        tmp_segments = []
        genetic_mat = 0
        anc_ind_list = []
        anc_chrsm_list = []

        for nc in range(self.pop.nb_chrsm):
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
                    genetic_mat += self.cur_segments[nc][i].b-self.cur_segments[nc][i].a+1
                    if not self.cur_segments[nc][i].indiv_id in anc_ind_list:
                        anc_ind_list.append(self.cur_segments[nc][i].indiv_id)
                    if not (self.cur_segments[nc][i].indiv_id, self.cur_segments[nc][i].chrsm) in anc_chrsm_list:
                        anc_chrsm_list.append((self.cur_segments[nc][i].indiv_id, self.cur_segments[nc][i].chrsm))
        self.nb_ind_genetic_ancestors.append(len(anc_ind_list))
        self.nb_chr_genetic_ancestors.append(len(anc_chrsm_list))


        if genetic_mat > self.genetic_mat[-1]:
            print(self.back_time, "we increased the genetic material we follow !", genetic_mat, self.genetic_mat)
            print("future segs")
            for seg in tmp_segments:
                seg.print()
            print("last segs")
            for seg in self.cur_segments:
                seg.print()

        self.cur_segments = tmp_segments

        self.nb_segments.append(sum(len(segnc) for segnc in self.cur_segments))
        self.nb_ind_genealogical_ancestors.append(len(self.pop.individuals))
        self.genetic_mat.append(int(genetic_mat))

        if genetic_mat < self.pop.chrsm_len:
            print(self.back_time, "genetic mat is ", genetic_mat, "while chrsm len", self.pop.chrsm_len)
            for seg in self.cur_segments:
                seg.print()   
    

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


