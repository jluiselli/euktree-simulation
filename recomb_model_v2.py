import numpy as np
import random as random
import json


#apparently the fastest way to rand (faster would be to do them in batch)
def my_randint(m1, m2):
    return int((m2 + 1 - m1) * random.random() + m1)




class Individual:
    def __init__(self, id, parent1, parent2, chrsm_len, recomb_rate, child):
        nb_recomb = np.random.binomial(chrsm_len, recomb_rate)
        # Positions de recombinaisons
        self.recomb_par1 = list(set([my_randint(1,chrsm_len-1) for _ in range(nb_recomb)]))     #random.randint(1,chrsm_len)
        self.recomb_par2 = list(set([my_randint(1,chrsm_len-1) for _ in range(nb_recomb)]))     #random.randint(1,chrsm_len)
        
        self.par1_gave_chrsm = my_randint(0, 1) #random.randint(0,1)
        self.par2_gave_chrsm = my_randint(0, 1) #random.randint(0,1)
        if child != None:
            self.final_desc = child.final_desc
        else:
            self.final_desc = [id]
        
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
    def __init__(self, size, chrsm_len, recomb_rate):
        self.size = size
        self.individuals = {}   #dict that only contains active indivs, key = indiv_id, value = Individual object
        self.chrsm_len = chrsm_len
        self.r = recomb_rate

    #this is the previous forward_step, now it generates the parent indivs from the child indivs
    #side effect: child_pop indivs get parent1 and parent2 assigned
    def generate_from_child_population(self, child_pop, to_check):
        
        for indiv in child_pop.individuals.values():
            par1_id = my_randint(0,self.size-1)             #random.randint(0,self.size-1)
            par2_id = my_randint(0,self.size-1)                                #random.randint(0,self.size-1)
            while par1_id == par2_id:
                par2_id = my_randint(0, self.size-1)        #random.randint(0,self.size-1)
                
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
            for indiv_id in range(self.pop.size):
                self.cur_segments.append(Segment(indiv_id, 0, 0, self.pop.chrsm_len-1))
                self.cur_segments.append(Segment(indiv_id, 1, 0, self.pop.chrsm_len-1))

        self.nb_segments = [len(self.cur_segments)]
        anc_list = []
        nb_bases = 0
        for seg in self.cur_segments:
            nb_bases += seg.b - seg.a +1
            if not seg.indiv_id in anc_list:
                anc_list.append(seg.indiv_id)
        self.nb_ancestors = [len(anc_list)]
        self.genetic_mat = [nb_bases]

        self.has_separated = False
        self.first_common_anc = 0
        self.all_common_anc = 0


    def backward_step(self):
        next_segments = []  #next going backwards in time
        next_pop = Population(self.pop.size, self.pop.chrsm_len, self.pop.r)
        next_pop.generate_from_child_population(self.pop, self.all_common_anc == 0)
        
        #print(self.pop.individuals.keys())
        #at this point, each indiv in self.pop must has parent1_id and parent2_id set
        for segment in self.cur_segments:
            
            if segment.chrsm == 0: # chrsm A : from parent 1
                par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par1_gave_chrsm
                par_id =  self.pop.individuals[segment.indiv_id].parent1_id
                recomb_pos_list =  self.pop.individuals[segment.indiv_id].recomb_par1
            else: # chrsm B : from parent 2
                par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par2_gave_chrsm
                par_id =  self.pop.individuals[segment.indiv_id].parent2_id
                recomb_pos_list =  self.pop.individuals[segment.indiv_id].recomb_par2
                
            
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
                next_segments.append(Segment(par_id, start_chrsm, segment.a, segment.b))
            else: # le segment est coup en x morceaux
                self.has_separated = True
                positions = [segment.a]
                for pos in [pos for pos in recomb_pos_list if (pos > segment.a and pos <= segment.b)]:
                    positions.append(pos)
                positions.append(segment.b+1)
                positions.sort()
                for i in range(nb_recomb_in_seg + 1):
                    next_segments.append(Segment(par_id, start_chrsm, positions[i], positions[i+1]-1))
                    start_chrsm = abs(start_chrsm - 1)
            
        self.back_time += 1
        
        
        self.cur_segments = next_segments
        
        self.pop.individuals = {}   #clear pop
        #next pop has only indivs with a segment - you could argue they should not have been created in the first place...
        for s in self.cur_segments:
            if s.indiv_id not in self.pop.individuals:
                self.pop.individuals[s.indiv_id] = next_pop.individuals[s.indiv_id]
            
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
                print("Most recent common ancestor of the population at time ", self.back_time)
                self.first_common_anc = self.back_time
            if min_nb_desc == self.pop.size:
                print("All common ancestor of the population at time ", self.back_time)
                self.all_common_anc = self.back_time
            print("min", min_nb_desc, "max", max_nb_desc)

    
    #side effect: sorts self.cur_segments
    def check_fused_segments(self):
        
        self.cur_segments.sort(key=lambda s: (s.indiv_id, s.chrsm, s.a, s.b))  # sorts in place
        
        indices_to_delete = set()
        
        for i in range(1, len(self.cur_segments)):
            
            #if seg[i] is right next to seg[i-1], on the same indiv and chrsm, fuse it --> index i - 1 will get deleted later on
            if (self.cur_segments[i].indiv_id == self.cur_segments[i - 1].indiv_id and 
               self.cur_segments[i].chrsm == self.cur_segments[i - 1].chrsm and 
               self.cur_segments[i].a <= self.cur_segments[i - 1].b+1):
               
                self.cur_segments[i].a = self.cur_segments[i - 1].a
                self.cur_segments[i].b = max(self.cur_segments[i - 1].b,self.cur_segments[i].b)
                indices_to_delete.add(i - 1)
        
        #it is apparently faster to delete a list of indices by recreating the array but omitting deleted indices
        tmp_segments = []
        genetic_mat = 0
        for i in range(len(self.cur_segments)):
            if i not in indices_to_delete:
                tmp_segments.append(self.cur_segments[i])
                genetic_mat += self.cur_segments[i].b-self.cur_segments[i].a+1

        self.cur_segments = tmp_segments

        self.nb_segments.append(len(self.cur_segments))
        self.nb_ancestors.append(len(self.pop.individuals))
        self.genetic_mat.append(genetic_mat)

        if genetic_mat < self.pop.chrsm_len:
            for seg in self.cur_segments:
                seg.print()
        
        
    
    
    
    #TODO: make this faster if things are too slow, should be feasible in a single loop over the segments
    def check_fused_segments_old(self):
        # We need to fuse if we are following several times the same segment
        for indiv in self.pop.individuals.values():
            id = indiv.id
            id_is_ancestor = False
            for chrsm in [0,1]:
                segments = []
                positions = []
                for seg in self.cur_segments:
                    if seg.indiv_id == id and seg.chrsm == chrsm:
                        segments += [seg]
                        positions += [(seg.a, seg.b)]
                if len(segments) > 1:
                    positions.sort(key=lambda tup: (tup[0], tup[1]))  # sorts in place
                    i = 1
                    while i < len(positions):
                        if positions[i][0] <= positions[i-1][1]:
                            # print(positions)
                            if positions[i][1] <= positions[i-1][1]:
                                # le 2e segment est inutile
                                positions.pop(i)
                            else:
                                new_a = positions[i-1][0]
                                new_b = positions[i][1]
                                # il faut fusionner les 2 segments
                                positions = positions[:i-1] + [(new_a,new_b)] + positions[i+1:]
                        else:
                            i+=1
                for seg in segments:
                    self.cur_segments.remove(seg)
                for tup in positions:
                    self.cur_segments.append(Segment(id, chrsm, tup[0], tup[1]))
                if len(segments) != 0:
                    id_is_ancestor = True
            #if id_is_ancestor:
            #    self.nb_ancestors[self.back_time] += 1
        #self.nb_segments[self.back_time] = len(self.lineage[self.back_time])
    
    
    
    def write_data(self,filename):
        d = {"nb_ancestors" : list(self.nb_ancestors),
             "nb_segments" : list(self.nb_segments),
             "genetic_mat" : list(self.genetic_mat)}
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
            
            
chrsm_len = 2000
final_time = 100000

pop = Population(1024, chrsm_len)
pop.individuals[0] = Individual(0, None, None, chrsm_len)   #now pop indvid don't exist until needed

lin = Lineage(indiv_id = 0, chrsm = 0, a = 0, b = 2000, initial_pop = pop)
while (lin.back_time < final_time):
    #and
    #((not lin.has_separated) or (lin.has_separated and len(lin.lineage[lin.back_time])>1))):
    
    if lin.back_time%100==0:
        print(f"T={lin.back_time}  S={len(lin.cur_segments)}   N={len(lin.pop.individuals)}")
    lin.backward_step()
    lin.check_fused_segments()


