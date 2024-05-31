import numpy as np
import random as random
import json


#apparently the fastest way to rand (faster would be to do them in batch)
def my_randint(m1, m2):
    return int((m2 + 1 - m1) * random.random() + m1)




class Individual:
    def __init__(self, id, parent1, parent2, chrsm_len):
        # Positions de recombinaisons
        self.recomb_par1 = my_randint(1,chrsm_len)     #random.randint(1,chrsm_len)
        self.recomb_par2 = my_randint(1,chrsm_len)     #random.randint(1,chrsm_len)
        
        self.par1_gave_chrsm = my_randint(0, 1) #random.randint(0,1)
        self.par2_gave_chrsm = my_randint(0, 1) #random.randint(0,1)
        
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
    def __init__(self, size, chrsm_len):
        self.size = size
        self.individuals = {}   #dict that only contains active indivs, key = indiv_id, value = Individual object
        self.chrsm_len = chrsm_len

    #this is the previous forward_step, now it generates the parent indivs from the child indivs
    #side effect: child_pop indivs get parent1 and parent2 assigned
    def generate_from_child_population(self, child_pop):
        
        for indiv in child_pop.individuals.values():
            par1_id = my_randint(1,self.size-1)             #random.randint(0,self.size-1)
            par2_id = my_randint(1,self.size-1)                                #random.randint(0,self.size-1)
            while par1_id == par2_id:
                par2_id = my_randint(1, self.size-1)        #random.randint(0,self.size-1)
                
            if par1_id not in self.individuals:
                self.individuals[par1_id] = Individual(par1_id, None, None, self.chrsm_len)
            if par2_id not in self.individuals:
                self.individuals[par2_id] = Individual(par2_id, None, None, self.chrsm_len)
            
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
    def __init__(self,indiv_id, chrsm, a, b, initial_pop):
        self.start_indiv = indiv_id 
        self.start_chrsm = chrsm    #useful?
        
        self.pop = initial_pop
        self.back_time = 0
        
        self.cur_segments = [ Segment(indiv_id, chrsm, a, b)  ]
        self.has_separated = False
        
        #self.nb_segments = np.zeros(pop.T)
        #self.nb_ancestors = np.zeros(pop.T)

    def backward_step(self):
        next_segments = []  #next going backwards in time
        next_pop = Population(self.pop.size, self.pop.chrsm_len)
        next_pop.generate_from_child_population(self.pop)
        
        #print(self.pop.individuals.keys())
        #at this point, each indiv in self.pop must has parent1_id and parent2_id set
        for segment in self.cur_segments:
            
            if segment.chrsm == 0: # chrsm A : from parent 1
                par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par1_gave_chrsm
                par_id =  self.pop.individuals[segment.indiv_id].parent1_id
                recomb_pos =  self.pop.individuals[segment.indiv_id].recomb_par1
            else: # chrsm B : from parent 2
                par_gave_chrsm =  self.pop.individuals[segment.indiv_id].par2_gave_chrsm
                par_id =  self.pop.individuals[segment.indiv_id].parent2_id
                recomb_pos =  self.pop.individuals[segment.indiv_id].recomb_par2
                
            
            if par_gave_chrsm == 0:
                # Le parent a donné le chrsm de gauche après la recombinaison
                if recomb_pos <= segment.a:
                    # La coupure était au-dessus du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "droite" d’avant
                    # a et b restent inchangés
                    next_segments.append( Segment(par_id, 1, segment.a, segment.b) )
                elif recomb_pos > segment.b:
                    # La coupure était en-dessous du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "gauche"
                    # a et b restent inchangés
                    next_segments.append( Segment(par_id, 0, segment.a, segment.b) )
                else: # Coupure dans le segment d’intérêt : on doit suivre 2 segments
                    # De a à la coupure, ça vient de "gauche"
                    next_segments.append( Segment(par_id, 0, segment.a, recomb_pos-1) )
                    # De la coupure à b, ça vient de "droite"
                    next_segments.append( Segment(par_id, 1, recomb_pos, segment.b) )
                    self.has_separated = True
            else:
                # Le parent a donné le chrsm de droite après la recombinaison
                if recomb_pos <= segment.a:
                    # La coupure était au-dessus du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "gauche" d’avant
                    # a et b restent inchangés
                    next_segments.append( Segment(par_id, 0, segment.a, segment.b) )
                elif recomb_pos > segment.b:
                    # La coupure était en-dessous du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "droite"
                    # a et b restent inchangés
                    next_segments.append( Segment(par_id, 1, segment.a, segment.b) )
                else: # Coupure dans le segment d’intérêt : on doit suivre 2 segments
                    # De a à la coupure, ça vient de "droite"
                    next_segments.append( Segment(par_id, 1, segment.a, recomb_pos-1) )
                    # De la coupure à b, ça vient de "gauche"
                    next_segments.append( Segment(par_id, 0, recomb_pos, segment.b) )
                    self.has_separated = True

            
        self.back_time += 1
        
        
        self.cur_segments = next_segments
        
        self.pop.individuals = {}   #clear pop
        #next pop has only indivs with a segment - you could argue they should not have been created in the first place...
        for s in self.cur_segments:
            if s.indiv_id not in self.pop.individuals:
                self.pop.individuals[s.indiv_id] = next_pop.individuals[s.indiv_id]
            
    
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
                indices_to_delete.add(i - 1)
        
        #it is apparently faster to delete a list of indices by recreating the array but omitting deleted indices
        tmp_segments = []
        for i in range(len(self.cur_segments)):
            if i not in indices_to_delete:
                tmp_segments.append(self.cur_segments[i])
        self.cur_segments = tmp_segments
        
        
    
    
    
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
        d = {"population_size" : self.pop.size,
             "chrsm_length" : self.pop.chrsm_len,
             "seg_len" : self.start_segment[1] - self.start_segment[0],
             "start_indiv" : self.start_indiv,
             "start_chrsm" : self.start_chrsm,
             "start_position" : self.start_segment[0],
             "total_generations" : self.pop.T,
             "nb_ancestors" : list(self.nb_ancestors),
             "nb_segments" : list(self.nb_segments)}
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


