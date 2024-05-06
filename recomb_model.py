import numpy as np
import random as random
import json

class Individual:
    def __init__(self, id, parent1, parent2, chrsm_len):
        # Positions de recombinaisons
        self.recomb_par1 = random.randint(1,chrsm_len)
        self.recomb_par2 = random.randint(1,chrsm_len)

        self.par1_gave_chrsm = random.randint(0,1)
        self.par2_gave_chrsm = random.randint(0,1)
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
    def __init__(self, size, final_time, chrsm_len):
        self.size = size
        self.T = 0
        self.simulation = np.empty((final_time+1, size), dtype=object)
        self.simulation[0] = [Individual(i, None, None, chrsm_len) for i in range(size)]
        self.chrsm_len = chrsm_len

    def forward_step(self):
        self.T+=1
        for i in range(self.size):
            par1_id = random.randint(0,self.size-1)
            par2_id = random.randint(0,self.size-1)
            while par1_id == par2_id:
                par2_id = random.randint(0,self.size-1)
            self.simulation[self.T][i] = Individual(i,self.simulation[self.T-1][par1_id],
                                                    self.simulation[self.T-1][par2_id],
                                                    self.chrsm_len)


class Segment:
    def __init__(self,id,chrsm,a,b):
        self.indiv_id = id
        self.chrsm = chrsm
        self.a = a
        self.b = b


class Lineage:
    def __init__(self,indiv_id, chrsm, a, b, pop):
        self.start_indiv = indiv_id
        self.start_chrsm = chrsm
        self.start_segment = [a,b]
        self.pop = pop
        self.end_time = pop.T
        self.back_time = 0
        self.lineage = [[Segment(indiv_id,chrsm,a,b)]]
        self.has_separated = False
        self.nb_segments = np.zeros(pop.T)
        self.nb_ancestors = np.zeros(pop.T)

    def backward_step(self):
        self.lineage += [[]]
        for segment in self.lineage[self.back_time]:
            # print(segment.indiv_id, segment.chrsm, segment.a, segment.b)
            if segment.chrsm == 0: # chrsm A : from parent 1
                par_gave_chrsm =  self.pop.simulation[self.end_time-self.back_time][segment.indiv_id].par1_gave_chrsm
                par_id =  self.pop.simulation[self.end_time-self.back_time][segment.indiv_id].parent1_id
                recomb_pos =  self.pop.simulation[self.end_time-self.back_time][segment.indiv_id].recomb_par1
            else: # chrsm B : from parent 2
                par_gave_chrsm =  self.pop.simulation[self.end_time-self.back_time][segment.indiv_id].par2_gave_chrsm
                par_id =  self.pop.simulation[self.end_time-self.back_time][segment.indiv_id].parent2_id
                recomb_pos =  self.pop.simulation[self.end_time-self.back_time][segment.indiv_id].recomb_par2
            if par_gave_chrsm == 0:
                # Le parent a donné le chrsm de gauche après la recombinaison
                if recomb_pos <= segment.a:
                    # La coupure était au-dessus du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "droite" d’avant
                    # a et b restent inchangés
                    self.lineage[self.back_time+1] += [Segment(par_id, 1, segment.a, segment.b)]
                elif recomb_pos >= segment.b:
                    # La coupure était en-dessous du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "gauche"
                    # a et b restent inchangés
                    self.lineage[self.back_time+1] += [Segment(par_id, 0, segment.a, segment.b)]
                else: # Coupure dans le segment d’intérêt : on doit suivre 2 segments
                    # De a à la coupure, ça vient de "gauche"
                    self.lineage[self.back_time+1] += [Segment(par_id, 0, segment.a, recomb_pos)]
                    # De la coupure à b, ça vient de "droite"
                    self.lineage[self.back_time+1] += [Segment(par_id, 1, recomb_pos, segment.b)]
                    self.has_separated = True
            else:
                # Le parent a donné le chrsm de droite après la recombinaison
                if recomb_pos <= segment.a:
                    # La coupure était au-dessus du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "gauche" d’avant
                    # a et b restent inchangés
                    self.lineage[self.back_time+1] += [Segment(par_id, 0, segment.a, segment.b)]
                elif recomb_pos >= segment.b:
                    # La coupure était en-dessous du segment qu’on suit
                    # Le chromosome qu’on va suivre est donc celui de "droite"
                    # a et b restent inchangés
                    self.lineage[self.back_time+1] += [Segment(par_id, 1, segment.a, segment.b)]
                else: # Coupure dans le segment d’intérêt : on doit suivre 2 segments
                    # De a à la coupure, ça vient de "droite"
                    self.lineage[self.back_time+1] += [Segment(par_id, 1, segment.a, recomb_pos)]
                    # De la coupure à b, ça vient de "gauche"
                    self.lineage[self.back_time+1] += [Segment(par_id, 0, recomb_pos, segment.b)]
                    self.has_separated = True

        self.back_time += 1
        # print(len(self.lineage[self.back_time]))
    

    def check_fused_segments(self):
        # We need to fuse if we are following several times the same segment
        for id in range(self.pop.size):
            id_is_ancestor = False
            for chrsm in [0,1]:
                segments = []
                positions = []
                for seg in self.lineage[self.back_time]:
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
                    self.lineage[self.back_time].remove(seg)
                for tup in positions:
                    self.lineage[self.back_time].append(Segment(id, chrsm, tup[0], tup[1]))
                if len(segments) != 0:
                    id_is_ancestor = True
            if id_is_ancestor:
                self.nb_ancestors[self.back_time] += 1
        self.nb_segments[self.back_time] = len(self.lineage[self.back_time])
    
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