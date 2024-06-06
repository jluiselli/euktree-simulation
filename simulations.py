from recomb_model_v2 import *
import random as random
import os as os
import multiprocessing
from execo_engine import sweep

final_time = 800

parameters = {
    "popsize" : [20, 100, 1000, 4000],
    "chrsm_len" : [10000, 100000, 1000000],
    "recomb_rate" : [1,10,36],
    "run_id" : [10,11,12,13,14]
}


sweeps = sweep(parameters)
print(sweeps)



def simulation(d):
    popsize = d["popsize"]
    chrsm_len = d["chrsm_len"]
    recomb_rate = d["recomb_rate"]
    run_id = d["run_id"]

    filename = "ghosts/pop_size-"+str(popsize)+"-chrsm_len-"+str(chrsm_len)+"-recomb-"+str(recomb_rate)+\
    "-rep-"+str(run_id)+".json"
    print(filename)
    if not os.path.isfile(filename):
        pop = Population(popsize, chrsm_len, recomb_rate/chrsm_len)
        for i in range(pop.size):
            pop.individuals[i] = Individual(i, None, None, chrsm_len, pop.r, None)

        lin = Lineage(pop)
        # lin  = Lineage(pop, [Segment(0, random.randint(0,1), start_seg, start_seg+seg_len-1)])
        while (lin.back_time<final_time-1):
            if lin.back_time%100==0:
                print(f"T={lin.back_time}  S={len(lin.cur_segments)}   N={len(lin.pop.individuals)}  Nb_bases={lin.genetic_mat[-1]}")
            lin.backward_step()
            lin.check_fused_segments()
            if len(lin.cur_segments) == 1:
                print("COALESCENT at time ", lin.back_time)
                lin.write_data(filename)
                break
        lin.write_data(filename)
        #lin.write_all_distrib(filename2)
    else:
        print("already computed")


pool = multiprocessing.Pool(8)
pool.map(simulation, sweeps)
