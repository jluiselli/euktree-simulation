from recomb_multichrsm import *
import random as random
import os as os
import multiprocessing
from execo_engine import sweep


parameters = {
    "popsize" : [1000 ],
    "chrsm_len" : [5000],
    "recomb_rate" : [1,10,36],
    "run_id" : [0,1,2,3,4],
    "nb_chrsm" : [1,10,36],
    "final_time" : [100]
}


sweeps = sweep(parameters)
print(sweeps)



def simulation(d):
    popsize = d["popsize"]
    chrsm_len = d["chrsm_len"]
    recomb_rate = d["recomb_rate"]
    run_id = d["run_id"]
    nb_chrsm = d["nb_chrsm"]
    final_time = d["final_time"]

    filename = "multi_chrsm/nbchrsm-"+str(nb_chrsm)+"-pop_size-"+str(popsize)+"-chrsm_len-"+\
        str(chrsm_len)+"-recomb-"+str(recomb_rate)+"-final_time-"+str(final_time)+\
    "-rep-"+str(run_id)+".json"
    print(filename)
    if not os.path.isfile(filename):
        pop = Population(popsize, chrsm_len, recomb_rate/chrsm_len, nb_chrsm)
        # for i in range(pop.size):
        #     pop.individuals[i] = Individual(i, None, None, chrsm_len, pop.r, None)

        lin = Lineage(pop)
        # lin  = Lineage(pop, [Segment(0, random.randint(0,1), start_seg, start_seg+seg_len-1)])
        while (lin.back_time<final_time-1):
            if lin.back_time%50==0:
               print(f"T={lin.back_time}  S={len(lin.cur_segments)}   N={len(lin.pop.individuals)}  Nb_bases={lin.genetic_mat[-1]}")
            lin.backward_step()
            lin.check_fused_segments()
            if sum(len(lc) for lc in lin.cur_segments) == 1:
                print("COALESCENT at time ", lin.back_time)
                lin.write_data(filename)
                break
        lin.write_data(filename)
        #lin.write_all_distrib(filename2)
    else:
        print("already computed")


pool = multiprocessing.Pool(4)
pool.map(simulation, sweeps)
