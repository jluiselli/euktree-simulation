from recomb_model import *
import random as random
import os as os

final_time = 50000
population_sizes = [500, 200, 1000]
chrsm_lengths = [10000, 5000, 15000]
segment_lengths = [100, 20, 500]

for j in range(5): #We repeat each conditions 5 times
    for popsize in [1000]:
        for chrsm_len in [5000,15000]:
            if popsize==1000 and chrsm_len==15000:
                continue
            for seg_len in segment_lengths:
                pop = Population(popsize,final_time, chrsm_len)
                while pop.T < final_time:
                    if (pop.T%500)==0:
                        print("current forward time", pop.T)
                    pop.forward_step()

                for indiv_id in range(10): # We take 10 random individuals to start our reconstructions
                    filename = "data/pop_size-"+str(popsize)+"-chrsm_len-"+str(chrsm_len)+"-seg_len-"+str(seg_len)+"-rep-"+str(j*10+indiv_id)+".json"
                    print(filename)
                    if os.path.isfile(filename):
                        start_seg = random.randint(0, chrsm_len - seg_len)
                        lin = Lineage(indiv_id, random.randint(0,1), start_seg, start_seg+seg_len, pop)
                        while (lin.back_time<final_time-1 and
                            ((not lin.has_separated)
                            or (lin.has_separated and len(lin.lineage[lin.back_time])>1))):
                            lin.backward_step()
                            lin.check_fused_segments()
                        lin.write_data(filename)
                    else:
                        print("already computed")

