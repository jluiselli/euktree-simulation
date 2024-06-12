from recomb_multichrsm import *
import random as random
import os as os
import multiprocessing
from execo_engine import sweep


parameters = {
    "popsize" : [1000 ],
    "chrsm_len" : [5000],
    "recomb_rate" : [1/5000,2/5000],
    "seed" : [47098438,3149838],
    "run_id" : [0,1,2,3,4],
    "nb_chrsm" : [1,10,36],
    "final_time" : [100]
    "verbose" : [True]
}


sweeps = sweep(parameters)
print(sweeps)

def launch(d):
    sim = Simulation(d)
    sim.run()

pool = multiprocessing.Pool(4)
pool.map(launch, sweeps)
