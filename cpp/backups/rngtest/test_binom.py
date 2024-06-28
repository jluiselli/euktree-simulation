import numpy as np
import random as random
#import json
import os
import sys

seed = 47098439


np.random.seed(seed)

nbtrials = 500000
prob_success = 1/nbtrials
size = 200000000
nb_runs = 10


for i in range(nb_runs):
	arr = np.random.binomial(nbtrials, prob_success, size)
	print(f"run {i + 1}/{nb_runs} done")


