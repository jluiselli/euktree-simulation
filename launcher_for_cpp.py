import os as os
import argparse
from execo_engine import sweep


parameters = {
    "pop_size" : [500],
    "chrlen" : [10000],
    "nb_recomb" : [1],
    "seed" : [4334892,987321],
    "nbchr" : [3],
    "nb_gen" : [800],
}


sweeps = sweep(parameters)
print("The following combinations will be run")
print(sweeps)


parser = argparse.ArgumentParser()
parser.add_argument("--bin", help="compiled program to run")
parser.add_argument("--outdir", help="output directory to write the data")

args = parser.parse_args()


def launch(d, bin_com):
    d["recomb_rate"] = d["nb_recomb"]/d["chrlen"]
    del d["nb_recomb"]   

    launching_command = bin_com
    for key in d.keys():
        launching_command = launching_command + " --" + key + " " + str(d[key])

    print(launching_command)
    os.system(launching_command)


bin_com = os.path.abspath(args.bin)
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)
os.chdir(args.outdir)

for d in sweeps:
    launch(d, bin_com)