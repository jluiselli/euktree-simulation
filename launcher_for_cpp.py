import os as os
import argparse
from execo_engine import sweep
import json as json



parser = argparse.ArgumentParser()
parser.add_argument("--bin", help="compiled program to run")
parser.add_argument("--outdir", help="output directory to write the data")
parser.add_argument("--params", help="json file with the combinations of parameters to run")
args = parser.parse_args()


with open(args.params, mode="r") as json_data:
    parameters = json.load(json_data)
sweeps = sweep(parameters)
print("The following combinations will be run")
print(sweeps)


def launch(d, bin_com):
    if "nb_recomb" in d:
        d["recomb_rate"] = d["nb_recomb"]/d["chrlen"]
        del d["nb_recomb"]
    if "total_len" in d:
        d["chrlen"] = int(d["total_len"]/d["nbchr"])
        del d["total_len"]

    data_file = "nbchr-"+str(d["nbchr"])+"-chrlen-"+str(d["chrlen"])+"-nb_gen-"+str(d["nb_gen"])+"-pop_size-"
    data_file += str(d["pop_size"])+"-recomb_rate-"+str(d["recomb_rate"])+"-seed-"+str(d["seed"])
    exact_ghosts = 0 if d["pop_size"] > 4000 else 1
    data_file += "-exact_ghosts-" + str(exact_ghosts) +".csv"

    launching_command = bin_com
    for key in d.keys():
        launching_command = launching_command + " --" + key + " " + str(d[key])

    print(launching_command)
    print(data_file)
    if os.path.exists(data_file):
        print("already computed, skipping combination")
    else:
        os.system(launching_command)


bin_com = os.path.abspath(args.bin)
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)
os.chdir(args.outdir)

for d in sweeps:
    launch(d, bin_com)