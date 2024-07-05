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