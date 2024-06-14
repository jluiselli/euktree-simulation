import pandas as pd
import numpy as np
import getopt
import sys
import os
import glob
import json

def main(argv):
    result_dir="defaultdir"
    output_csv="default.csv"

    try:
        opts, extraparams = getopt.getopt(sys.argv[1:],"hi:o:g:",["input-dir=","output-csv="]) 
    except getopt.GetoptError:
        print('data2csv.py -i <inputdirectory> -o <outputcsv>')
        sys.exit(2)

    for o,p in opts:
        if o == '-h':
            print('data2csv.py -i <inputdir> -o <outputcsv> -g <last_generation>')
            print('inputdir: a input directory where the results are stored')
            print('outputcsv: an output directory where the results will be stored')
            print('last_generation: check if all data are available up to last_generation')
            sys.exit(2)
        elif o in ['-i','--input-dir']:
            result_dir = p
        elif o in ['-o','--output-csv']:
            output_csv = p

    if result_dir == "defaultdir":
        print("You must define a input directory where the results are stored")
        sys.exit()

    if output_csv == "default.csv":
        print("You must define an output directory where the results will be stored")
        sys.exit()

    print('Input result directory:',result_dir)
    print('Output CSV directory: ',output_csv)

    list_comb = os.listdir(result_dir)


    total_df = pd.DataFrame()

    for comb_dir in list_comb:
        k = comb_dir.replace('/',' ').split('-')
        i = iter(k)
        params = dict(zip(i,i))
        try:
            params['verbose'] = params['verbose'][:-5] # "remove.json of last param"
        except:
            print(params)
            continue
        if params["recomb_rate"]== "0002":
            params["recomb_rate"] = 0.002
        elif params["recomb_rate"]=="00002":
            params["recomb_rate"] = 0.0002
        # print(params)

        # Verify stats files are computed
        List = glob.glob(result_dir+'/'+comb_dir)

        if len(List) > 0:
            if os.path.isfile(List[0]):
                # Opening JSON file
                f = open(List[0])
                
                # returns JSON object as a dictionary
                data = json.load(f)
                df = pd.DataFrame()
                df["nb_segments"] = data["nb_segments"]
                df["first_commom_anc"] = data["first_commom_anc"]
                df["all_common_anc"] = data["all_common_anc"]
                df["back_time"] = [i for i in range(len(df["nb_segments"]))]
                df["nb_ind_genealogical_ancestors"] = data["nb_ind_genealogical_ancestors"]
                df["nb_ind_genetic_ancestors"] = data["nb_ind_genetic_ancestors"]
                df["nb_chr_genetic_ancestors"] = data["nb_chr_genetic_ancestors"]
                df["genetic_mat"] = data["genetic_mat"]
                            
                # Closing file
                f.close()

                # df = pd.DataFrame.from_dict(data)
                for p in params:
                    df[p] = params[p]
                # print(df)
                total_df = pd.concat([total_df,df])
                print(total_df.shape)
        else:
            print("empty")
    
    total_df.to_csv(output_csv,sep=';')



if __name__ == "__main__":
   main(sys.argv[1:])
