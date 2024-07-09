import pandas as pd
import numpy as np
import getopt
import sys
import os
import glob

def main(argv):
    result_dir="defaultdir"
    output_csv="default.csv"
    last_generation=-1

    try:
        opts, extraparams = getopt.getopt(sys.argv[1:],"hi:o:g:",["input-dir=","output-csv="]) 
    except getopt.GetoptError:
        print('data2csv.py -i <inputdirectory> -o <outputcsv> -g <last_generation>')
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
        k = comb_dir[:-4].replace('/',' ').split('-')
        for no,elt in enumerate(k):
            if elt == '1e' or elt == '2e' or elt=='5e':
                k[no] = k[no] + '-' + k[no+1]
                del k[no+1]

        i = iter(k)
        params = dict(zip(i,i))
        

        # Verify stats files are computed
        if comb_dir[:4] != 'coal':
            print(params)
            List = glob.glob(result_dir+'/'+comb_dir)

            if len(List) > 0:
                if os.path.isfile(List[0]):
                    
                    # returns JSON object as a dictionary
                    data = pd.read_csv(List[0])


                    df = pd.DataFrame.from_dict(data)
                    for p in params:
                        df[p] = params[p]
                    total_df = pd.concat([total_df,df])
                    print(total_df.shape)
            else:
                print("empty")
        else:
            print("not doing coal data")
    
    
    total_df.to_csv(output_csv,sep=';')



if __name__ == "__main__":
   main(sys.argv[1:])
