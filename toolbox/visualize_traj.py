
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import argparse
import os
import glob

# Map names to file patterns
FILE_PATTERNS = {
    "entanglement":  "traj_entanglement.dat", 
    "norm":  "traj_norm.dat",
    "density":  "traj_densities.dat",
}
 

def find_files(name: str, dic):
    if name not in dic:
        print(f"Unknown name '{name}'. Available: {', '.join(dic.keys())}")
        return
    
    pattern = dic[name]
   
    
    return pattern

    

def main():
    parser = argparse.ArgumentParser(description="¨Plot specific quantity (works only for one trajectory)")

    parser.add_argument(
        "name", 
        type=str,
        help="Name of the quantity to visualize"
    )
    parser.add_argument(
        "-n", "--number", 
        type=int,
        default=None,
        help="Argument for a specific size or site"
    )

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    path_data = '/../build/res/'
    path = script_dir + path_data + find_files(args.name,FILE_PATTERNS)
    
    if(args.name == "entanglement" or args.name == "density"):
        quant_df = pd.read_csv( path, sep=',',skipinitialspace=True,header=0 )   
        quant_df.columns = quant_df.columns.str.strip() 

        time = np.loadtxt(script_dir + path_data + "traj_time.dat") 
        quant = quant_df['S_l='+str(args.number)]
 
    else:
        quant = np.loadtxt(path)  
        time = np.loadtxt(script_dir + path_data + "traj_time.dat")  

    fig, axs = plt.subplots( 1,1,  figsize=(5,3), dpi = 200, sharex=True)  

    ax = axs
    ax.plot(time,quant) 
    ax.set_xlabel("t")
    ylabel=""
    if args.name == "entanglement" :
        ylabel = ("S_"+ str(args.number) + "(t)")if args.number >0 else "S(t)"
    if args.name == "density" :
        ylabel = ("n_"+ str(args.number) + "(t)")if args.number >0 else "n(t)"
    if args.name == "norm" :
        ylabel = "N(t)" 
    
    ax.set_ylabel(ylabel)


    fig.suptitle("Quantity :  " + args.name)
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()





# Clean column names (strip whitespace)



# sent = pd.read_csv("SimQJ_project/build/res/traj_entanglement.dat")['       S_l=2']
# time = np.loadtxt("SimQJ_project/build/res/traj_time.dat") 
# # obs = np.loadtxt("SimQJ_project/build/res/traj_densities.dat")['       n_i=0']
                 
# # sent =  np.loadtxt('test/traj_0_ent_array_2.txt')
# # time = np.loadtxt("test/traj_0_time_array.txt") 
# # obs = np.loadtxt("test/traj_0_obs.txt") 

# fig, axs = plt.subplots( 1,2,  figsize=(8,3), dpi = 200, sharex=True)  

# ax = axs[0]
# ax.plot(time,sent) 



# ax = axs[1]
# ax.plot(time,obs_df["n_i=0"]) 
