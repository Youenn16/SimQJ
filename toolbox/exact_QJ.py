from scipy.linalg import expm 
from tqdm import tqdm
from tqdm import trange 
import numpy as np
import random as random 
import matplotlib.pyplot as plt 
import itertools
from functools import reduce
from bisect import bisect
import pandas as pd 
import argparse

 
Z_m = np.array([[1,0],[0,-1]],dtype=np.complex64)
X_m = np.array([[0,1],[1,0]],dtype=np.complex64)
Y_m = np.array([[0,-1j],[1j,0]],dtype=np.complex64)
I_m = np.eye(2,dtype=np.complex64) 

def C_op_dagger(x_pos, L):
    return  reduce(np.kron,  [Z_m]*(x_pos)+ [(X_m-1j*Y_m)/2]+[I_m]*(L-1-x_pos))
def C_op(x_pos, L):
    return  reduce(np.kron, [Z_m]*(x_pos) +[(X_m+1j*Y_m)/2]+[I_m]*(L-1-x_pos))    


def get_exp_ED(H,dt,N_prec): 
    expH_list = []  
    for i in range(N_prec):
        expH_list.append(expm(-1j*H*dt/(2**i)))
    return expH_list

def get_index(delta_proba_array,sum_delta_p, r_rand):
    if sum_delta_p > 0 :
        cum_delta_p = np.cumsum(delta_proba_array/sum_delta_p)
    x_pos = np.searchsorted(cum_delta_p, r_rand, side='left')
    return x_pos

def observable_value(Obs_operator, psi):
    return  (psi.T.conj()@Obs_operator@psi)[0][0].real  


def entropy_wavefunction(psi,l):  
    rho = np.reshape(psi,(2**(l),-1))
    rho = rho@rho.T.conj() 
    eigs = np.linalg.eigvalsh(rho)
    eigs = eigs[ eigs>0]
    return -np.dot(np.log(eigs),eigs).real

def get_init_binary_state(binary_psi0,N):
    psi = np.zeros((N,1),dtype=np.complex64) 
    psi[int(binary_psi0,2),:] = 1   
    return psi/np.linalg.norm(psi)



def get_H_ED_Ising(gamma,h,kappa, L, pbc=False):
    N = 2**L # Hilbert space size 
    
    H_ed = np.zeros((N,N),dtype=complex)
    c_xx = (1+kappa)/2
    c_yy = (1-kappa)/2
    J=1

    # Tunnelling terms
    for i in range(0,L-1):
        H_ed += -J*C_op_dagger(i+1,L)@C_op(i,L)
        H_ed += -J*C_op_dagger(i,L)@C_op(i+1,L) # -J c_{m+1}^\dagger c_{m} -J c_{m}^\dagger c_{m+1} 
        H_ed += -kappa*C_op_dagger(i,L)@C_op_dagger(i+1,L)
        H_ed += -kappa*C_op(i+1,L)@C_op(i,L)


    for i in range(0,L):
        H_ed += -h*C_op_dagger(i,L)@C_op(i,L) 
        H_ed += -1j*(gamma/2.0)*C_op_dagger(i,L)@C_op(i,L) # -\frac{gamma}{2} c_{m}^\dagger c_{m}  

    # Periodic conditions
    if pbc:
        H_ed += -J*C_op_dagger(0,L)@C_op(L-1,L)-J*C_op_dagger(L-1,L)@C_op(0,L)
        H_ed += -kappa*C_op_dagger(L-1,L)@C_op_dagger(0,L)-kappa*C_op(0,L)@C_op(L-1,L)

    return H_ed



def trajectory_evolution_ED(gamma,h,kappa,L,l,dt,t_max,Nb_part,pb_r, pb_r_site, binary_psi0="", pbc = False):
    N = 2**L 
    print("------------------------------------------------------------------")
    print("gamma : " , gamma, ", L  : ", L , ", t_max : ", t_max , ", dt : ", int(dt*100)/100 )
    print("------------------------------------------------------------------")

    if binary_psi0=="":  
        psi = 0 
    else:
        psi = get_init_binary_state(binary_psi0,N) # Wavefunction
         
    H = get_H_ED_Ising(gamma,h,kappa,L,pbc=pbc) # Hamiltonian 

    Meas_op = []
    for x in range(L):
        Meas_op.append(C_op_dagger(x,L)@C_op(x,L))

    Obs_op = C_op_dagger(0,L)@C_op(0,L)

    proba_jump_array_r = pb_r
    proba_jump_site_array_r = pb_r_site 
    # Outputs
    time_array = []
    ent_array = [] 
    obs_array = [] 
    norm_array = [] 
    index_jump_array = []
    time_jump_array = []

    proba_jump_array = [] 
    proba_jump_site_array = [] 
    
    N_prec = 15
    expH =  get_exp_ED(H,dt,N_prec)

    
    t = 0 
    delta_proba_array = np.zeros((L))
    sum_delta_p = 0 
    jump_count = 0 
    error_sampling=1e-5
 
    t_count=0
    r = proba_jump_array_r[jump_count]
    r_site =proba_jump_site_array_r[jump_count] 


    while t < t_max : 
        psi_next = expH[0]@psi
        p_jump = np.linalg.norm(psi_next)**2
        index_prec = 0
        dt_eff = 0 

        if r > p_jump :  
            epsilon = np.abs(r-p_jump)
            index_prec += 1  
            while epsilon > error_sampling  and index_prec < N_prec: 
                if r < p_jump : #under shot
                    psi = psi_next
                    dt_eff += dt/(2**(index_prec-1))
                    p_jump = np.linalg.norm(psi_next)**2
                else:
                    index_prec += 1   
                psi_next = expH[index_prec]@psi
                p_jump = np.linalg.norm(psi_next)**2
                epsilon = np.abs(r-p_jump)


            t += dt_eff + dt/(2**(index_prec-1))
            
            obs_array.append(observable_value(Obs_op, psi/np.linalg.norm(psi)))
            ent_array.append(entropy_wavefunction(psi/np.linalg.norm(psi),l))
            time_array.append(t)
            norm_array.append(np.linalg.norm(psi))

            for j in range(L): 
                delta_proba_array[j] = observable_value(Meas_op[j], psi/np.linalg.norm(psi))

            sum_delta_p = np.sum(delta_proba_array) 
            jump_pos = get_index(delta_proba_array, sum_delta_p, r_site)
            
            psi = Meas_op[jump_pos]@psi 
            psi = psi /np.linalg.norm(psi)  

            time_jump_array.append(t)
            norm_array.append(np.linalg.norm(psi))
            ent_array.append(entropy_wavefunction(psi/np.linalg.norm(psi),l))
            obs_array.append(observable_value(Obs_op, psi/np.linalg.norm(psi)))
            
            t += 0.001/(2**(N_prec))
            time_array.append(t)
 
            if jump_count + 1 < len(proba_jump_array_r):
                jump_count += 1 
            r = proba_jump_array_r[jump_count]
            r_site = proba_jump_site_array_r[jump_count]   

        else :  

            psi = psi_next 
            norm_array.append(np.linalg.norm(psi))
            obs_array.append(observable_value(Obs_op, psi/np.linalg.norm(psi)))
            
            ent_array.append(entropy_wavefunction(psi/np.linalg.norm(psi),l))
            t += dt
            time_array.append(t)
 
        
    return np.array(time_array),np.array(ent_array),np.array(obs_array),np.array(norm_array), np.array(proba_jump_array),np.array(proba_jump_site_array), np.array(time_jump_array)



import os
 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=str, required=True, help="Model name")
    args = parser.parse_args()



    script_dir = os.path.dirname(os.path.abspath(__file__))

    if args.model == "IsingModel" : 
        path_data = '/../build/res/'
        pbr=np.loadtxt(script_dir+path_data+'traj_proba_jump.dat') 
        pbr_site=np.loadtxt(script_dir+path_data+'traj_proba_jump_site.dat') 

        gamma = 0.3
        h = 0.7*2
        kappa = 1.0

        t_max = 10
        L = 8
        l = 2
        dt = 0.001
        nb_part= 4

        t_e, e_e, o_e,n_e, pe, pje, tj_e  = trajectory_evolution_ED(gamma, h, kappa,L,l,dt,t_max,nb_part,pbr, pbr_site, "01010101") 

        obs_df = pd.read_csv(script_dir+path_data+'traj_densities.dat', header=0, skipinitialspace=True)
        obs_df.columns = obs_df.columns.str.strip().str.rstrip(',') 
        ent_df = pd.read_csv( script_dir+path_data+'traj_entanglement.dat', sep=',', skipinitialspace=True,  header=0 )
        ent_df.columns = ent_df.columns.str.strip() 
        time_df = np.loadtxt( script_dir+path_data+'traj_time.dat')


        fig, axs = plt.subplots( 1,2,  figsize=(8,3), dpi = 200, sharex=True)  
        ax = axs[0]
        ax.plot(time_df,ent_df["S_l=2"], label="C++")
        ax.plot(t_e,e_e, color="orange", linestyle="--" , label="ED python")
        ax.set_xlabel("t")
        ax.set_ylabel("S(t)")
        ax.legend(frameon=False)

        ax = axs[1]
        ax.plot(time_df,obs_df["n_i=0"].values, label="C++")
        ax.plot(t_e,o_e, color="orange", linestyle="--" , label="ED python")
        ax.set_xlabel("t")
        ax.set_ylabel("n_0(t)")
        ax.legend(frameon=False)
        
        fig.suptitle("Benchmark " + args.model)
        fig.tight_layout()
        plt.show()

    else:
        print(f"Unsupported model !")
    

if __name__ == "__main__":
    main()