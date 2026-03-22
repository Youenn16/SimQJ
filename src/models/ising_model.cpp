#pragma once
#include "models/ising_model.h"
#include "simulation.h"
#include "logger.h"
#include "utils.h"
#include <chrono> 



// IsingModel  IsingModel::create(Config& config) 
//     {
//         if (config.checkParam("h") & config.checkParam("gamma") & config.checkParam("kappa")) {
//             // Default values 
//             double J = config.checkParam("J") ? config.getParam("J") : IsingModel::DEFAULT_J;
//             return IsingModel(config.getParam("gamma"),config.getParam("h"),config.getParam("kappa"), J ); 
        
//         }else{
//             throw std::invalid_argument("Some parameters are missing for Ising Model");
//         }
        
//     }

 

IsingModel::IsingModel(Config& config): Model<IsingModel>("ising_chain","density", false), gamma_(0), J_(0), h_(0), kappa_(0){
        // Default values 
    J_ = config.checkParam("J") ? config.getParam<double>("J") : IsingModel::DEFAULT_J; 

    if(!config.checkParam("gamma") || !config.checkParam("h") || !config.checkParam("kappa")){
        Logger::error("Missing at least one of model parameters (gamma, h, kappa)" );
        exit(EXIT_FAILURE);
    }


    gamma_ = config.getParam<double>("gamma");
    h_ = config.getParam<double>("h");
    kappa_ = config.getParam<double>("kappa"); 
    createLogs();
    
    
} 
 

IsingModel::IsingModel(double gamma,  double h, double kappa, double J)
    : Model<IsingModel>("ising_chain","density", false), gamma_(gamma), J_(J), h_(h), kappa_(kappa) {
        createLogs();
}

void IsingModel::createLogs(){
        Logger::info("[Init]");
        Logger::info_big("Ising Model", MAGENTA, "h=",h_,"gamma=",gamma_,"kappa=", kappa_, "J=", J_); 
}


void IsingModel::init(Simulation& simu){
    if(!simu.has_been_init() ){
        Logger::error("Simulation should be initialized before model initialization !" );
        exit(EXIT_FAILURE);
    } 


    int L = simu.L();
     
    H = std::unique_ptr<sim_math::Matrix>(new sim_math::Matrix(2*L,2*L));
    

    zero_matrix(*H); 


    for (size_t i = 0; i < L; i++)
    {
        (*H)(i,i) += -h_;  
        (*H)(i+L,i+L) += h_;     
    }

    for (size_t i = 0; i < L-1; i++)
    {
        (*H)(i,i+1) += 0.5;
        (*H)(i+1,i) += 0.5;     

        (*H)(i+1,i+L) += -0.5*kappa_;     
        (*H)(i,i+L+1) += 0.5*kappa_;  
        
        (*H)(i+L+1,i) += 0.5*kappa_;     
        (*H)(i+L,i+1) += -0.5*kappa_;    

        (*H)(i+L,i+L+1) += -0.5;
        (*H)(i+1+L,L+i) += -0.5; 
    }
 

    // if(pbc){
    //     H[(L_total-1)*2*L_total + 0] += 0.5;
    //     H[(0)*2*L_total + L_total-1] += 0.5;     

    //     H[(0)*2*L_total + L_total + L_total-1] += -0.5*kappa;     
    //     H[(L_total-1)*2*L_total + L_total + 0] += 0.5*kappa;  
        
    //     H[(0)*2*L_total + 2*L_total*L_total + L_total-1] += 0.5*kappa;     
    //     H[(L_total-1)*2*L_total + 2*L_total*L_total + 0] += -0.5*kappa;    

    //     H[(L_total-1)*2*L_total + 2*L_total*L_total +L_total+ 0] += -0.5;
    //     H[(0)*2*L_total + 2*L_total*L_total +L_total+ L_total-1] += -0.5; 
    // } 

    Logger::info("[Init] Simulation and Model initialized ! ");

}

void IsingModel::simulate(Simulation& simu){
        int L = simu.L();
        double dt_eff = simu.dt();
        double epsilon = 0;
        double dt_curr = simu.dt();
        simu.set_time(0);
       

        // Temp to move 
        int N_prec = 15 ;
        double error_sampling = 0.00001;
        double N = 1;
        double N_next = 1;
        simu.set_norm(1);

        int index_prec = 1;
        int count_jump = 0;


        double r = (double)rand() / RAND_MAX;
        double r_site = (double)rand() / RAND_MAX;
        simu.collectProbas(r,r_site);
        simu.collectData();

        int count_message = 0;
 
        
        auto start = std::chrono::high_resolution_clock::now();

        
        Logger::info("[Simulation] Init state : |", utils::to_string(simu.get_densities(), "") , ">" ); 
        Logger::info("[Simulation] Progress"); 

        while (simu.time() < simu.time_max() + simu.dt() ){

            Logger::loadingBar( simu.time()/simu.time_max()*100,100, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.0);
 
 
            simu.jump_loc("before_jump");
            simu.RK4_step(*H, simu.dt(), gamma_); 

            N_next = N - N *gamma_ *simu.dt()*simu.get_delta_norm() ;

 
            if (r > N_next && !simu.is_no_click()){ 
                dt_eff = 0;
                epsilon = abs(r-N_next);
                index_prec = 1; 
                while (epsilon > error_sampling && index_prec <  N_prec ){
                    if(r < N_next ){
                        dt_eff += simu.dt()/(pow(2,index_prec-1)) ;
                        simu.update_C(); // Copy C_next to C 
                        
                        N = N_next;
                    }else{
                        index_prec += 1;   
                    }
                    dt_curr = simu.dt()/(pow(2,index_prec));
                    simu.RK4_step(*H, dt_curr, gamma_);
                    N_next = N - N *gamma_ *dt_curr*simu.get_delta_norm() ;
                    epsilon = abs(r-N_next);
                }
                
                simu.update_C();
                simu.set_norm(N_next);
                simu.set_time(simu.time()+ dt_eff+ simu.dt()/(pow(2.0,index_prec-1))); 

                simu.collectData();
                simu.collectEntanglementStatistics();
 
                simu.jump(r_site);
                simu.jump_loc("after_jump");

                N = 1; 
                simu.set_norm(1);
                simu.set_time(simu.time()+ 0.00001);
                simu.collectData();
                simu.collectEntanglementStatistics();

                count_jump+=1;
                r = (double)rand() / RAND_MAX;
                r_site = (double)rand() / RAND_MAX;  

                simu.collectProbas(r,r_site);
                
                simu.jump_loc("non_hermi");
            }else{
                N = N_next;
                simu.update_C();
                simu.set_norm(N_next);
                simu.set_time(simu.time()+ simu.dt());

                simu.collectData();

            }
    }

}