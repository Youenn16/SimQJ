#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "sim_math.h"
#include "logger.h"
#include <lapacke.h> 
#include <map>
#include <array>


class Config;
class OutputData; 

class Simulation {
public:
    constexpr static bool DEFAULT_NO_CLICK=false; 
    constexpr static bool DEFAULT_HIGH_ORDER=true; 
    constexpr static bool RECORD_PROBA=false; 
    constexpr static int DEFAULT_L=8; 
    constexpr static double DEFAULT_DT=0.01; 
    constexpr static double DEFAULT_TIME_MAX=10; 
    constexpr static const char* DEFAULT_INIT_STATE = "neel"; 
    constexpr static const char* RANDOM_INIT_STATE = "random"; 

    explicit Simulation(
                    OutputData& output,
                    Config& config); 

    
    // Jumps probas recording variables
    bool record_probas_=false; 
    bool record_entanglement_statistics_=false;

    std::unique_ptr<std::vector<std::vector<double>>> probas_list_;

    // =========================
    //    Simulation functions 
    // =========================
    void setInitState(std::string binary_psi);
    void setRandomInitState(int nb_part);
    void setRandomInitState() { setRandomInitState(L_ / 2); } 
    void setNeelInitState(bool flip=false);


    void dk( lapack_complex_double* C,lapack_complex_double* k,const lapack_complex_double* H,double gamma);
    void RK4_step(const sim_math::Matrix& H, double dt_eff, double gamma);
    void linComb(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double alpha, lapack_complex_double *B, int n);

    void update_C(); // copy C_next in C 
    double get_delta_norm(); // Norm differential computation 
    int get_index(double* distrib, int n, double r_site); // 
    int get_jump_index(double r_site ); // 
    void apply_jump(int index); //
    void jump(double r_site ); // 


    // Generic methods 
    void init();
    void clean();
    void reinit();
    void createLogs(int batch_id);

    // Saving data methods
    void collectProbas(double r, double r_site);
    void collectData();   
    void collectEntanglementStatistics();   

    // =================
    //    Getters/Setters
    // =================
    double get_entanglement( int l );
    double compute_entanglement( int l );
    double get_density( int index );
    std::vector<double>  get_densities();

    void print_MatC(lapack_complex_double *A){
        std::cout << "[" << std::endl;
        for (int i = 0; i < 2*L_; i++) {
            for (int j = 0; j <  2*L_; j++) {
                std::cout << sim_math::creal(A[i*2*L_+j])<<" + "<< sim_math::cimag(A[i*L_*2+j]) << "j    ";
            }
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << std::endl;
    }

     void print_MatC(double *A){
        std::cout << "[" << std::endl;
        for (int i = 0; i < 2*L_; i++) {
            for (int j = 0; j <  2*L_; j++) {
                std::cout << A[i*2*L_+j]<< " " ;
            }
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << std::endl;
    }

     void print_MatCC(){
        print_MatC(C_);
    }


    const int L() const {
        return L_;
    }

    const double time_max() const {
        return time_max_;
    }

    const double dt() const {
        return dt_;
    }

    const double time() const {
        return time_;
    }

    const double norm() const {
        return norm_;
    }

    void set_time(double time){
        time_ = time;
    }

    void set_norm(double norm){
        norm_ = norm;
    }

    void jump_loc(std::string tag){
        if(tag == "before_jump"){
            jump_case_ = true;
            before_jump_ = true;
        }else if(tag == "after_jump"){
            jump_case_ = true;
            before_jump_ = false;
        }else if(tag == "non_hermi"){
            jump_case_ = false;
            before_jump_ = false;
        }
    }

    const bool jump_case() const {
        return jump_case_;
    }
    
    const bool before_jump() const {
        return before_jump_;
    }

    double get_current_entanglement(int l, int pos){
        return current_entanglement[l][pos];
    }

    int get_index_jump(){
        return index_jump_;
    }
    
    bool has_been_init(){
        return is_init;
    }

     bool is_no_click(){
        return no_click_;
    }

 
protected: 
    OutputData& output_;
    Config& config_;
    bool no_click_=false; 
    bool high_order_=true;

    bool jump_case_ = false;
    bool before_jump_ = false;

    int L_=DEFAULT_L;
    double time_max_=DEFAULT_TIME_MAX;
    double dt_=DEFAULT_DT;
    int index_jump_= 0;

    double time_;
    double norm_;
    
    bool is_init = false;
    

    lapack_complex_double *C_;
    lapack_complex_double *C_next;
    double *C_ent;
    double *delta_pm;

    lapack_complex_double *k1,*k2,*k3,*k4, *C_acc, *CL, *CR;
    lapack_complex_double *GG, *FF;

    std::map<int, std::array<double, 2>> current_entanglement;


   

};