#include "simulation.h"
#include "logger.h"
#include "output_data.h"
#include "config.h"
#include <cassert>


Simulation::Simulation(
        OutputData& output,
        Config& config)
        : output_(output),
        config_(config)
{}


void Simulation::createLogs(int batch_id){
    Logger::spacer();
    Logger::info("[Simulation] : "); 
    Logger::info("***************************************************************");
    Logger::info("                          Start Simulation                          ");
    if(batch_id > -1 ){
        Logger::info(" batch=",batch_id,"       time_max=", time_max(), "      dt=", dt(), "      L=", L());
    }else{
        Logger::info("    time_max=", time_max(), "      dt=", dt(), "      L=", L());
    }
    Logger::info("***************************************************************"); 
    Logger::spacer();
}



void Simulation::init()  { 
        
        
        // Main Simulation parameters 
        no_click_ = config_.checkParam(NO_CLICK_KEY) ? config_.getParam<bool>(NO_CLICK_KEY) : DEFAULT_NO_CLICK; 
        high_order_ = config_.checkParam(HIGH_ORDER_KEY) ? config_.getParam<bool>(HIGH_ORDER_KEY) : DEFAULT_HIGH_ORDER; 
        L_ = config_.checkParam(SIZE_KEY) ? (int)config_.getParam<double>(SIZE_KEY) : DEFAULT_L; 
        dt_ = config_.checkParam(TIME_STEP_KEY) ? config_.getParam<double>(TIME_STEP_KEY) : DEFAULT_DT; 
        time_max_ = config_.checkParam(TIME_MAX_KEY) ? config_.getParam<double>(TIME_MAX_KEY) : DEFAULT_TIME_MAX; 
         

        // Initialize useful matrices
        C_ = new lapack_complex_double[L_*L_*4];
        C_next = new lapack_complex_double[L_*L_*4];
        C_ent = new double[L_*L_*4];
        delta_pm = new double[L_];

        // tmp var for high-order Runge Kutta integration
        k1 = new lapack_complex_double[L_*L_*4];
        k2 = new lapack_complex_double[L_*L_*4];
        k3 = new lapack_complex_double[L_*L_*4];
        k4 = new lapack_complex_double[L_*L_*4];
        C_acc = new lapack_complex_double[L_*L_*4];
        CL = new lapack_complex_double[L_*L_*4];
        CR = new lapack_complex_double[L_*L_*4];
 
        GG = new lapack_complex_double[L_*L_];
        FF = new lapack_complex_double[L_*L_]; 



        
        // Init state configuration
        std::string init_state = config_.checkParam(INIT_STATE_KEY) ? config_.getParam<std::string>(INIT_STATE_KEY) : DEFAULT_INIT_STATE; 
 
        if(init_state == DEFAULT_INIT_STATE || !config_.checkParam(INIT_STATE_KEY) ){
            setNeelInitState();
        }else if(init_state == RANDOM_INIT_STATE){ 
            if(!config_.checkParam(NB_PART_KEY)){
                std::cerr << RED << "Error: "   << RESET << "For a random init state, the number of particle should be specified in nb_part parameter" << std::endl;
                exit(EXIT_FAILURE);
            }
            int nb_part = (int)config_.getParam<double>(NB_PART_KEY) ;
            setRandomInitState(nb_part);
        }else{
            setInitState(init_state);
        }

        
        record_probas_ = config_.checkParam(RECORD_PROBA_KEY) ? config_.getParam<bool>(RECORD_PROBA_KEY) : RECORD_PROBA;
        
        record_entanglement_statistics_ = config_.checkParam(RECORD_ENTANGLEMENT_STATISTIC_KEY) ? config_.getParam<bool>(RECORD_ENTANGLEMENT_STATISTIC_KEY) : false;
        if(record_entanglement_statistics_){

            for (auto &var : config_.monitoring)
            {
                if(var.name == ENTANGLEMENT_STATISTIC_KEY){
                    current_entanglement[var.param] = {0,0};

                }
            } 


        }

        // Check for init order
        is_init = true;

        
} 




void Simulation::reinit()  { 
    if(probas_list_){
        for (auto& inner : *probas_list_) {
            inner.clear();
        } 
    }

} 


void Simulation::clean(){ 
    delete[] C_next;
    delete[] C_;
    delete[] delta_pm;

    delete[] GG;
    delete[] FF;  

    delete[] k1;
    delete[] k2;
    delete[] k3;    
    delete[] k4; 
    delete[] C_acc; 
    delete[] CL; 
    delete[] CR;
}
 

void Simulation::collectData(){
    output_.collectData(*this);
}  

void Simulation::collectProbas(double r, double r_site){
    if(record_probas_){
         output_.collectProbas( *this, r, r_site);
    }
}

void Simulation::collectEntanglementStatistics(){
    if(record_entanglement_statistics_){
         output_.collectEntanglementStatistics( *this);
    }
}

void Simulation::setInitState(std::string binary_psi){ 
    if(binary_psi.length() != L_ ){
        std::cerr << RED << "Error: "   << RESET << "Init state length not matching the system size" << std::endl;
        exit(EXIT_FAILURE);
    }else if(binary_psi.find_first_not_of("01") != std::string::npos){
        std::cerr << RED << "Error: "   << RESET << "Init state description should be only composed of 0 and 1" << std::endl;
        exit(EXIT_FAILURE);
    }else{

        for (int i = 0; i < 2*L_; i++) {
            for (int j = 0; j < 2*L_; j++) {
                C_[i*2*L_ + j ] = 0;
            }
        }
        for (int i = 0; i < L_; i++) {
            if (binary_psi[i] == '0') {
                C_[i*2*L_ + i ] = 0;
                C_[(i+L_)*2*L_ + (i+L_) ] = 1;
            } else {
                C_[i*2*L_ + i ] = 1;
                C_[(i+L_)*2*L_ + (i+L_) ] = 0;
            }
        }
    }
}


void Simulation::setRandomInitState(int nb_part){
    int count_part = 0 ; 
    std::string binary_psi(L_, '0');

    while (count_part < nb_part)
    {
        double rn = (double)rand() / RAND_MAX;
        if (binary_psi[int(rn*L_)] == '0'){
            binary_psi[int(rn*L_)] = '1';
            count_part += 1 ;
        }
    }

    setInitState(binary_psi);
}




void Simulation::setNeelInitState(bool flip){
    std::string binary_psi(L_, '0');


    for (int i = 0; i < L_; i++) {
        if(i %2 == 0 && flip || i%2 != 0 && !flip ){
            binary_psi[i] = '1';
        } 
    }

    setInitState(binary_psi);
}


void Simulation::dk(lapack_complex_double* C, lapack_complex_double* k, const lapack_complex_double* H,double gamma){
    
    
    int size=L_;
    for (size_t i = 0; i < 2*size; i++){
        for (size_t j = 0; j < 2*size; j++){
            int start = 0 ;
            int stop = 3 ;
            lapack_complex_double r = sim_math::zero_c;
            if(i==0 || i == size){start=1;}
            if(i==size-1 || i == 2*size-1){stop =2;}
 
            for (size_t k = start; k < stop; k++){ 
                r +=  H[i*2*size + (i%size)+ k- 1] * C[2*size*( (i%size)+k- 1)+j]+ H[i*2*size+size +  (i%size)+k- 1] * C[2*size*( (i%size)+k- 1+size) + j] ; 
            } 
            CL[i*2*size + j] =  r; 
     }
    }
    

    for (size_t i = 0; i < 2*size; i++){
        for (size_t j = 0; j < 2*size; j++){
            int start = 0 ;
            int stop = 3 ;
            lapack_complex_double r = sim_math::zero_c;
            if(j==0 || j == size){start=1;}
            if(j==size-1 || j == 2*size-1){stop =2;}
 
            for (size_t k = start; k < stop; k++){ 
                r +=  H[((j%size)+ k- 1)*2*size + j] * C[2*size*(i)+ (j%size)+k- 1]+ H[((j%size)+size+ k- 1)*2*size + j] * C[2*size*(i)+ size+(j%size)+k- 1] ; 
            } 
            CR[i*2*size + j] =  r; 
     }
    }

   
    
    for (size_t i = 0; i < 2*size; i++){
        for (size_t j = 0; j < 2*size; j++){
            k[i*2*size + j] =  -2*sim_math::i_c*( CL[i*2*size + j]-CR[i*2*size + j] );  
     }
    } 

    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, &sim_math::one_c,  C, 2*size, C, 2*size, &sim_math::zero_c, CL , 2*size);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, &sim_math::one_c,  C, 2*size, C+size, 2*size, &sim_math::zero_c, CL+size , 2*size);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, &sim_math::one_c,  C+2*size*size, 2*size, C, 2*size, &sim_math::zero_c, CL+2*size*size , 2*size);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, &sim_math::one_c,  C+2*size*size, 2*size, C+size, 2*size, &sim_math::zero_c, CL+2*size*size+size , 2*size);


    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, &sim_math::one_c,  C+size, 2*size, C+2*size*size, 2*size, &sim_math::zero_c, CR , 2*size);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, &sim_math::one_c,  C+size, 2*size, C +2*size*size+size, 2*size, &sim_math::zero_c, CR +size, 2*size);
    cblas_zgemm(CblasRowMajor,CblasNoTrans,  CblasNoTrans, size, size, size, &sim_math::one_c,  C+2*size*size+size, 2*size, C +2*size*size, 2*size, &sim_math::zero_c, CR+2*size*size , 2*size);
    cblas_zgemm(CblasRowMajor,CblasNoTrans,  CblasNoTrans, size, size, size, &sim_math::one_c,  C+2*size*size+size, 2*size, C +2*size*size+size, 2*size, &sim_math::zero_c, CR+2*size*size+size , 2*size);
    
  
    
    for (size_t i = 0; i < 2*size; i++){
        for (size_t j = 0; j < 2*size; j++){
            k[i*2*size+j] +=  gamma*(CL[i*2*size+j]  - CR[i*2*size+j]);
            if(i < size && j < size) {
                k[i*2*size+j] -=  gamma*C[i*2*size+j] ;
            } 
            if( i >= size && j >= size){
                k[i*2*size+j] +=  gamma*C[i*2*size+j] ;
            } 
        }
    }

    // print_MatC(k);
    // std::exit(0); 
    // print_MatC(k);
}


// Move to math ?? 
void Simulation::linComb(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double alpha, lapack_complex_double *B, int n){
    for (size_t i = 0; i < n; i++){
            for (size_t j = 0; j < n; j++){
                Res[i*n+j] = A[i*n+j] + alpha*B[i*n+j] ; 
            }
        }
}


void Simulation::RK4_step(const sim_math::Matrix& H, double dt_eff, double gamma)  {

    dk(C_, k1,H.mat.get(), gamma ) ;
    linComb(C_acc,C_, dt_eff/2.0,k1,2*L_);
    dk(C_acc,k2,H.mat.get(), gamma );
    linComb(C_acc,C_, dt_eff/2.0,k2,2*L_);
    dk(C_acc,k3,H.mat.get(), gamma );
    linComb(C_acc,C_, dt_eff,k3,2*L_);
    dk(C_acc, k4,H.mat.get(), gamma );
    for (size_t i = 0; i < 2*L_; i++){
        for (size_t j = 0; j < 2*L_; j++){
            C_next[i* 2*L_+j] = C_[i* 2*L_+j] + (1.0/6.0)* dt_eff*(k1[i* 2*L_+j]+2*k2[i* 2*L_+j]+2*k3[i* 2*L_+j]+k4[i* 2*L_+j]) ; 
        }
    }

    //  print_MatC(C_next);
    // std::exit(0); 

}
 

double Simulation::get_delta_norm(){
  	double delta = 0.0 ;
	for  (size_t j = 0; j < L_; j++){
		delta +=  sim_math::creal(C_next[j*2*L_+j]); 
	} 

    return delta;
}

void Simulation::update_C(){
  	sim_math::copy_mat(C_next,C_,2*L_,2*L_); // To update with matrices
}

/**
 * @brief Returns the index corresponding to a random site selection
 *        based on a probability distribution.
 *
 * Converts the input array into a normalized cumulative distribution,
 * then performs a binary search to find the index matching the given
 * random value.
 *
 * @warning Modifies @p distrib in-place (cumulative sum + normalization).
 *
 * @param distrib Array of weights/probabilities (will be modified in-place)
 * @param n       Number of elements in @p distrib
 * @param r_site  Random value in [0, 1] used to sample the distribution
 * @return        Index of the selected site
 */
int Simulation::get_index(double* distrib, int n, double r_site){
    double sum = 0 ;
    double temp = 0;
    for (size_t i = 0; i < n; i++)
    {
        temp = distrib[i];
        distrib[i] = distrib[i] + sum; //cumulative sum computation
        sum += temp;
    }
    sim_math::divide_list(distrib,sum,n); //normalization
    
    return sim_math::array_binarysearch(distrib,n,r_site);
    
}

/**
 * @brief Selects a jump site index by sampling the diagonal of the
 *        correlation matrix C_.
 *
 * Extracts the real part of the diagonal elements of C_ as jump
 * probabilities, then delegates to get_index() to sample the resulting
 * distribution.
 *
 * @note The diagonal elements C_[i*2*L_+i] represent the site occupation
 *       probabilities used to determine the next jump site.
 *
 * @param r_site Random value in [0, 1] used to sample the distribution
 * @return       Index of the selected jump site
 */
int Simulation::get_jump_index(double r_site ){
    for (size_t i = 0; i < L_; i++)
    {
        delta_pm[i] = sim_math::creal(C_[i*2*L_+i]);
    }
    return get_index(delta_pm, L_, r_site);
}


// To refactor
void Simulation::apply_jump(int index){

    double Cll = sim_math::creal(C_[index*L_*2+ index]);

    sim_math::copy_mat_part(C_,GG,2*L_,L_, 0,0,0,0,L_,L_); 

    for (size_t i = 0; i < L_; i++) {
        for (size_t j = 0; j < L_; j++) {
            if (i != index && j!=index){
                GG[i*L_ + j] -= C_[i*L_*2 + index]*C_[index*L_*2 + j]/Cll + C_[i*L_*2 + L_ + index]*sim_math::conj(C_[index*L_*2  + L_ +j])/Cll;
            }else{
                GG[i*L_ + j] = 0  ;
            }
        }
    }
    GG[index*L_+ index] = 1;
    
    sim_math::copy_mat_part(C_,FF,2*L_,L_, 0,L_,0,0,L_,L_); 
    for (size_t i = 0; i < L_; i++) {
        for (size_t j = 0; j < L_; j++) {
            if(i != index && j!=index){
                FF[i*L_ + j] += C_[j*L_*2 + index]*C_[index*L_*2 + L_ + i]/Cll - C_[i*L_*2 + index]*C_[index*L_*2 + L_+ j]/Cll;
            }else{
                FF[i*L_ + j] = 0  ;
            }
            // F[i*L_total + j] -= C[i*L_total*2 + L_total + index]*C[j*L_total*2 + index]/Cll + C[i*L_total*2 + index]*C[index*L_total*2 + L_total +j]/Cll;
        }
    }

     for (size_t i = 0; i < L_; i++) {
        for (size_t j = 0; j < L_; j++) {
            C_[i*L_*2 + j] = GG[i*L_ + j];
            if (i == j ){
                C_[(i+L_)*L_*2 + L_ + j] = 1 - GG[j*L_ + i] ;
            }else{
                C_[(i+L_)*L_*2 + L_ + j] =  - GG[j*L_ + i] ;
            }
            C_[(L_+i)*L_*2 + j] = sim_math::conj(FF[j*L_ + i]);
            C_[(i)*L_*2 + L_ + j] =  FF[i*L_ + j];
        }
    }
    
}


void Simulation::jump(double r_site ){
    index_jump_ = get_jump_index( r_site);    
    apply_jump(index_jump_); 
    
}

double Simulation::get_entanglement( int l){ 
    if(record_entanglement_statistics_ && jump_case_){

        if(before_jump_){
            current_entanglement[l][0] = compute_entanglement(l);
            return  current_entanglement[l][0];
        }else{
            current_entanglement[l][1] = compute_entanglement(l);
            return  current_entanglement[l][1];
        }

    }else{
        return compute_entanglement(l);
    }
}

double Simulation::compute_entanglement( int l ){ 
    int n = 2 * L_;
    for (int i = 0; i < 2*L_; i++) {
        for (int j = 0; j <  2*L_; j++) { 
            C_ent[i*2*L_+j] = 0;
        }
    } 
    
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < l; j++) {
            C_ent[i*2*L_+j] = (2*sim_math::cimag(C_[i*n+j]) + 2*sim_math::cimag(C_[i*n+L_+j]));
            C_ent[(l+i)*2*L_+j+l] = (2*sim_math::cimag(C_[i*n+j]) - 2*sim_math::cimag(C_[i*n+L_+j]));
            if(i == j ){
                C_ent[(i)*2*L_+l+j] = (2*sim_math::creal(C_[i*n+j]) - 2*sim_math::creal(C_[i*n+L_+j]) - 1 );
                C_ent[(i+l)*2*L_+j] = (-2*sim_math::creal(C_[i*n+j]) - 2*sim_math::creal(C_[i*n+L_+j]) + 1 );

            }else{
                C_ent[(i)*2*L_+l+j] = (2*sim_math::creal(C_[i*n+j]) - 2*sim_math::creal(C_[i*n+L_+j]));
                C_ent[(i+l)*2*L_+j] = (-2*sim_math::creal(C_[i*n+j]) - 2*sim_math::creal(C_[i*n+L_+j]));
            }
        }
    }
 

    // double wr[2*l];
    // double wi[2*l]; 
    
    std::vector<double> wr(2*l);
    std::vector<double> wi(2*l);
    std::vector<double> dummy(1);

    assert(L_ >= l);  // catch lda mismatch early

    int info = LAPACKE_dgeev(
        LAPACK_ROW_MAJOR, 'N', 'N',
        2*l,
        C_ent, 2*L_,
        wr.data(), wi.data(),
        dummy.data(), 1,
        dummy.data(), 1
    );

    if (info != 0)
        std::cerr << "dgeev failed w²ith info = " << info << "\n";

    // int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', 2*l, C_ent, 2*L_, wr, wi, nullptr, 2*L_,nullptr, 2*L_); 
    double entropy = 0.0;
    for (int i = 0; i < 2*l; i++) {
        if (wi[i] > 0.0 && wi[i] < 1.0) {
            entropy += -((wi[i]+1.0)/2.0)*std::log(((wi[i]+1.0)/2.0)) - ((-wi[i]+1.0)/2.0)*std::log(((-wi[i]+1.0)/2.0));
        }
    } 

    return entropy; 
}


double Simulation::get_density( int index ){ 
    return sim_math::creal(C_[index*2*L_+index]);
}


std::vector<double>  Simulation::get_densities(){ 
    std::vector<double> densities(L_); 
    for (int i = 0; i < L_; i++) {
        densities[i] = sim_math::creal(C_[i*2*L_+i]);
    }

    return densities;
}