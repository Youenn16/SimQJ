#include "sim_math.h"
#include <stdexcept>
#include <iostream>

namespace sim_math {

const lapack_complex_double zero_c  = {0.0,0.0};
const lapack_complex_double i_c = {0.0,1.0};
const lapack_complex_double one_c  = {1.0,0.0};



double add(double a, double b) {
return a + b;
}


double sub(double a, double b) {
return a - b;
}


double mul(double a, double b) {
return a * b;
}


double div(double a, double b) {
if (b == 0.0) throw std::runtime_error("division by zero");
return a / b;
}


double creal(lapack_complex_double z ){
    return ((double*)&z)[0];
}

double cimag(lapack_complex_double z ){
    return ((double*)&z)[1];
}

double cabs(lapack_complex_double z ){
    return sqrt(((double*)&z)[1]*((double*)&z)[1] + ((double*)&z)[0]*((double*)&z)[0]);
} 


// Conjugation of complex number
lapack_complex_double conj(lapack_complex_double a){
    return creal(a)*one_c - i_c* cimag(a);
}


// Compute the frobenius norm of a matrix
double norm_1(lapack_complex_double *A,int n,int m){
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            norm +=  creal(A[i*m + j])*creal(A[i*m + j])+cimag(A[i*m + j])*cimag(A[i*m + j]);
        }
    }
    return sqrt(norm);
}

// Compute the frobenius norm of a square matrix
double norm_1(lapack_complex_double *A,int n){
    return norm_1(A,n,n);
}

// MATRIX OPERATIONS 

// Square zero matrix filling (double)
void zero_matrix(lapack_complex_double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i*n+j] = zero_c;
        }
    }
}
void zero_matrix(Matrix& A) { 
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            A(i,j) = zero_c;
        }
    }
}

// Square identity matrix filling (complex)
void identity_matrix(lapack_complex_double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j){
                A[i*n+j] = one_c;
            }else{
                A[i*n+j] = zero_c;
            }
        }
    }
}

// Square random matrix filling (complex)
void random_matrix(lapack_complex_double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            lapack_complex_double r =  ((double)rand() / RAND_MAX)*one_c+i_c*((double)rand() / RAND_MAX) ;
            A[i*n+j] = r;
        }
    }
}


// Matrix copy (complex) 
void copy_mat(lapack_complex_double *Res, lapack_complex_double *A, int n,int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A[i*m + j] = Res[i*m + j];
        }
    }
}

// Square matrix copy (complex) 
void copy_mat(lapack_complex_double *Res, lapack_complex_double *A, int n){
    copy_mat(Res,A,n,n);
}

// Matrix copy part (complex) 
void copy_mat_part(lapack_complex_double *Res, lapack_complex_double *A, int n_ref, int n_targ , int r_ref, int c_ref,int r_targ, int c_targ, int r_size, int c_size ){
    for (int i = 0; i < r_size; i++) {
        for (int j = 0; j < c_size; j++) { 
            A[(i+r_targ)*n_targ + (j+c_targ)] = Res[(i+r_ref)*n_ref + (j+c_ref)];
        }
    }
}

//Matrix Conjugation (complex)
void conj(lapack_complex_double *A, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            A[i*m+j] = creal(A[i*m+j])*one_c - i_c* cimag(A[i*m+j]);
        }
    }
}

// Square Matrix Conjugation (complex)
void conj(lapack_complex_double *A, int n)
{
    conj(A,n,n);
}

// Square Matrix Transpose (complex)
void transpose(lapack_complex_double* A, int n)
{
    lapack_complex_double temp;
    for(int i = 0; i < n; i++)
    {
        for(int j = i+1; j < n; j++)
        {
            temp = A[i*n+j];
            A[i*n+j] = A[j*n+i];
            A[j*n+i] = temp;
        }
    }
}


// Square Matrix multiplication (complex) A*B = Res 
void multmat(lapack_complex_double* Res, lapack_complex_double* A,lapack_complex_double* B, int n){
     cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one_c, A, n, B, n, &zero_c, Res, n);
}
 
// Matrix multiplication (complex) with A: nxm , B: m*n and Res : nxn =>  A*B = Res 
void multmat(lapack_complex_double* Res, lapack_complex_double* A,lapack_complex_double* B, int n, int m){
     cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, m, &one_c, A, n, B, m, &zero_c, Res, n);
}
 
// Get the eigenvalues (w) and eigenvectors (vr) of a matrix (complex)
void eig(lapack_complex_double* A, int n, lapack_complex_double* w, lapack_complex_double* vr){
    char jobvl = 'N';  // Compute eigenvalues only (left)
    char jobvr = 'V';  // Compute eigenvalues and right eigenvectors 

    lapack_complex_double vl[n]; // left eigenvectors
    int ldvl = 1; //leading dimension of the array VL

    int ldvr = n; //leading dimension of the array VR
    int info;
    //zgeev parameters ::  row size , matrix , column size 
    info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, jobvl, jobvr, n, A, n, w, nullptr, ldvr, vr, ldvr); 
}
 
// Ordering the eigenvalues (w) and eigenvectors (vr) of a matrix (complex)
void eig_ordered(lapack_complex_double* A, int n, lapack_complex_double* w, lapack_complex_double* vr){
    eig(A, n, w, vr);
    // sort eigenvalues and eigenvectors     

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if ( creal(w[i]) > creal(w[j]) || cabs(w[i]) > cabs(w[j])) { 
                lapack_complex_double temp = w[i];
                w[i] = w[j];
                w[j] = temp;
                for (int k = 0; k < n; k++) {
                    temp = vr[k*n+i];
                    vr[k*n+i] = vr[k*n+j];
                    vr[k*n+j] = temp;
                }
            }
        }
    }
}


// Get the eigenvalues (w) and eigenvectors (vr) of a matrix (complex) with eigenvalues bigger that the splitter
// Done for a hermitian A matrix, the result is stored in res, n_h is size of the rows of the new matrix res 
// The use case is typically a projector with eigenvalues 0 and 1, the splitter is 0.5
void eig_splitter(lapack_complex_double* A, int n, int n_h, double* w, lapack_complex_double* res, int count_jump){
    char jobz = 'V';  // Compute eigenvectors 
    char uplo = 'U';  // Store Upper part of A 
    int info;
    info = LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, n, A, n, w); 

    int count = 0;
    double splitter = 0.5;  

    for (int i = 0; i < n; i++) {
        if ( (creal(w[i]) > splitter || cabs(w[i]) > splitter) && count < n_h) { 
            for (int j = 0; j < n; j++) {
                res[j*n_h + count ] =  A[j*n+i];
            }
            count += 1;
        }
    }
}

// Product of a double with a complex matrix 
void prod(lapack_complex_double *A, double factor, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A[i*m + j] = A[i*m + j]*factor;
        }
    }
}

void prod(lapack_complex_double *A, double factor, int n){
    prod(A,factor,n,n);
}

// Product of a complex with a complex matrix 
void prod(lapack_complex_double *A, lapack_complex_double factor, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A[i*m + j] = A[i*m + j]*factor;
        }
    }
}

void prod(lapack_complex_double *A, lapack_complex_double factor, int n){
    prod(A,factor,n,n);
}

// Division of a complex matrix by a double
void divide(lapack_complex_double *A, double factor, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A[i*m + j] = A[i*m + j]/factor;
        }
    }
}

void divide(lapack_complex_double *A, double factor, int n){
    divide(A,factor,n,n);
}

// Addition of two complex matrices 
void add(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            Res[i*m + j] = A[i*m + j]+ B[i*m + j];
        }
    }
}

void add(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n){
    add(Res,A,B,n,n);
}

// Substraction of two complex matrices 
void substract(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            Res[i*m + j] = A[i*m + j]- B[i*m + j];
        }
    }
}

void substract(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n){
    substract(Res,A,B,n,n);
}


// Addition of a double to the diagonal of a complex sqaure matrix
void add_diag(lapack_complex_double *A, double factor, int n){
    for (int i = 0; i < n; i++) { 
            A[i*n + i] = A[i*n + i]+ factor;
    }
}

// Compute the inverse matrix of complex square matrix 
void inverse_mat(lapack_complex_double *A,int n){
    int ipiv[n];

    // Call the zgetrf function to compute the LU decomposition of A
    int info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
    if (info > 0) {
        std::cout << "The matrix is singular." << std::endl;
        // return 1;
    }

    // Call the zgetri function to compute the inverse of A
    info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, A, n, ipiv );
    if (info > 0) {
        std::cout << "The matrix is singular." << std::endl;
        // return 1;
    }
}

// PADE METHOD 
const double coeff_pade_full[5][14] = {
    {120., 60., 12., 1.,0,0,0,0,0,0,0,0,0,0},
    {30240., 15120., 3360., 420., 30., 1.,0,0,0,0,0,0,0,0},
    {17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.,0,0,0,0,0,0},
    {17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                2162160., 110880., 3960., 90., 1.,0,0,0,0}, 
    {64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.}};

const double thetap[5] = {0.01495585217958292, 0.25393983300632, 0.950417899616293,2.0978479612570,5.3719203511481};
//function for exponential of a complex matrix (using Pade approximation and squaring)
void expm(lapack_complex_double *A, int n, lapack_complex_double *expA){
    lapack_complex_double *A2,*U, *V,*UVp, *UVm, *M_temp, *Res;
    A2 = new lapack_complex_double[n*n];
    U = new lapack_complex_double[n*n];
    V = new lapack_complex_double[n*n];
    UVp = new lapack_complex_double[n*n];
    UVm = new lapack_complex_double[n*n];
    M_temp = new lapack_complex_double[n*n];
    Res = new lapack_complex_double[n*n];
    lapack_complex_double *temp;

    int s = 0; 
    
    copy_mat(A,expA,n); 
    double A_L1 = norm_1(expA,n);

    int m = 3 ; 
    int nb_pade = 0;
    if (A_L1 > thetap[0] && A_L1 <  thetap[1] ){
        m = 5 ; 
        nb_pade = 1;
    }else if(A_L1 > thetap[1] && A_L1 < thetap[2]){
        m = 7 ; 
        nb_pade = 2;
    }else if(A_L1 > thetap[2] && A_L1 < thetap[3]){
        m = 7 ; 
        nb_pade = 3;
    }else if(A_L1 > thetap[2] && A_L1 < thetap[3]){
        m = 9 ; 
        nb_pade = 4;
    }else{
        s = 1; 
        divide(expA,2.0,n); 
        while (norm_1(expA,n) > thetap[4]){
            s += 1;
            divide(expA,2.0,n);
        } 
        m = 13 ; 
        nb_pade = 4;
    }
    
    
    // identity_matrix(U,n);


    multmat(A2,expA,expA,n); // A2 = A*A 
    copy_mat(A2,U,n);
    prod(U,coeff_pade_full[nb_pade][m],n);  // U = A * b[9]

    copy_mat(A2,V,n);
    prod(V,coeff_pade_full[nb_pade][m-1],n); // V = 1 * b[8]
    

    int m_ = int(m / 2.0); 

    for (int i = 0; i < m_; i++) {
        add_diag(U,coeff_pade_full[nb_pade][2*(m_-i)-1],n);

        if( i < m_ -1){
            multmat(M_temp,U,A2,n);
            temp = U;
            U = M_temp;
            M_temp = temp;
        }
       
    }
    multmat(M_temp,U,expA,n);
    temp = U;
    U = M_temp;
    M_temp = temp; 

    for (int i = 0; i < m_; i++) {
        add_diag(V,coeff_pade_full[nb_pade][2*(m_-i-1)],n);

        if( i < m_ -1){
            multmat(M_temp,V,A2,n);
            temp = V;
            V = M_temp;
            M_temp = temp;
        }
    }
    
    add(UVp,U,V,n); // P = U + V 
    substract(UVm,U,V,n); // Q = -U +V
    inverse_mat(UVm,n);
    multmat(Res,UVm,UVp,n); // res = Q^{-1} P 
    for (size_t i = 0; i < s; i++)
    {
        multmat(M_temp,Res,Res,n);
        temp = Res;
        Res = M_temp;
        M_temp = temp;
    }
    copy_mat(Res,expA,n);
    
}


// LIST METHODS 

// Division of a double list by a double
void divide_list(double *A, double factor, int n){
    for (int j = 0; j < n; j++) {
            A[j] = A[j]/factor;
    }
}

int binarySearch(double* sortedArray, int left, int right, double target) {
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (sortedArray[mid] < target) {
            left = mid + 1;
        } else if (sortedArray[mid] > target && (mid == 0 || sortedArray[mid-1] < target)) {
            // Found the index where the double value is greater than the value at this index and smaller than the one at index+1
            return mid;
        } else {
            right = mid - 1;
        }
    } 
    return -1;
}
int array_binarysearch(double* sortedArray, int size, double target) {
    return binarySearch(sortedArray, 0, size-1, target);
}

}