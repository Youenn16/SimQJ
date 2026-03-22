#pragma once 
#include <iostream> 
#include <cstdlib>
#include <cmath>
#include <complex>

// Now safe to include OpenBLAS headers
#include <lapacke.h>
#include <cblas.h>
#include <memory>


namespace sim_math {

extern const lapack_complex_double i_c ;
extern const lapack_complex_double one_c ;
extern const lapack_complex_double zero_c ;



double add(double a, double b);
double sub(double a, double b);
double mul(double a, double b);
double div(double a, double b);

double creal(lapack_complex_double z );
double cimag(lapack_complex_double z );
double cabs(lapack_complex_double z );


// MATRIX STRUCTURE 

struct Matrix {
    int rows, cols;
    // lapack_complex_double* mat;
     std::unique_ptr<lapack_complex_double[]> mat;

    Matrix() = default;

    Matrix(size_t r, size_t c)
        : rows(r), cols(c),
          mat(new lapack_complex_double[r * c]) {}

    Matrix(Matrix&&) noexcept = default;
    Matrix& operator=(Matrix&&) noexcept = default;

    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;

    lapack_complex_double& operator()(size_t i, size_t j) {
        return mat[i*cols + j];
    }

    const lapack_complex_double& operator()(size_t i, size_t j) const {
        return mat[i*cols + j ];
    }

    lapack_complex_double* data() { return mat.get(); }
    const lapack_complex_double* data() const { return mat.get(); }


    void print() const { 
        std::cout << "[" << std::endl;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << sim_math::creal(mat[i*cols+j])<<" + "<< sim_math::cimag(mat[i*cols+j]) << "j    ";
            }
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << std::endl;
    }
    
};



double norm_1(lapack_complex_double *A,int n,int m);
double norm_1(lapack_complex_double *A,int n);

void zero_matrix(lapack_complex_double *A, int n);
void zero_matrix(Matrix& A);

void identity_matrix(lapack_complex_double *A, int n);
void random_matrix(lapack_complex_double *A, int n) ;
void copy_mat(lapack_complex_double *Res, lapack_complex_double *A, int n,int m);
void copy_mat(lapack_complex_double *Res, lapack_complex_double *A, int n);
void copy_mat_part(lapack_complex_double *Res, lapack_complex_double *A, int n_ref, int n_targ , int r_ref, int c_ref,int r_targ, int c_targ, int r_size, int c_size );
void conj(lapack_complex_double *A, int n, int m);
void conj(lapack_complex_double *A, int n);

lapack_complex_double conj(lapack_complex_double a);
void transpose(lapack_complex_double* A, int n);


void multmat(lapack_complex_double* Res, lapack_complex_double* A,lapack_complex_double* B, int n);
void multmat(lapack_complex_double* Res, lapack_complex_double* A,lapack_complex_double* B, int n, int m);
void eig(lapack_complex_double* A, int n, lapack_complex_double* w, lapack_complex_double* vr);
void eig_ordered(lapack_complex_double* A, int n, lapack_complex_double* w, lapack_complex_double* vr);
void eig_splitter(lapack_complex_double* A, int n, int n_h, double* w, lapack_complex_double* res, int count_jump);


void prod(lapack_complex_double *A, double factor, int n, int m);
void prod(lapack_complex_double *A, double factor, int n);
void prod(lapack_complex_double *A, lapack_complex_double factor, int n, int m);
void prod(lapack_complex_double *A, lapack_complex_double factor, int n);
void divide(lapack_complex_double *A, double factor, int n, int m);
void divide(lapack_complex_double *A, double factor, int n);
void add(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n, int m);
void add(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n);
void substract(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n, int m);
void substract(lapack_complex_double *Res,lapack_complex_double *A, lapack_complex_double *B, int n);

void add_diag(lapack_complex_double *A, double factor, int n);
void inverse_mat(lapack_complex_double *A,int n);

// PADE METHOD  


extern const double coeff_pade_full[5][14] ;
extern const double thetap[5];

void expm(lapack_complex_double *A, int n, lapack_complex_double *expA);




// LIST METHODS 

void divide_list(double *A, double factor, int n);
int binarySearch(double* sortedArray, int left, int right, double target);
int array_binarysearch(double* sortedArray, int size, double target);


}