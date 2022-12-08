//#pragma once
//#include <cstdio>
//#include <cstdlib>
//#include <cstring> // for memset
//#include <limits>
//#include <iostream>
//#include <vector>
//
//#include <math.h>
//
//#include "vector.h"
//
//class Matrix {
//
//public:
//    // default constructor (don't allocate)
//    Matrix() : m(0), n(0), data(nullptr) {}
//
//    // constructor with memory allocation, initialized to zero
//    Matrix(int m_, int n_) : Matrix() {
//        m = m_;
//        n = n_;
//        allocate(m_, n_);
//    }
//
//    // copy constructor
//    Matrix(const Matrix& mat) : Matrix(mat.m, mat.n) {
//
//        for (int i = 0; i < m; i++)
//            for (int j = 0; j < n; j++)
//                (*this)(i, j) = mat(i, j);
//    }
//
//    // constructor from array
//    template<int rows, int cols>
//    Matrix(double(&a)[rows][cols]) : Matrix(rows, cols) {
//
//        for (int i = 0; i < m; i++)
//            for (int j = 0; j < n; j++)
//                (*this)(i, j) = a[i][j];
//    }
//
//    // destructor
//    ~Matrix() {
//        deallocate();
//    }
//
//
//    // access data operators
//    double& operator() (int i, int j) {
//        return data[i + m * j];
//    }
//    double  operator() (int i, int j) const {
//        return data[i + m * j];
//    }
//
//    // operator assignment
//    Matrix& operator=(const Matrix& source) {
//
//        // self-assignment check
//        if (this != &source) {
//            if ((m * n) != (source.m * source.n)) { // storage cannot be reused
//                allocate(source.m, source.n);          // re-allocate storage
//            }
//            // storage can be used, copy data
//            std::copy(source.data, source.data + source.m * source.n, data);
//        }
//        return *this;
//    }
//
//    // compute minor
//    void compute_minor(const Matrix& mat, int d) {
//
//        allocate(mat.m, mat.n);
//
//        for (int i = 0; i < d; i++)
//            (*this)(i, i) = 1.0;
//        for (int i = d; i < mat.m; i++)
//            for (int j = d; j < mat.n; j++)
//                (*this)(i, j) = mat(i, j);
//
//    }
//
//    // Matrix multiplication
//    // c = a * b
//    // c will be re-allocated here
//    void mult(const Matrix& a, const Matrix& b) {
//
//        if (a.n != b.m) {
//            std::cerr << "Matrix multiplication not possible, sizes don't match !\n";
//            return;
//        }
//
//        // reallocate ourself if necessary i.e. current Matrix has not valid sizes
//        if ((a.m != m )|| (b.n != n))
//            allocate(a.m, b.n);
//
//        memset(data, 0, m * n * sizeof(double));
//
//        for (int i = 0; i < a.m; i++)
//            for (int j = 0; j < b.n; j++)
//                for (int k = 0; k < a.n; k++)
//                    (*this)(i, j) += a(i, k) * b(k, j);
//
//    }
//
//    void transpose() {
//        for (int i = 0; i < m; i++) {
//            for (int j = 0; j < i; j++) {
//                double t = (*this)(i, j);
//                (*this)(i, j) = (*this)(j, i);
//                (*this)(j, i) = t;
//            }
//        }
//    }
//
//    // take c-th column of m, put in v
//    void extract_column(Vector& v, int c);
//
//    // memory allocation
//    void allocate(int m_, int n_) {
//
//        // if already allocated, memory is freed
//        deallocate();
//
//        // new sizes
//        m = m_;
//        n = n_;
//
//        data = new double[m_ * n_];
//        memset(data, 0, m_ * n_ * sizeof(double));
//
//    } // allocate
//
//    // memory free
//    void deallocate() {
//
//        if (data)
//            delete[] data;
//
//        data = nullptr;
//
//    }
//
//    // take c-th column of a matrix, put results in Vector v
//    void extract_column(Vector& v, int c) {
//        if (m != v.size) {
//            std::cerr << "[Matrix::extract_column]: Matrix and Vector sizes don't match\n";
//            return;
//        }
//
//        for (int i = 0; i < m; i++)
//            v(i) = (*this)(i, c);
//    }
//
//
//    int m, n;
//
//private:
//    double* data;
//
//}; // struct Matrix
//
//
////double in[][3] = {
////  { 12, -51,   4},
////  {  6, 167, -68},
////  { -4,  24, -41},
////  { -1,   1,   0},
////  {  2,   0,   3},
////};
//
////int main()
////{
////    Matrix A(in);
////    Matrix Q, R;
////
////    matrix_show(A, "A");
////
////    // compute QR decompostion
////    householder(A, R, Q);
////
////    matrix_show(Q, "Q");
////    matrix_show(R, "R");
////
////    // compare Q*R to the original matrix A
////    Matrix A_check;
////    A_check.mult(Q, R);
////
////    // compute L2 norm ||A-A_check||^2
////    double l2 = matrix_compare(A, A_check);
////
////    // display Q*R
////    matrix_show(A_check, l2 < 1e-12 ? "A == Q * R ? yes" : "A == Q * R ? no");
////
////    return EXIT_SUCCESS;
////}
