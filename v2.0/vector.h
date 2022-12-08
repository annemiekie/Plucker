//#include <cstdio>
//#include <cstdlib>
//#include <cstring> // for memset
//#include <limits>
//#include <iostream>
//#include <vector>
//
//// column vector
//class Vector {
//
//public:
//    // default constructor (don't allocate)
//    Vector() : size(0), data(nullptr) {}
//
//    // constructor with memory allocation, initialized to zero
//    Vector(int size_) : Vector() {
//        size = size_;
//        allocate(size_);
//    }
//
//    // destructor
//    ~Vector() {
//        deallocate();
//    }
//
//    // access data operators
//    double& operator() (int i) {
//        return data[i];
//    }
//    double  operator() (int i) const {
//        return data[i];
//    }
//
//    // operator assignment
//    Vector& operator=(const Vector& source) {
//
//        // self-assignment check
//        if (this != &source) {
//            if (size != (source.size)) {   // storage cannot be reused
//                allocate(source.size);         // re-allocate storage
//            }
//            // storage can be used, copy data
//            std::copy(source.data, source.data + source.size, data);
//        }
//        return *this;
//    }
//
//    // memory allocation
//    void allocate(int size_) {
//
//        deallocate();
//
//        // new sizes
//        size = size_;
//
//        data = new double[size_];
//        memset(data, 0, size_ * sizeof(double));
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
//    //   ||x||
//    double norm() {
//        double sum = 0;
//        for (int i = 0; i < size; i++) sum += (*this)(i) * (*this)(i);
//        return sqrt(sum);
//    }
//
//    // divide data by factor
//    void rescale(double factor) {
//        for (int i = 0; i < size; i++) (*this)(i) /= factor;
//    }
//
//    void rescale_unit() {
//        double factor = norm();
//        rescale(factor);
//    }
//
//    int size;
//
//private:
//    double* data;
//
//}; // class Vector