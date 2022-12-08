//#pragma once
//#include "matrix.h"
//#include "vector.h"
//namespace QrDecomposition {
//
//
//    // c = a + b * s
//    static void vmadd(const Vector& a, const Vector& b, double s, Vector& c)
//    {
//        if ((c.size != a.size) || (c.size != b.size)) {
//            std::cerr << "[vmadd]: vector sizes don't match\n";
//            return;
//        }
//
//        for (int i = 0; i < c.size; i++)
//            c(i) = a(i) + s * b(i);
//    }
//
//    // mat = I - 2*v*v^T
//    // !!! m is allocated here !!!
//    static void compute_householder_factor(Matrix& mat, const Vector& v)
//    {
//
//        int n = v.size;
//        mat.allocate(n, n);
//        for (int i = 0; i < n; i++)
//            for (int j = 0; j < n; j++)
//                mat(i, j) = -2 * v(i) * v(j);
//        for (int i = 0; i < n; i++)
//            mat(i, i) += 1;
//    }
//
//
//    static void matrix_show(const Matrix& m, const std::string& str = "")
//    {
//        std::cout << str << "\n";
//        for (int i = 0; i < m.m; i++) {
//            for (int j = 0; j < m.n; j++) {
//                printf(" %8.3f", m(i, j));
//            }
//            printf("\n");
//        }
//        printf("\n");
//    }
//
//    // L2-norm ||A-B||^2
//    static double matrix_compare(const Matrix& A, const Matrix& B) {
//        // matrices must have same size
//        if ((A.m != B.m) || (A.n != B.n))
//            return std::numeric_limits<double>::max();
//
//        double res = 0;
//        for (int i = 0; i < A.m; i++) {
//            for (int j = 0; j < A.n; j++) {
//                res += (A(i, j) - B(i, j)) * (A(i, j) - B(i, j));
//            }
//        }
//
//        res /= A.m * A.n;
//        return res;
//    }
//
//    static void householder(Matrix& mat,
//        Matrix& R,
//        Matrix& Q)
//    {
//
//        int m = mat.m;
//        int n = mat.n;
//
//        // array of factor Q1, Q2, ... Qm
//        std::vector<Matrix> qv(m);
//
//        // temp array
//        Matrix z(mat);
//        Matrix z1;
//
//        for (int k = 0; k < n && k < m - 1; k++) {
//
//            Vector e(m), x(m);
//            double a;
//
//            // compute minor
//            z1.compute_minor(z, k);
//
//            // extract k-th column into x
//            z1.extract_column(x, k);
//
//            a = x.norm();
//            if (mat(k, k) > 0) a = -a;
//
//            for (int i = 0; i < e.size; i++)
//                e(i) = (i == k) ? 1 : 0;
//
//            // e = x + a*e
//            vmadd(x, e, a, e);
//
//            // e = e / ||e||
//            e.rescale_unit();
//
//            // qv[k] = I - 2 *e*e^T
//            compute_householder_factor(qv[k], e);
//
//            // z = qv[k] * z1
//            z.mult(qv[k], z1);
//
//        }
//
//        Q = qv[0];
//
//        // after this loop, we will obtain Q (up to a transpose operation)
//        for (int i = 1; i < n && i < m - 1; i++) {
//
//            z1.mult(qv[i], Q);
//            Q = z1;
//
//        }
//
//      //  R.mult(Q, mat);
//        Q.transpose();
//    }
//
//}