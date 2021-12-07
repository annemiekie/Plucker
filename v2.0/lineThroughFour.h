#ifndef LINEFOUR_H
#define LINEFOUR_H

#include "ray.h"
#include "model.h"

#include <Eigen/Dense>
#include <Eigen/Jacobi>

namespace LineThroughFour {

	static Eigen::MatrixXd makePluckerMat(std::vector<Ray>& lines) {

		Eigen::MatrixXd m(4, 6);

		for (int i = 0; i < lines.size(); i++) {
			lines[i].normalize();
			std::vector<double> plucker = lines[i].plucker();
			for (int j = 0; j < plucker.size(); j++) {
				m(i, j) = plucker[j];
			}
		}
		return m;
	}

	static int singularValueDecomp(Eigen::MatrixXd m, Ray& F, Ray& G, int& dim) {
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
		int rank = svd.nonzeroSingularValues();
		dim = 5 - rank;
		//if (rank < 4) {
		//	std::cout << "rank: " << rank << std::endl;
		//}
		//if (rank > 4) {
		//	std::cout << "rank: " << rank <<  std::endl;
		//}
		Eigen::VectorXd singularvalues = svd.singularValues();

		Eigen::MatrixXd fg = svd.matrixV();
		F = Ray(fg.col(4), true);
		G = Ray(fg.col(5), true);
		return rank;
	}

	static int intersectLinesQuadric3(Ray& F, Ray& G, Ray& H, std::vector<Ray>& intersectLines) {
		return 0;
	}


	static int intersectLinesQuadric(Ray& F, Ray& G, std::vector<Ray>& intersectLines) {
		double a = F.sideVal(F);
		double b = F.sideVal(G);
		double c = G.sideVal(G);
		if (abs(a) < 1E-10 && abs(c) < 1E-10) {
			intersectLines.push_back(F);
			intersectLines.push_back(G);
			return 2;
		}
		double delta = b * b - a * c;
		if (delta < 0) return 0;
		double sq = sqrt(delta);
		double t0 = (-b - sq) / a;
		double t1 = (-b + sq) / a;
		Ray l0, l1;
		l0.u = F.u * t0 + G.u;
		l0.v = F.v * t0 + G.v;
		l1.v = F.v * t1 + G.v;
		l1.u = F.u * t1 + G.u;
		l0.normalize();
		l1.normalize();
		intersectLines.push_back(l0);
		if (delta > 1E-6) intersectLines.push_back(l1);
		return delta == 0.0 ? 1 : 2;
	}

	static void findNullSpace(Eigen::MatrixXd& m, Ray &F, Ray &G, Ray &H, int rank = 4) {
		Eigen::FullPivLU<Eigen::MatrixXd> lu(m);
		Eigen::MatrixXd A_null_space = lu.kernel();

		F = Ray(A_null_space.col(0), true);
		G = Ray(A_null_space.col(1), true);
		if (rank == 3) H = Ray(A_null_space.col(2), true);
	}

	static void findNullSpace2(Eigen::MatrixXd& m, Ray& F, Ray& G, Ray &H) {
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
		cod.compute(m);
		// Find URV^T
		Eigen::MatrixXd V = cod.matrixZ().transpose();
		Eigen::MatrixXd A_null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
		Eigen::MatrixXd P = cod.colsPermutation();
		A_null_space = P * A_null_space; // Unpermute the columns
		F = Ray(A_null_space.col(0), true);
		G = Ray(A_null_space.col(1), true);
	}

	static std::vector<Ray> find(std::vector<Ray>& lines, Model* model) {
		std::vector<Ray> intersectLines;
		//int dim;
		Ray F, G, H;
		Eigen::MatrixXd m = makePluckerMat(lines);

		//int rank = 4;// singularValueDecomp(m, F, G, dim);
		findNullSpace(m, F, G, H);
		//int num;
		//if (rank == 4) 
		intersectLinesQuadric(F, G, intersectLines);
		//else if (rank == 5) num = intersectLinesQuadric3(F, G, H, intersectLines);
		//for (Ray i : intersectLines) {
		//	float sideval = i.sideVal(i);
		//	if (fabs(sideval) > 1E-6) std::cout << "Not on plucker surface " << sideval << std::endl;
		//	for (Ray l : lines) {
		//		sideval = l.sideVal(i);
		//		if (fabs(sideval) > 1E-6) std::cout << "Not really intersecting " << sideval << std::endl;
		//	}
		//}
		return intersectLines;
	}
};

#endif