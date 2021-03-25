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
			std::vector<double> plucker = lines[i].plucker();
			for (int j = 0; j < plucker.size(); j++) {
				m(i, j) = plucker[j];
			}
		}
		return m;
	}

	static void singularValueDecomp(Eigen::MatrixXd m, glm::dvec3& Fu, glm::dvec3& Fv, glm::dvec3& Gu, glm::dvec3& Gv, int& dim) {//& F, Ray& G, int& dim) {
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
		int rank = svd.nonzeroSingularValues();
		dim = 5 - rank;
		if (dim > 1) {
			std::cout << "dim>1" << std::endl;
			return;
		}
		if (dim < 1) {
			std::cout << "dim<1" << std::endl;
		}
		Eigen::VectorXd singularvalues = svd.singularValues();

		Eigen::MatrixXd fg = svd.matrixV();
		Eigen::Vector3d fCol = fg.col(4);
		Eigen::Vector3d gCol = fg.col(5);
		Fu = glm::dvec3(fCol[3], fCol[4], fCol[5]);
		Fv = glm::dvec3(fCol[0], fCol[1], fCol[2]);
		Gu = glm::dvec3(gCol[3], gCol[4], gCol[5]);
		Gv = glm::dvec3(gCol[0], gCol[1], gCol[2]);
	}

	static int intersectLinesQuadric(Ray &F, Ray &G, std::vector<Ray>& intersectLines) {
		float a = F.sideVal(F);
		float b = F.sideVal(G);
		float c = G.sideVal(G);
		float delta = b * b - a * c;
		if (delta < 0) return 0;
		float sq = sqrt(delta);
		float t0 = (-b - sq) / a;
		float t1 = (-b + sq) / a;
		Ray l0, l1;
		l0.u = F.u * t0 + G.u;
		l0.v = F.v * t0 + G.v;
		l1.v = F.v * t1 + G.v;
		l1.u = F.u * t1 + G.u;
		intersectLines.push_back(l0);
		if (delta > 1E-6) intersectLines.push_back(l1);
		return delta == 0.0 ? 1 : 2;
	}

	static int intersectLinesQuadricA(glm::dvec3 Fu, glm::dvec3 Fv, glm::dvec3 Gu, glm::dvec3 Gv, std::vector<Ray>& intersectLines) {
		double a = 2. * glm::dot(Fu, Fv);//F.sideVal(F);
		double b = glm::dot(Fu, Gv) + glm::dot(Fv, Gu);// F.sideVal(G);
		double c = 2. * glm::dot(Gu, Gv); //G.sideVal(G);
		double delta = b * b - a * c;
		if (delta < 0) return 0;
		double sq = sqrt(delta);
		double t0 = (-b - sq) / a;
		double t1 = (-b + sq) / a;
		Ray l0, l1;
		l0.u = Fu * t0 + Gu;
		l0.v = Fv * t0 + Gv;
		l1.v = Fv * t1 + Gv;
		l1.u = Fu * t1 + Gu;
		
		//std::cout << glm::dot(Fu * t0 + Gu, Fv * t0 + Gv) << std::endl;

		intersectLines.push_back(l0);
		if (delta > 1E-6) intersectLines.push_back(l1);
		return delta == 0.0 ? 1 : 2;
	}

	static void findNullSpace(Eigen::MatrixXd m, glm::dvec3 &Fu, glm::dvec3 &Fv, glm::dvec3 &Gu, glm::dvec3 &Gv) {
		Eigen::FullPivLU<Eigen::MatrixXd> lu(m);
		Eigen::MatrixXd A_null_space = lu.kernel();

		//Eigen::Vector3d fCol = A_null_space.col(0);
		//Eigen::Vector3d gCol = A_null_space.col(1);

		Fu = glm::dvec3(A_null_space.col(0)[3], A_null_space.col(0)[4], A_null_space.col(0)[5]);
		Fv = glm::dvec3(A_null_space.col(0)[0], A_null_space.col(0)[1], A_null_space.col(0)[2]);
		Gu = glm::dvec3(A_null_space.col(1)[3], A_null_space.col(1)[4], A_null_space.col(1)[5]);
		Gv = glm::dvec3(A_null_space.col(1)[0], A_null_space.col(1)[1], A_null_space.col(1)[2]);
		//l1 = Ray(A_null_space.col(0), true);
		//l2 = Ray(A_null_space.col(1), true);
	}

	static void findNullSpace2(Eigen::MatrixXd m, glm::dvec3& Fu, glm::dvec3& Fv, glm::dvec3& Gu, glm::dvec3& Gv) {//Ray& l1, Ray& l2) {
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
		cod.compute(m);
		// Find URV^T
		Eigen::MatrixXd V = cod.matrixZ().transpose();
		Eigen::MatrixXd A_null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
		Eigen::MatrixXd P = cod.colsPermutation();
		A_null_space = P * A_null_space; // Unpermute the columns

		Eigen::Vector3d fCol = A_null_space.col(0);
		Eigen::Vector3d gCol = A_null_space.col(1);

		Fu = glm::dvec3(fCol[3], fCol[4], fCol[5]);
		Fv = glm::dvec3(fCol[0], fCol[1], fCol[2]);
		Gu = glm::dvec3(gCol[3], gCol[4], gCol[5]);
		Gv = glm::dvec3(gCol[0], gCol[1], gCol[2]);

		//l1 = Ray(Null_space.col(0), true);
		//l2 = Ray(Null_space.col(1), true);
	}

	static std::vector<Ray> find(std::vector<Ray>& lines, Model* model) {
		std::vector<Ray> intersectLines;
		int dim;
		Ray F, G;
		Eigen::MatrixXd m = makePluckerMat(lines);
		glm::dvec3 Fu, Fv, Gu, Gv;

		//singularValueDecomp(m, Fu, Fv, Gu, Gv, dim);
		findNullSpace(m, Fu, Fv, Gu, Gv);
		int num = intersectLinesQuadricA(Fu, Fv, Gu, Gv, intersectLines);
		
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