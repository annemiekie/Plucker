#ifndef LINEFOUR_H
#define LINEFOUR_H

#include "lineThroughFourLines.h"
#include "model.h"
#include "qrDecomposition.h"

#include <Eigen/Dense>
#include <Eigen/Jacobi>

namespace Lines4Finder {

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

	static int singularValueDecomp(Eigen::MatrixXd m, Line4& F, Line4& G) {
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
		int rank = svd.nonzeroSingularValues();
		//dim = 5 - rank;
		Eigen::VectorXd singularvalues = svd.singularValues();

		Eigen::MatrixXd fg = svd.matrixV();
		F = Line4(fg.col(4), true);
		G = Line4(fg.col(5), true);
		return rank;
	}

	static int intersectLinesQuadric(Line4& F, Line4& G, std::vector<Line4>& intersectLines) {
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
		Line4 l0, l1;
		l0.u = F.u * t0 + G.u;
		l0.v = F.v * t0 + G.v;
		l1.v = F.v * t1 + G.v;
		l1.u = F.u * t1 + G.u;
		l0.normalize();
		l1.normalize();
		intersectLines.push_back(l0);
		if (delta > 1E-15) intersectLines.push_back(l1);
		return delta == 0.0 ? 1 : 2;
	}

	static void findNullSpaceLU(Eigen::MatrixXd& m, Line4 &F, Line4 &G, Line4 &H, int rank = 4) {
		Eigen::FullPivLU<Eigen::MatrixXd> lu(m);
		Eigen::MatrixXd A_null_space = lu.kernel();

		F = Line4(A_null_space.col(0), true);
		G = Line4(A_null_space.col(1), true);
		if (rank == 3) H = Line4(A_null_space.col(2), true);
	}

	static void findNullSpaceCOD(Eigen::MatrixXd& m, Line4& F, Line4& G, Line4 &H) {
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
		cod.compute(m);
		// Find URV^T
		Eigen::MatrixXd V = cod.matrixZ().transpose();
		Eigen::MatrixXd A_null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
		Eigen::MatrixXd P = cod.colsPermutation();
		A_null_space = P * A_null_space; // Unpermute the columns
		F = Line4(A_null_space.col(0), true);
		G = Line4(A_null_space.col(1), true);
	}

	static void findNullSpaceQR(Eigen::MatrixXd& m, Line4& F, Line4& G, Line4& H) {
		Eigen::MatrixXd mt = m.transpose();
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(mt);
		Eigen::MatrixXd q = qr.matrixQ();

		F = Line4(q.col(4), true);
		G = Line4(q.col(5), true);
	}

	static std::vector<Line4> find(std::vector<Ray>& lines) {
		std::vector<Line4> intersectLines;
		Line4 F, G, H;
		Eigen::MatrixXd m = makePluckerMat(lines);
		//int rank = singularValueDecomp(m, F, G);
		//if (rank != 4) std::cout << "dim " << rank << std::endl;
		findNullSpaceLU(m, F, G, H);
		//if (rank == 4) 
		intersectLinesQuadric(F, G, intersectLines);
		//else if (rank == 5) num = intersectLinesQuadric3(F, G, H, intersectLines);
		return intersectLines;
	}
};

#endif