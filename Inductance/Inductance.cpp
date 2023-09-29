#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "matrix.h"

#define gnuplot "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist setting.plt"

#define IN
#define OUT
#define Pi 3.14159265
#define mu0 (4 * Pi * 1e-7)

#define r_shield 0.015
#define w_tape 0.012
#define t_tape 0.1e-3
#define ring_st 1
#define ring_sp 19
#define step 500

using namespace std;

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix PosVec_Rad(IN Matrix init, IN Matrix end, IN long n);
Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix TangLinVec_Rad(IN Matrix init, IN Matrix end, IN long n);

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	ofstream ofs_out("output.csv");

	double inductance = 0;
	double dist = 0;
	double dotPro = 0;

	Matrix pos0_center = PosVec_For(0, Pi, step, r_shield, w_tape);
	Matrix pos1_center = PosVec_Rev(-Pi, 0, step, r_shield, w_tape);
	Matrix tang0_center = TangLinVec_For(0, Pi, step, r_shield, w_tape);
	Matrix tang1_center = TangLinVec_Rev(-Pi, 0, step, r_shield, w_tape);

	Matrix pos0_edge = PosVec_For(0, Pi, step, r_shield + t_tape / 2, w_tape + w_tape / 2);
	Matrix pos1_edge = PosVec_Rev(-Pi, 0, step, r_shield + t_tape / 2, w_tape + w_tape / 2);
	Matrix tang0_edge = TangLinVec_For(0, Pi, step, r_shield + t_tape / 2, w_tape + w_tape / 2);
	Matrix tang1_edge = TangLinVec_Rev(-Pi, 0, step, r_shield + t_tape / 2, w_tape + w_tape / 2);

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos0_center[i][0] - pos0_edge[j][0], 2.) +
						pow(pos0_center[i][1] - pos0_edge[j][1], 2.) +
						pow(pos0_center[i][2] - pos0_edge[j][2], 2.));
			dotPro = tang0_center[i][0] * tang0_edge[j][0] +
					 tang0_center[i][1] * tang0_edge[j][1] +
					 tang0_center[i][2] * tang0_edge[j][2];
			inductance += dotPro / dist;
		}
	}

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos0_center[i][0] - pos1_edge[j][0], 2.) +
						pow(pos0_center[i][1] - pos1_edge[j][1], 2.) +
						pow(pos0_center[i][2] - pos1_edge[j][2], 2.));
			dotPro = tang0_center[i][0] * tang1_edge[j][0] +
					 tang0_center[i][1] * tang1_edge[j][1] +
					 tang0_center[i][2] * tang1_edge[j][2];
			inductance += dotPro / dist;
		}
	}

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_center[i][0] - pos0_edge[j][0], 2.) +
						pow(pos1_center[i][1] - pos0_edge[j][1], 2.) +
						pow(pos1_center[i][2] - pos0_edge[j][2], 2.));
			dotPro = tang1_center[i][0] * tang0_edge[j][0] +
					 tang1_center[i][1] * tang0_edge[j][1] +
					 tang1_center[i][2] * tang0_edge[j][2];
			inductance += dotPro / dist;
		}
	}

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_center[i][0] - pos1_edge[j][0], 2.) +
						pow(pos1_center[i][1] - pos1_edge[j][1], 2.) +
						pow(pos1_center[i][2] - pos1_edge[j][2], 2.));
			dotPro = tang1_center[i][0] * tang1_edge[j][0] +
					 tang1_center[i][1] * tang1_edge[j][1] +
					 tang1_center[i][2] * tang1_edge[j][2];
			inductance += dotPro / dist;
		}
	}

	cout << 0 << "/" << ring_sp << " : " << mu0 / (4 * Pi) * inductance << " [H]" << endl;
	ofs_out << 0 << "," << mu0 / (4 * Pi) * inductance << "\n";

	for (int ring = ring_st; ring <= ring_sp; ring++) {
		Matrix pos0_base = PosVec_For(0, Pi, step, r_shield, w_tape);
		Matrix pos1_base = PosVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape);
		Matrix pos_radial01_base = PosVec_Rad(pos0_base.extract_row(step).transposed(), pos1_base.extract_row(1).transposed(), step);
		Matrix pos_radial10_base = PosVec_Rad(pos1_base.extract_row(step).transposed(), pos1_base.extract_row(1).transposed(), step);
		Matrix tang0_base = TangLinVec_For(0, Pi, step, r_shield, w_tape);
		Matrix tang1_base = TangLinVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape);
		Matrix tang_radial01_base = TangLinVec_Rad(pos0_base.extract_row(step).transposed(), pos1_base.extract_row(1).transposed(), step);
		Matrix tang_radial10_base = TangLinVec_Rad(pos1_base.extract_row(step).transposed(), pos0_base.extract_row(1).transposed(), step);
		
		Matrix pos0 = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield, w_tape);
		Matrix pos1 = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + t_tape, w_tape);
		Matrix pos_radial01 = PosVec_Rad(pos0.extract_row(step).transposed(), pos1.extract_row(1).transposed(), step);
		Matrix pos_radial10 = PosVec_Rad(pos1.extract_row(step).transposed(), pos0.extract_row(1).transposed(), step);
		Matrix tang0 = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield, w_tape);
		Matrix tang1 = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + t_tape, w_tape);
		Matrix tang_radial01 = TangLinVec_Rad(pos0.extract_row(step).transposed(), pos1.extract_row(1).transposed(), step);
		Matrix tang_radial10 = TangLinVec_Rad(pos1.extract_row(step).transposed(), pos0.extract_row(1).transposed(), step);

		inductance = 0;
		dist = 0;
		dotPro = 0;

		for (int i = 0; i < step; i++) {
			for (int j = 0; j < step; j++) {
				dist = sqrt(pow(pos0_base[i][0] - pos0[j][0], 2.) +
							pow(pos0_base[i][1] - pos0[j][1], 2.) +
							pow(pos0_base[i][2] - pos0[j][2], 2.));
				dotPro = tang0_base[i][0] * tang0[j][0] +
						 tang0_base[i][1] * tang0[j][1] +
						 tang0_base[i][2] * tang0[j][2];
				inductance += dotPro / dist;
			}
			//cout << "1 : " << i + 1 << "/" << step << endl;
		}

		for (int i = 0; i < step; i++) {
			for (int j = 0; j < step; j++) {
				dist = sqrt(pow(pos0_base[i][0] - pos1[j][0], 2.) +
							pow(pos0_base[i][1] - pos1[j][1], 2.) +
							pow(pos0_base[i][2] - pos1[j][2], 2.));
				dotPro = tang0_base[i][0] * tang1[j][0] +
						 tang0_base[i][1] * tang1[j][1] +
						 tang0_base[i][2] * tang1[j][2];
				inductance += dotPro / dist;
			}
			//cout << "2 : " << i + 1 << "/" << step << endl;
		}

		for (int i = 0; i < step; i++) {
			for (int j = 0; j < step; j++) {
				dist = sqrt(pow(pos0_base[i][0] - pos1[j][0], 2.) +
							pow(pos0_base[i][1] - pos1[j][1], 2.) +
							pow(pos0_base[i][2] - pos1[j][2], 2.));
				dotPro = tang0_base[i][0] * tang1[j][0] +
					tang0_base[i][1] * tang1[j][1] +
					tang0_base[i][2] * tang1[j][2];
				inductance += dotPro / dist;
			}
			//cout << "3 : " << i + 1 << "/" << step << endl;
		}

		for (int i = 0; i < step; i++) {
			for (int j = 0; j < step; j++) {
				dist = sqrt(pow(pos1_base[i][0] - pos0[j][0], 2.) +
							pow(pos1_base[i][1] - pos0[j][1], 2.) +
							pow(pos1_base[i][2] - pos0[j][2], 2.));
				dotPro = tang1_base[i][0] * tang0[j][0] +
						 tang1_base[i][1] * tang0[j][1] +
						 tang1_base[i][2] * tang0[j][2];
				inductance += dotPro / dist;
			}
			//cout << "3 : " << i + 1 << "/" << step << endl;
		}

		for (int i = 0; i < step; i++) {
			for (int j = 0; j < step; j++) {
				dist = sqrt(pow(pos1_base[i][0] - pos1[j][0], 2.) +
							pow(pos1_base[i][1] - pos1[j][1], 2.) +
							pow(pos1_base[i][2] - pos1[j][2], 2.));
				dotPro = tang1_base[i][0] * tang1[j][0] +
						 tang1_base[i][1] * tang1[j][1] +
						 tang1_base[i][2] * tang1[j][2];
				inductance += dotPro / dist;
			}
			//cout << "4 : " << i + 1 << "/" << step << endl;
		}

		cout << ring << "/" << ring_sp << " : " << mu0 / (4 * Pi) * inductance << " [H]" << endl;
		ofs_out << ring << "," << mu0 / (4 * Pi) * inductance << "\n";
	}

	/*
#pragma omp parallel
{
	int i, j;
	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (i = 0; i < step; i++) {
		for (j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_For[i][0] - pos2_For[j][0], 2.) +
						pow(pos1_For[i][1] - pos2_For[j][1], 2.) +
						pow(pos1_For[i][2] - pos2_For[j][2], 2.));
			dotPro = tang1_For[i][0] * tang2_For[j][0] +
					 tang1_For[i][1] * tang2_For[j][1] +
					 tang1_For[i][2] * tang2_For[j][2];
			inductance += dotPro / dist;
		}
		cout << "1 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;

	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_For[i][0] - pos2_Rev[j][0], 2.) +
						pow(pos1_For[i][1] - pos2_Rev[j][1], 2.) +
						pow(pos1_For[i][2] - pos2_Rev[j][2], 2.));
			dotPro = tang1_For[i][0] * tang2_Rev[j][0] +
					 tang1_For[i][1] * tang2_Rev[j][1] +
					 tang1_For[i][2] * tang2_Rev[j][2];
			inductance += dotPro / dist;
		}
		cout << "2 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;

	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_Rev[i][0] - pos2_For[j][0], 2.) +
						pow(pos1_Rev[i][1] - pos2_For[j][1], 2.) +
						pow(pos1_Rev[i][2] - pos2_For[j][2], 2.));
			dotPro = tang1_Rev[i][0] * tang2_For[j][0] +
					 tang1_Rev[i][1] * tang2_For[j][1] +
					 tang1_Rev[i][2] * tang2_For[j][2];
			inductance += dotPro / dist;
		}
		cout << "3 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;

	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_Rev[i][0] - pos2_Rev[j][0], 2.) +
						pow(pos1_Rev[i][1] - pos2_Rev[j][1], 2.) +
						pow(pos1_Rev[i][2] - pos2_Rev[j][2], 2.));
			dotPro = tang1_Rev[i][0] * tang2_Rev[j][0] +
					 tang1_Rev[i][1] * tang2_Rev[j][1] +
					 tang1_Rev[i][2] * tang2_Rev[j][2];
			inductance += dotPro / dist;
		}
		cout << "4 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;
}	*/

	processTime = clock() - startTime;
	cout << static_cast<double>(processTime) / 1000 << " [s]" << endl;

	system(gnuplot);

	return 0;
}

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列であることに注意！！

	Matrix position(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		position[i][0] = rad * cos(t);		// x成分
		position[i][1] = rad * sin(t);		// y成分
		position[i][2] = wid * t / (2 * Pi);	// z成分

		t += dh;
	}

	return position;
}

Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列であることに注意！！

	Matrix position(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		position[i][0] = rad * cos(t);		// x成分
		position[i][1] = rad * sin(t);		// y成分
		position[i][2] = wid * -t / (2 * Pi);	// z成分

		t += dh;
	}

	return position;
}

Matrix PosVec_Rad(IN Matrix init, IN Matrix end, IN long n) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// 戻り値はn行3列の配列であることに注意！！

	if (init.row_size() != 3 || init.column_size() != 1) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : PosVec_Rad" << endl;
		exit(1);
	}
	if (end.row_size() != 3 || end.column_size() != 1) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : PosVec_Rad" << endl;
		exit(1);
	}

	Matrix position(n, 3);
	Matrix dh(3);

	*dh[0] = (*end[0] - *init[0]) / n;
	*dh[1] = (*end[1] - *init[1]) / n;
	*dh[2] = (*end[2] - *init[2]) / n;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		position[i][0] = *init[0] + *dh[0] * i;
		position[i][1] = *init[1] + *dh[1] * i;
		position[i][2] = *init[2] + *dh[2] * i;
	}

	return position;
}

Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列であることに注意！！

	Matrix tangent_line(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、微小接線ベクトルを計算
		tangent_line[i][0] = -dh * rad * sin(t);	// x成分
		tangent_line[i][1] = dh * rad * cos(t);	// y成分
		tangent_line[i][2] = dh * wid / (2 * Pi);// z成分

		t += dh;
	}

	return tangent_line;
}

Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列であることに注意！！

	Matrix tangent_line(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、微小接線ベクトルを計算
		tangent_line[i][0] = -dh * rad * sin(t);	// x成分
		tangent_line[i][1] = dh * rad * cos(t);	// y成分
		tangent_line[i][2] = dh * wid * -1 / (2 * Pi);// z成分

		t += dh;
	}

	return tangent_line;
}

Matrix TangLinVec_Rad(IN Matrix init, IN Matrix end, IN long n) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// 戻り値は3行1列の配列であることに注意！！

	if (init.row_size() != 3 || init.column_size() != 1) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : TangLinVec_Rad" << endl;
		exit(1);
	}
	if (end.row_size() != 3 || end.column_size() != 1) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : TangLinVec_Rad" << endl;
		exit(1);
	}

	Matrix tangent_line(3);

	*tangent_line[0] = (*end[0] - *init[0]) / n;
	*tangent_line[1] = (*end[1] - *init[1]) / n;
	*tangent_line[2] = (*end[2] - *init[2]) / n;

	return tangent_line;
}
