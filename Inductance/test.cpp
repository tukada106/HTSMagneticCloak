#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "matrix.h"

#define gnuplot "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist setting.plt"

#define IN
#define OUT
#define Pi 3.14159265
#define mu0 (4 * Pi * 1e-7)

#define r_shield 0.015
#define w_tape 0.012
#define sec_start 0
#define sec_stop (20 * Pi)
#define step 2000

using namespace std;

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid);

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	Matrix pos1_For = PosVec_For(0, 2 * Pi, step, 0.01, w_tape);
	Matrix pos1_Rev = PosVec_Rev(0, 2 * Pi, step, 0.01, w_tape);
	Matrix pos2_For = PosVec_For(0, 2 * Pi, step, 0.01, w_tape);
	Matrix pos2_Rev = PosVec_Rev(0, 2 * Pi, step, 0.01, w_tape);
	Matrix tang1_For = TangLinVec_For(0, 2 * Pi, step, 0.01, w_tape);
	Matrix tang1_Rev = TangLinVec_Rev(0, 2 * Pi, step, 0.01, w_tape);
	Matrix tang2_For = TangLinVec_For(0, 2 * Pi, step, 0.01, w_tape);
	Matrix tang2_Rev = TangLinVec_Rev(0, 2 * Pi, step, 0.01, w_tape);
	double inductance = 0;
	double dist = 0;
	double dotPro = 0;

/*	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_For[i][0] - pos2_For[j][0], 2.) +
						pow(pos1_For[i][1] - pos2_For[j][1], 2.) +
						pow(pos1_For[i][2] - pos2_For[j][2], 2.));
			dotPro = tang1_For[i][0] * tang2_For[j][0] +
					 tang1_For[i][1] * tang2_For[j][1] +
					 tang1_For[i][2] * tang2_For[j][2];
			inductance += dotPro / dist;
		}
	}
	*/
	//for (int i = 0; i < step; i++) {
	//	for (int j = 0; j < step; j++) {
	//		dist = sqrt(pow(pos1_For[i][0] - pos2_Rev[j][0], 2.) +
	//					pow(pos1_For[i][1] - pos2_Rev[j][1], 2.) +
	//					pow(pos1_For[i][2] - pos2_Rev[j][2], 2.));
	//		dotPro = tang1_For[i][0] * tang2_Rev[j][0] +
	//				 tang1_For[i][1] * tang2_Rev[j][1] +
	//				 tang1_For[i][2] * tang2_Rev[j][2];
	//		inductance += dotPro / dist;
	//	}
	//}
	/*
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
	}
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
	}
	*/
	cout << mu0 / (4 * Pi) * inductance << " [H]" << endl;

	string str_buf;
	string str_conma_buf;
	ofstream ofs_out("output.csv");

	for (int i = 0; i < step; i++) {
		ofs_out << pos1_For[i][0] << ",";
		ofs_out << pos1_For[i][1] << ",";
		ofs_out << pos1_For[i][2] << ",";
		ofs_out << pos1_For[i][0] << ",";
		ofs_out << pos1_For[i][1] << ",";
		ofs_out << pos1_For[i][2];
		ofs_out << "\n";
	}

	processTime = clock() - startTime;
	cout << static_cast<double>(processTime) / 1000 << " [s]" << endl;

	system(gnuplot);

	return 0;
}

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// arr	:戻り値の配列
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
		position[i][2] = 0.5;	// z成分

		t += dh;
	}

	return position;
}

Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// arr	:戻り値の配列
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
		position[i][2] = -0.5;	// z成分

		t += dh;
	}

	return position;
}

Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// arr:戻り値の配列
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
		tangent_line[i][2] = 0;// z成分

		t += dh;
	}

	return tangent_line;
}

Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// arr:戻り値の配列
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
		tangent_line[i][2] = 0;// z成分

		t += dh;
	}

	return tangent_line;
}
