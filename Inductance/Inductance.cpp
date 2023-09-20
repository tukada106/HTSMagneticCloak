#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "matrix.h"

#define gnuplot "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist setting.plt"

#define IN
#define OUT
#define Pi 3.14159265

#define r_shield 0.015
#define w_tape 0.012
#define sec_start 0
#define sec_stop 20 * Pi
#define step 200

using namespace std;

Matrix PosVec(IN double init, IN double end, IN long n, IN double rad, IN double wid);
Matrix TangLinVec(IN double init, IN double end, IN long n, IN double rad, IN double wid);

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	Matrix test = PosVec(sec_start, sec_stop, step, r_shield, w_tape);
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < 3; j++) {
			cout << test[i][j] << "\t";
		}
		cout << "\n";
	}

	cout << "\n";
	cout << "\n";
	Matrix vec1(3);
	Matrix vec2(3);
	Matrix ans(3);
	double answ = 0;
	vec1[0][0] = 1;
	vec1[1][0] = 8;
	vec1[2][0] = 6;
	vec2[0][0] = 5;
	vec2[1][0] = 10;
	vec2[2][0] = 7;
	ans = cross(vec1, vec2);
	answ = dot(vec1, vec2);
	for (int i = 0; i < ans.row_size(); i++) {
		for (int j = 0; j < ans.column_size(); j++) {
			cout << ans[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << answ << endl;

	string str_buf;
	string str_conma_buf;
	ofstream ofs_out("output.csv");

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < 3; j++) {
			ofs_out << test[i][j] << ",";
		}
		ofs_out << "\n";
	}

	system(gnuplot);

	return 0;
}

Matrix PosVec(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
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
		position[i][2] = wid * t / (2 * Pi);	// z成分

		t += dh;
	}

	return position;
}

Matrix TangLinVec(IN double init, IN double end, IN long n, IN double rad, IN double wid) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// arr:戻り値の配列
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列であることに注意！！

	Matrix tangent_line;
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
