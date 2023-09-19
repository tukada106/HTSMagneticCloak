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

void PosVec(IN double init, IN double end, IN long n, IN double rad, IN double wid, OUT double** arr);
void TangLinVec(IN double init, IN double end, IN long n, IN double rad, IN double wid, OUT double** arr);

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	// 位置ベクトル格納領域の確保
	double** dpPos = nullptr;
	dpPos = (double**)malloc(step * sizeof(double*));
	if (dpPos == nullptr) {
		cout << "メモリ確保失敗\n";
		return -1;
	}
	for (int i = 0; i < step; i++) {
		dpPos[i] = (double*)malloc(3 * sizeof(double));
		if (dpPos[i] == nullptr) {
			cout << "メモリ確保失敗";
			return -1;
		}
	}

	PosVec(sec_start, sec_stop, step, r_shield, w_tape, dpPos);
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < 3; j++) {
			cout << dpPos[i][j] << "\t";
		}
		cout << "\n";
	}

	string str_buf;
	string str_conma_buf;
	ofstream ofs_out("output.csv");

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < 6; j++) {
			ofs_out << dpPos[i][j] << ",";
		}
		ofs_out << "\n";
	}

	// ベクトル格納領域の解放
	for (int i = 0; i < step; i++) {
		free(dpPos[i]);
	}
	free(dpPos);

	system(gnuplot);

	return 0;
}

void PosVec(IN double init, IN double end, IN long n, IN double rad, IN double wid, OUT double** arr) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// arr	:戻り値の配列
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列の配列を渡すこと！！

	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		arr[i][0] = rad * cos(t);		// x成分
		arr[i][1] = rad * sin(t);		// y成分
		arr[i][2] = wid * t / (2 * Pi);	// z成分

		t += dh;
	}
}

void TangLinVec(IN double init, IN double end, IN long n, IN double rad, IN double wid, OUT double** arr) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// arr:戻り値の配列
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// 戻り値はn行3列の配列の配列を渡すこと！！

	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、微小接線ベクトルを計算
		arr[i][0] = -dh * rad * sin(t);	// x成分
		arr[i][1] = dh * rad * cos(t);	// y成分
		arr[i][2] = dh * wid / (2 * Pi);// z成分

		t += dh;
	}
}
