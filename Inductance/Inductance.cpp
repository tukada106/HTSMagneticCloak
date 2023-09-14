#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#define gnuplot "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist setting.plt"

#define IN
#define OUT
#define Pi 3.14159265

#define r_shield 0.015
#define w_tape 0.012
#define step 200

using namespace std;

void curveVector(IN double init, IN double end, IN long n, IN double rad, IN double wid, OUT double arr[][6]);

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	// ベクトル情報の取得
	long n = step;
	double vector[step][6];

	curveVector(0, 20 * Pi, n, r_shield, w_tape, vector);
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < 6; j++) {
			cout << vector[i][j] << "\t";
		}
		cout << "\n";
	}

	string str_buf;
	string str_conma_buf;
	ofstream ofs_out("output.csv");

	for (int i = 0; i < step; i++) {
		for (int j = 0; j < 6; j++) {
			ofs_out << vector[i][j] << ",";
		}
		ofs_out << "\n";
	}

	system(gnuplot);

	return 0;
}

void curveVector(IN double init, IN double end, IN long n, IN double rad, IN double wid, OUT double arr[][6]) {
	// init:区間の始点
	// end:区間の終点
	// n:分割数
	// arr:戻り値の配列
	// 戻り値はn行6列の配列の配列を渡すこと！！

	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		arr[i][0] = r_shield * cos(t);		// x成分
		arr[i][1] = r_shield * sin(t);		// y成分
		arr[i][2] = w_tape * t / (2 * Pi);	// z成分

		// 以下、微小接線ベクトルを計算
		arr[i][3] = -dh * r_shield * sin(t);// x成分
		arr[i][4] = dh * r_shield * cos(t);	// y成分
		arr[i][5] = dh * w_tape / (2 * Pi);	// z成分

		t += dh;
	}
}
