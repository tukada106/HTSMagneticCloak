#include <iostream>
#include "matrix.h"

#define IN
#define OUT
#define Pi 3.14159265
#define mu0 (4 * Pi * 1e-7)

#define r_shield 0.015
#define w_tape 0.012
#define t_tape 0.1e-3
#define n_layer 2	// 2layers' condition
#define n_turn 25	// 2layers' condition
#define n_ring (n_turn * 2)
#define n_loop (n_layer - 1) * n_ring

using namespace std;

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);

int main() {
	cout << "hello!" << endl;

	return 0;
}

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// thick:線材厚さ（同一経路を積分する場合に設定、それ以外は0にする）
	// 戻り値はn行3列の配列であることに注意！！

	Matrix position(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		position[i][0] = (rad + thick / 2) * cos(t);	// x成分
		position[i][1] = (rad + thick / 2) * sin(t);	// y成分
		position[i][2] = wid * t / (2 * Pi);			// z成分

		if (thick != 0.) position[i][2] += wid / 2;

		t += dh;
	}

	return position;
}

Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// thick:線材厚さ（同一経路を積分する場合に設定、それ以外は0にする）
	// 戻り値はn行3列の配列であることに注意！！

	Matrix position(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、位置ベクトルを計算
		position[i][0] = (rad + thick / 2) * cos(t);	// x成分
		position[i][1] = (rad + thick / 2) * sin(t);	// y成分
		position[i][2] = wid * -t / (2 * Pi);			// z成分

		if (thick != 0.) position[i][2] += wid / 2;

		t += dh;
	}

	return position;
}

Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// thick:線材厚さ（同一経路を積分する場合に設定、それ以外は0にする）
	// 戻り値はn行3列の配列であることに注意！！

	Matrix tangent_line(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、微小接線ベクトルを計算
		tangent_line[i][0] = -dh * (rad + thick / 2) * sin(t);	// x成分
		tangent_line[i][1] = dh * (rad + thick / 2) * cos(t);	// y成分
		tangent_line[i][2] = dh * wid / (2 * Pi);				// z成分

		t += dh;
	}

	return tangent_line;
}

Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick) {
	// init	:区間の始点
	// end	:区間の終点
	// n	:分割数
	// rad	:螺旋の半径
	// wid	:螺旋の上がり幅
	// thick:線材厚さ（同一経路を積分する場合に設定、それ以外は0にする）
	// 戻り値はn行3列の配列であることに注意！！

	Matrix tangent_line(n, 3);
	double dh, t;

	dh = (end - init) / n;
	t = init;

	for (int i = 0; i < n; i++) {
		// 以下、微小接線ベクトルを計算
		tangent_line[i][0] = -dh * (rad + thick / 2) * sin(t);	// x成分
		tangent_line[i][1] = dh * (rad + thick / 2) * cos(t);	// y成分
		tangent_line[i][2] = dh * wid * -1 / (2 * Pi);			// z成分

		t += dh;
	}

	return tangent_line;
}
