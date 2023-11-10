#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "matrix.h"

#define IN
#define OUT
#define Pi 3.14159265
#define mu0 (4 * Pi * 1e-7)

#define r_shield 0.015
#define w_tape 0.012
#define t_tape 0.1e-3
#define n_layer 4
#define n_turn 2
#define n_ring (n_turn * 2)

#define step 500

using namespace std;

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix PosVec_Rad(IN Matrix init, IN Matrix end, IN long n);
Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_Rad(IN Matrix init, IN Matrix end, IN long n);

int main() {
	// 時間計測開始
	clock_t startTime, processTime;
	startTime = clock();

	// ファイル書き込みオープン
	ofstream csv_out[n_layer - 1];
	for (int i = 0; i < n_layer - 1; i++) {
		ostringstream file_name;
		file_name << "inductance_layer_base_" << i << ".csv";
		new(&csv_out[i]) ofstream(file_name.str());
	}

	// インダクタンスの入れ物用意
	Matrix ind[n_layer - 1];
	for (int i = 0; i < n_layer - 1; i++) {
		new(&ind[i]) Matrix(n_ring, n_layer - 1);
	}

	// インダクタンス計算
	for (int layer_base = 0; layer_base < n_layer - 1; layer_base++) {
		Matrix pos_base_in(step, 3);
		Matrix pos_base_out(step, 3);
		Matrix tang_base_in(step, 3);
		Matrix tang_base_out(step, 3);
		if (layer_base % 2 == 0) {
			pos_base_in = PosVec_For(0, Pi, step, r_shield + layer_base * t_tape, w_tape, 0);
			pos_base_out = PosVec_Rev(-Pi, 0, step, r_shield + (layer_base + 1) * t_tape, w_tape, 0);
			tang_base_in = TangLinVec_For(0, Pi, step, r_shield + layer_base * t_tape, w_tape, 0);
			tang_base_out = TangLinVec_Rev(-Pi, 0, step, r_shield + (layer_base + 1) * t_tape, w_tape, 0);
		}
		else {
			pos_base_in = PosVec_Rev(-Pi, 0, step, r_shield + layer_base * t_tape, w_tape, 0);
			pos_base_out = PosVec_For(0, Pi, step, r_shield + (layer_base + 1) * t_tape, w_tape, 0);
			tang_base_in = TangLinVec_Rev(-Pi, 0, step, r_shield + layer_base * t_tape, w_tape, 0);
			tang_base_out = TangLinVec_For(0, Pi, step, r_shield + (layer_base + 1) * t_tape, w_tape, 0);
		}

		for (int layer = 0; layer < n_layer - 1; layer++) {
			for (int ring = 0; ring < n_ring; ring++) {
				double inductance = 0;
				double dist = 0;
				double dotPro = 0;
				double thickness;

				Matrix pos_in(step, 3);
				Matrix pos_out(step, 3);
				Matrix tang_in(step, 3);
				Matrix tang_out(step, 3);

				thickness = 0;
				if (layer == layer_base && ring == 0) thickness = t_tape;
				if (layer % 2 == 0) {
					pos_in = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
					tang_in = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
				}
				else {
					pos_in = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
					tang_in = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
				}
				for (int i = 0; i < step; i++) {
					for (int j = 0; j < step; j++) {
						dist = sqrt(pow(pos_base_in[i][0] - pos_in[j][0], 2.) +
									pow(pos_base_in[i][1] - pos_in[j][1], 2.) +
									pow(pos_base_in[i][2] - pos_in[j][2], 2.));
						dotPro = tang_base_in[i][0] * tang_in[j][0] +
								 tang_base_in[i][1] * tang_in[j][1] +
								 tang_base_in[i][2] * tang_in[j][2];
						inductance += dotPro / dist;
					}
				}

				thickness = 0;
				if (layer == layer_base - 1 && ring == 0) thickness = t_tape;
				if (layer % 2 == 0) {
					pos_out = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
					tang_out = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
				}
				else {
					pos_out = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
					tang_out = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
				}
				for (int i = 0; i < step; i++) {
					for (int j = 0; j < step; j++) {
						dist = sqrt(pow(pos_base_in[i][0] - pos_out[j][0], 2.) +
									pow(pos_base_in[i][1] - pos_out[j][1], 2.) +
									pow(pos_base_in[i][2] - pos_out[j][2], 2.));
						dotPro = tang_base_in[i][0] * tang_out[j][0] +
								 tang_base_in[i][1] * tang_out[j][1] +
								 tang_base_in[i][2] * tang_out[j][2];
						inductance += dotPro / dist;
					}
				}

				thickness = 0;
				if (layer == layer_base + 1 && ring == 0) thickness = t_tape;
				if (layer % 2 == 0) {
					pos_in = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
					tang_in = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
				}
				else {
					pos_in = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
					tang_in = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
				}
				for (int i = 0; i < step; i++) {
					for (int j = 0; j < step; j++) {
						dist = sqrt(pow(pos_base_out[i][0] - pos_in[j][0], 2.) +
									pow(pos_base_out[i][1] - pos_in[j][1], 2.) +
									pow(pos_base_out[i][2] - pos_in[j][2], 2.));
						dotPro = tang_base_out[i][0] * tang_in[j][0] +
								 tang_base_out[i][1] * tang_in[j][1] +
								 tang_base_out[i][2] * tang_in[j][2];
						inductance += dotPro / dist;
					}
				}

				thickness = 0;
				if (layer == layer_base && ring == 0) thickness = t_tape;
				if (layer % 2 == 0) {
					pos_out = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
					tang_out = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
				}
				else {
					pos_out = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
					tang_out = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + (layer + 1) * t_tape, w_tape, thickness);
				}
				for (int i = 0; i < step; i++) {
					for (int j = 0; j < step; j++) {
						dist = sqrt(pow(pos_base_out[i][0] - pos_out[j][0], 2.) +
									pow(pos_base_out[i][1] - pos_out[j][1], 2.) +
									pow(pos_base_out[i][2] - pos_out[j][2], 2.));
						dotPro = tang_base_out[i][0] * tang_out[j][0] +
								 tang_base_out[i][1] * tang_out[j][1] +
								 tang_base_out[i][2] * tang_out[j][2];
						inductance += dotPro / dist;
					}
				}

				ind[layer_base][ring][layer] = mu0 / (4 * Pi) * inductance;
				cout << layer_base << "," << ring << ", " << layer << endl;
			}
		}
	}

	// 時間計測終了・表示
	processTime = clock() - startTime;
	cout << static_cast<double>(processTime) / 1000 << " [s]" << endl;

	// 確認等用csv書き出し
	{
		for (int layer_base = 0; layer_base < n_layer - 1; layer_base++) {
			for (int ring = 0; ring < n_ring; ring++) {
				for (int layer = 0; layer < n_layer - 1; layer++) {
					csv_out[layer_base] << ind[layer_base][ring][layer];
					if (layer != n_layer - 2)csv_out[layer_base] << ",";
				}
				csv_out[layer_base] << endl;
			}
		}
	}

	// ファイル書き込みクローズ・gnuplot書き出し
	for (int i = 0; i < n_layer - 1; i++) {
		csv_out[i].close();
	}
	//FILE* gnuplot = _popen("gnuplot -persist", "w");
	//fprintf(gnuplot, "load 'setting.plt'\n");
	//_pclose(gnuplot);

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
