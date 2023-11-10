﻿#include <iostream>
#include <fstream>
#include <sstream>

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
#define n_loop (n_layer - 1) * n_ring

#define B_apply 8.48e-3
#define t_sweep 0.5e-3
#define R_contact 10e-6

#define t_init 0
#define t_end 1e-3
#define dh_max 1e-6
#define dh_min 1e-9
#define interval 10
#define tolerance 1000

using namespace std;

int main() {
	// 時間計測開始
	clock_t startTime, processTime;
	startTime = clock();

	// インダクタンスの入れ物用意
	Matrix ind[n_layer - 1];
	for (int i = 0; i < n_layer - 1; i++) {
		new(&ind[i]) Matrix(n_ring, n_layer - 1);
	}

	// ファイル読み取りオープン
	ifstream ifs_csv_in[n_layer - 1];
	for (int i = 0; i < n_layer - 1; i++) {
		ostringstream oss_file_name;
		oss_file_name << "../Inductance/";
		oss_file_name << "inductance_layer_base_" << i << ".csv";
		new(&ifs_csv_in[i]) ifstream(oss_file_name.str());
	}

	// インダクタンス読み込み
	for (int layer_base = 0; layer_base < n_layer - 1; layer_base++) {
		string str_csv_in_line;
		string str_csv_in_comma;
		for (int ring = 0; ring < n_ring; ring++) {
			if (getline(ifs_csv_in[layer_base], str_csv_in_line)) {
				istringstream iss_csv_in_line(str_csv_in_line);
				for (int layer = 0; layer < n_layer - 1; layer++) {
					if (getline(iss_csv_in_line, str_csv_in_comma, ',')) {
						ind[layer_base][ring][layer] = stod(str_csv_in_comma);
					}
					else {
						cerr << "The specified ring is missing!" << endl;
						cerr << "getline(iss_csv_in_line, str_csv_in_comma, ',')" << endl;
						exit(1);
					}
				}
			}
			else {
				cerr << "The specified layer_base is missing!" << endl;
				cerr << "getline(ifs_csv_in[layer_base], str_csv_in_line)" << endl;
				exit(1);
			}
		}
	}

	// ファイル読み取りクローズ
	for (int i = 0; i < n_layer - 1; i++) {
		ifs_csv_in[i].close();
	}

	// インダクタンス行列の作成
	Matrix mat_ind(n_loop, n_loop);
	for (int mat_ind_row = 0; mat_ind_row < n_loop; mat_ind_row++) {
		static int second_row = 0;
		if (mat_ind_row % 3 == 0 && mat_ind_row != 0) second_row -= 1;
		int second_col = 0;
		for (int mat_ind_col = 0; mat_ind_col < n_loop; mat_ind_col++) {
			int first = mat_ind_row % (n_layer - 1);
			if (mat_ind_col % 3 == 0 && mat_ind_col != 0) second_col += 1;
			int second = abs(second_row + second_col);
			int third = mat_ind_col % (n_layer - 1);
			mat_ind[mat_ind_row][mat_ind_col] = ind[first][second][third];
		}
	}

	// インダクタンス行列の逆行列を用意
	Matrix mat_ind_inverse = mat_ind.inverse();

	// 抵抗行列の作成
	Matrix mat_res(n_loop, n_loop);
	for (int mat_res_row = 0; mat_res_row < n_loop; mat_res_row++) {
		for (int mat_res_col = 0; mat_res_col < n_loop; mat_res_col++) {
			if (mat_res_row == mat_res_col) mat_res[mat_res_row][mat_res_col] = -2 * R_contact;
			if (mat_res_col == mat_res_row + n_layer - 1) mat_res[mat_res_row][mat_res_col] = R_contact;
			if (mat_res_col == mat_res_row - (n_layer - 1)) mat_res[mat_res_row][mat_res_col] = R_contact;
		}
	}

	// RK法　ベクトル準備
	Matrix vec_current_4th(n_loop);
	Matrix vec_current_5th(n_loop);
	Matrix vec_alpha(n_loop);
	Matrix vec_result(n_loop + 1);
	Matrix vec_K0(n_loop);
	Matrix vec_K1(n_loop);
	Matrix vec_K2(n_loop);
	Matrix vec_K3(n_loop);
	Matrix vec_K4(n_loop);
	Matrix vec_K5(n_loop);
	Matrix vec_K6(n_loop);
	Matrix vec_error(n_loop);
	Matrix vec_delta(n_loop);

	// RK法　電流ベクトル初期値代入
	for (int loop = 0; loop < n_loop; loop++) {
		*vec_current_4th[loop] = 0.;
		*vec_current_5th[loop] = 0.;
	}
	double t = t_init;
	double dh = dh_max;
	bool flag_calculate = true;

	// RK法　計算結果　ファイル書き込みオープン
	ofstream csv_out_result("result_current.csv");

	// RK法　RK-DPで計算
	while (flag_calculate) {
		static int count = 0;
		//cout << "\r" << static_cast<double>((t - t_init) / (t_end - t_init)) * 100. << "%" << "\t\t" << count++;
		cout << "t:" << t << "[s], " << dh << endl;

		if (t <= t_sweep) {
			for (int loop = 0; loop < n_loop; loop++) {
				*vec_alpha[loop] = 1 * B_apply * Pi * pow(r_shield, 2.) / t_sweep;
			}
		}
		else {
			for (int loop = 0; loop < n_loop; loop++) {
				*vec_alpha[loop] = 0;
			}
		}

		vec_K0 = mat_ind_inverse * (vec_alpha + mat_res * vec_current_4th);
		vec_K1 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_4th + dh * (1. / 5.) * vec_K0));
		vec_K2 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_4th + dh * ((3. / 40.) * vec_K0 +
																				   (9. / 40.) * vec_K1)));
		vec_K3 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_4th + dh * ((44. / 45.) * vec_K0 +
																				  (-56. / 15.) * vec_K1 +
																					(32. / 9.) * vec_K2)));
		vec_K4 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_4th + dh * ((19372. / 6561.) * vec_K0 +
																				  (-25360. / 2187.) * vec_K1 +
																				   (64448. / 6561.) * vec_K2 +
																					 (-212. / 729.) * vec_K3)));
		vec_K5 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_4th + dh * ((9017. / 3168.) * vec_K0 +
																					 (-355. / 33.) * vec_K1 +
																				  (46732. / 5247.) * vec_K2 +
																					 (-49. / 176.) * vec_K3 +
																				 (-5103. / 18656.) * vec_K4)));
		vec_current_4th = vec_current_4th + dh * ((35. / 384.) * vec_K0 +
														//(0.) * vec_K1 +
												(500. / 1113.) * vec_K2 +
												 (125. / 192.) * vec_K3 +
											  (-2187. / 6784.) * vec_K4 +
												   (11. / 84.) * vec_K5);
		vec_K6 = mat_ind_inverse * (vec_alpha + mat_res * vec_current_4th);
		vec_current_5th = vec_current_5th + dh * ((5179. / 57600.) * vec_K0 +
															//(0.) * vec_K1 +
												  (7571. / 16695.) * vec_K2 +
													 (393. / 640.) * vec_K3 +
											   (-92097. / 339200.) * vec_K4 +
													(187. / 2100.) * vec_K5 +
														(1. / 40.) * vec_K6);

		for (int loop = 0; loop < n_loop; loop++) {
			*vec_error[loop] = abs(*vec_current_4th[loop] - *vec_current_5th[loop]);
		}

		bool output = true;
		for (int loop = 0; loop < n_loop; loop++) {
			if (*vec_error[loop] > tolerance) {
				output = false;
			}
		}
		if (output) {
			t += dh;
			static int count_interval = 0;
			if (count_interval++ % interval == 0) {
				csv_out_result << t << ",";
				for (int loop = 0; loop < n_loop; loop++) {
					csv_out_result << *vec_current_4th[loop];
					if (loop != n_loop - 1) {
						csv_out_result << ",";
					}
					else {
						csv_out_result << endl;
					}
				}
			}
		}

		double delta_min = 1e12;
		for (int loop = 0; loop < n_loop; loop++) {
			*vec_delta[loop] = pow(tolerance * dh / (2 * *vec_error[loop]), 1. / 5.);
			if (*vec_delta[loop] < delta_min) {
				delta_min = *vec_delta[loop];
			}
		}

		dh = dh * delta_min;
		//if (delta_min <= 0.1) {
		//	dh *= 0.1;
		//}
		//else if (delta_min >= 4.) {
		//	dh *= 4;
		//}
		//else {
		//	dh *= delta_min;
		//}

		if (dh > dh_max) {
			dh = dh_max;
		}

		if (t >= t_end) {
			flag_calculate = false;
		}
		else if (t + dh > t_end) {
			dh = t_end - t;
			if (dh < dh_min) {
				flag_calculate = false;
				cout << "Minimum dh exceeded!" << endl;
			}
		}
	}

	// RK法　計算結果　ファイル書き込みクローズ
	csv_out_result.close();

	// 時間計測終了・表示
	processTime = clock() - startTime;
	cout << endl << static_cast<double>(processTime) / 1000 << " [s]" << endl;
}
