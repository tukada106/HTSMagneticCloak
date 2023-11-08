#include <iostream>
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
#define n_turn 10
#define n_ring (n_turn * 2)
#define n_loop (n_layer - 1) * n_ring

#define B_apply 8.48e-3
#define t_sweep 0.5e-3
#define R_contact 100e-6

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

	// 抵抗行列の作成
	Matrix mat_res(n_loop, n_loop);
	for (int mat_res_row = 0; mat_res_row < n_loop; mat_res_row++) {
		for (int mat_res_col = 0; mat_res_col < n_loop; mat_res_col++) {
			if (mat_res_row == mat_res_col) mat_res[mat_res_row][mat_res_col] = -2 * R_contact;
			if (mat_res_col == mat_res_row + n_layer - 1) mat_res[mat_res_row][mat_res_col] = R_contact;
			if (mat_res_col == mat_res_row - (n_layer - 1)) mat_res[mat_res_row][mat_res_col] = R_contact;
		}
	}

	// 4次ルンゲクッタ
	Matrix vec_current(n_loop);
	Matrix vec_alpha(n_loop);
	Matrix vec_result(n_loop + 1);
	Matrix vec_K0(n_loop);
	Matrix vec_K1(n_loop);
	Matrix vec_K2(n_loop);
	Matrix vec_K3(n_loop);

	// ファイル書き込みオープン
	ofstream csv_out[2];
	for (int i = 0; i < 2; i++) {
		ostringstream file_name;
		file_name << "matrix_" << i << ".csv";
		new(&csv_out[i]) ofstream(file_name.str());
	}

	// 確認等用csv書き出し
	{
		for (int i = 0; i < 2; i++) {
			for (int mat_row = 0; mat_row < n_loop; mat_row++) {
				for (int mat_col = 0; mat_col < n_loop; mat_col++) {
					if (i == 0) {
						csv_out[i] << mat_ind[mat_row][mat_col];
					}
					else {
						csv_out[i] << mat_res[mat_row][mat_col];
					}
					if (mat_col != n_loop - 1) csv_out[i] << ",";
				}
				csv_out[i] << endl;
			}
		}
	}

	// ファイル書き込みクローズ・gnuplot書き出し
	for (int i = 0; i < 2; i++) {
		csv_out[i].close();
	}

	// 時間計測終了・表示
	processTime = clock() - startTime;
	cout << static_cast<double>(processTime) / 1000 << " [s]" << endl;
}
