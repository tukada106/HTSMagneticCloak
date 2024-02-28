#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <queue>
#include <omp.h>

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

#define B_apply 8.48e-3
#define t_sweep 0.5e-3
#define R_contact 4.99e-6 // or 7.36e-6 // 2layers' condition

#define t_init 0
#define t_end 2	// 2layers' condition
#define dh_init 1e-12
#define dh_max 0.001
#define dh_min 1e-12
#define interval 10
#define tolerance 1e-2
#define que_max 0

using namespace std;

int RKDPupdate(Matrix& current_4th, Matrix& current_5th,
			   Matrix& const inductance_inv, Matrix& const resistance, Matrix& const alpha, double dh,
			   int thread_num, int n_threads);

int main() {
	// Matrixクラスの足し算掛け算並列化の検証
	clock_t old;
	const int row_max_test = 5000;
	Matrix a(row_max_test, row_max_test);
	cout << "0" << endl;
	Matrix b(row_max_test, row_max_test);
	cout << "1" << endl;
	double c = 3.;
	Matrix calc(row_max_test, row_max_test);
	cout << "2" << endl;
	Matrix correct(row_max_test, row_max_test);
	cout << "3" << endl;
	for (int i = 0; i < row_max_test; i++) {
		for (int j = 0; j < row_max_test; j++) {
			a[i][j] = 1;
			b[i][j] = 15;
		}
	}
	cout << "4" << endl;
	old = clock();
	correct = a * b;
	cout << "5 " << clock() - old << "[ms]" << endl;
	//omp_set_num_threads(8);
	old = clock();
#pragma omp parallel
{
	omp_multi(calc, a, b, omp_get_thread_num(), omp_get_num_threads());
}
	cout << "6 " << clock() - old << "[ms]" << endl;
	for (int i = 0; i < row_max_test; i++) {
		for (int j = 0; j < row_max_test; j++) {
			if (calc[i][j] != correct[i][j]) {
				cout << "calc[" << i << "][" << j << "] = " << calc[i][j] << endl;
			}
		}
	}
	return 0;

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
		oss_file_name << "master_thesis_condition_test/";	// 2layers' condition
		oss_file_name << "2layers.csv";						// 2layers' condition
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
		if (mat_ind_row % (n_layer - 1) == 0 && mat_ind_row != 0) second_row -= 1;
		int second_col = 0;
		for (int mat_ind_col = 0; mat_ind_col < n_loop; mat_ind_col++) {
			int first = mat_ind_row % (n_layer - 1);
			if (mat_ind_col % (n_layer - 1) == 0 && mat_ind_col != 0) second_col += 1;
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

	// 確認等用csv書き出し
	{
		ofstream csv_out("matrix_test.csv");
		for (int row = 0; row < n_loop; row++) {
			for (int col = 0; col < n_loop; col++) {
				csv_out << scientific << setprecision(15) << uppercase;
				csv_out << mat_ind[row][col];
				if (col != n_loop - 1)csv_out << ",";
			}
			csv_out << endl;
		}
	}

	// RK法　ベクトル準備
	Matrix vec_current_4th(n_loop);
	Matrix vec_current_5th(n_loop);
	Matrix vec_K0(n_loop);
	Matrix vec_K1(n_loop);
	Matrix vec_K2(n_loop);
	Matrix vec_K3(n_loop);
	Matrix vec_K4(n_loop);
	Matrix vec_K5(n_loop);
	Matrix vec_K6(n_loop);
	Matrix vec_alpha(n_loop);
	Matrix vec_error(n_loop);

	// RK法　結果ロールバック用キュー準備
	queue<double> que_t_old;
	queue<Matrix> que_current_5th_old;

	// RK法　電流ベクトル、時間、刻み幅初期値代入
	// （電流は基本的に0でOK、デバッグ時に0以外に変更）
	for (int loop = 0; loop < n_loop; loop++) {
		*vec_current_4th[loop] = 0.;
		*vec_current_5th[loop] = 0.;
	}
	double t = t_init;
	double dh = dh_init;
	bool flag_calculate = true;

	// RK法　計算結果　ファイル書き込みオープン
	ofstream csv_out_result("result_current.csv");
	csv_out_result << scientific << setprecision(15) << uppercase;

	// RK法　RK-DPで計算
	while (flag_calculate) {
		// コンソールに解析内の時間を表示
		cout << "t:" << t << "[s],\t" << dh << "\t->\t";

		// 外部磁場の印加条件を定義
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

		// 直前の計算結果を一時退避
		que_t_old.push(t);
		que_current_5th_old.push(vec_current_5th);

		// RK-DP法の係数を計算
		vec_K0 = mat_ind_inverse * (vec_alpha + mat_res * vec_current_5th);
		vec_K1 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_5th + dh * (1. / 5.) * vec_K0));
		vec_K2 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_5th + dh * ((3. / 40.) * vec_K0 +
																				   (9. / 40.) * vec_K1)));
		vec_K3 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_5th + dh * ((44. / 45.) * vec_K0 +
																				  (-56. / 15.) * vec_K1 +
																					(32. / 9.) * vec_K2)));
		vec_K4 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_5th + dh * ((19372. / 6561.) * vec_K0 +
																				  (-25360. / 2187.) * vec_K1 +
																				   (64448. / 6561.) * vec_K2 +
																					 (-212. / 729.) * vec_K3)));
		vec_K5 = mat_ind_inverse * (vec_alpha + mat_res * (vec_current_5th + dh * ((9017. / 3168.) * vec_K0 +
																					 (-355. / 33.) * vec_K1 +
																				  (46732. / 5247.) * vec_K2 +
																					  (49. / 176.) * vec_K3 +
																				 (-5103. / 18656.) * vec_K4)));
		vec_current_4th = vec_current_5th + dh * ((35. / 384.) * vec_K0 +
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
		//RKDPupdate(vec_current_4th, vec_current_5th, mat_ind_inverse, mat_res, vec_alpha, dh, 0, 0);

		// 4・5次の局所誤差を計算
		for (int loop = 0; loop < n_loop; loop++) {
			*vec_error[loop] = abs(*vec_current_4th[loop] - *vec_current_5th[loop]) / dh;
		}

		// 4・5次の局所誤差ノルムがトレランス以上か評価、CSVに書き出せる値か判断
		bool output = true;
		double norm_error = vec_error.norm();
		if (norm_error > tolerance) output = false;
		cout << norm_error << ",\t";

		// 4・5次誤差評価から次の時間ステップの推測値を計算
		double delta = pow(tolerance / norm_error, 1. / 5.);
		//delta = pow(delta, 1. / 1.);	// 多分最適、無しでいい
		cout << delta << ",\t";

		// 次ステップの刻み幅を変化
		if (delta < 1.) {
			dh *= 0.2;	// 多分最適
		}
		else if (delta >= 2.) {
			dh *= 2.;
		}
		else {
			dh *= delta;
		}
		if (dh > dh_max) {
			dh = dh_max;
		}
		else if (dh < dh_min) {
			dh = dh_min;
		}

		// トレランス以下ならCSV書き出し、そうでないならque_max個前の結果に戻す
		if (output) {
			t += dh;
			static int count_interval = 0;
			if (count_interval++ % interval == 0) {
				csv_out_result << t << ",";
				for (int loop = 0; loop < n_loop; loop++) {
					//csv_out_result << *vec_current_4th[loop];
					csv_out_result << *vec_current_5th[loop];
					if (loop != n_loop - 1) {
						csv_out_result << ",";
					}
					else {
						csv_out_result << endl;
					}
				}
			}
			cout << "true" << endl;
		}
		else {
			t = que_t_old.front();
			vec_current_5th = que_current_5th_old.front();
			//dh = dh_init;
			cout << "False" << endl;
			//return 0;
		}

		// キューの最大要素数を超えていればpopする
		if (que_current_5th_old.size() > que_max) {
			que_t_old.pop();
			que_current_5th_old.pop();
		}

		// 計算終了か判断
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

int RKDPupdate(Matrix& current_4th, Matrix& current_5th,
			   Matrix& const inductance_inv, Matrix& const resistance, Matrix& const alpha, double dh,
			   int thread_num, int n_threads) {
	// RKDPの係数の入れ物用意
	Matrix vec_K[7];
	Matrix vec_temp[7];
	for (int i = 0; i < 7; i++) {
		new(&vec_K[i]) Matrix(n_loop);
		new(&vec_temp[i]) Matrix(n_loop);
	}

	// K0の計算
	//vec_K[0] = inductance_inv * (alpha + resistance * current_5th);
	omp_multi(vec_K[0], resistance, current_5th, thread_num, n_threads);	// vec_K[0] = resistance * current_5th
	omp_plus(vec_K[0], alpha, vec_K[0], thread_num, n_threads);				// vec_K[0] = alpha + vec_K[0]
	omp_multi(vec_K[0], inductance_inv, vec_K[0], thread_num, n_threads);	// vec_K[0] = inductance_inv * vec_K[0]

	// K1の計算
	//vec_K[1] = inductance_inv * (alpha + resistance * (current_5th + dh * (1. / 5.) * vec_K[0]));
	omp_multi(vec_temp[0], dh * (1. / 5.), vec_K[0], thread_num, n_threads);// vec_temp[0] = dh * (1. / 5.) * vec_K[0]
	omp_plus(vec_temp[0], current_5th, vec_temp[0], thread_num, n_threads);	// vec_temp[0] = current_5th + vec_temp[0]
	omp_multi(vec_K[1], resistance, vec_temp[0], thread_num, n_threads);	// vec_K[1] = resistance * vec_temp[0]
	omp_plus(vec_temp[0], alpha, vec_K[1], thread_num, n_threads);			// vec_temp[0] = alpha + vec_K[1]
	omp_multi(vec_K[1], inductance_inv, vec_temp[0], thread_num, n_threads);// vec_K[1] = inductance_inv * vec_temp[0]

	// K2の計算
	//vec_K[2] = inductance_inv * (alpha + resistance * (current_5th + dh * ((3. / 40.) * vec_K[0] +
	//																	 (9. / 40.) * vec_K[1])));
	omp_multi(vec_temp[0], dh * (3. / 40.), vec_K[0], thread_num, n_threads);	// vec_temp[0] = dh * (3. / 40.) * vec_K[0]
	omp_multi(vec_temp[1], dh * (9. / 40.), vec_K[1], thread_num, n_threads);	// vec_temp[1] = dh * (9. / 40.) * vec_K[1]
	omp_plus(vec_temp[1], vec_temp[0], vec_temp[1], thread_num, n_threads);		// vec_temp[1] = vec_temp[0] + vec_temp[1]
	omp_plus(vec_temp[1], current_5th, vec_temp[1], thread_num, n_threads);		// vec_temp[1] = current_5th + vec_temp[1]
	omp_multi(vec_K[2], resistance, vec_temp[1], thread_num, n_threads);		// vec_K[2] = resistance * vec_temp[1]
	omp_plus(vec_temp[1], alpha, vec_K[2], thread_num, n_threads);				// vec_temp[1] = alpha + vec_K[2]
	omp_multi(vec_K[2], inductance_inv, vec_temp[1], thread_num, n_threads);	// vec_K[2] = inductance_inv * vec_temp[1]

#pragma omp barrier
	//vec_K[3] = inductance_inv * (alpha + resistance * (current_5th + dh * ((44. / 45.) * vec_K[0] +
	//																	(-56. / 15.) * vec_K[1] +
	//																	  (32. / 9.) * vec_K[2])));
	omp_multi(vec_temp[0], dh * (44. / 45.), vec_K[0], thread_num, n_threads);	// vec_temp[0] = dh * (44. / 45.) * vec_K[0]
	omp_multi(vec_temp[1], dh * (-56. / 15.), vec_K[1], thread_num, n_threads);	// vec_temp[1] = dh * (-56. / 15.) * vec_K[1]
	omp_multi(vec_temp[2], dh * (32. / 9.), vec_K[2], thread_num, n_threads);	// vec_temp[2] = dh * (32. / 9.) * vec_K[2]
	omp_plus(vec_temp[2], vec_temp[0], vec_temp[2], thread_num, n_threads);		// vec_temp[2] = vec_temp[0] + vec_temp[2]
	omp_plus(vec_temp[2], vec_temp[1], vec_temp[2], thread_num, n_threads);		// vec_temp[2] = vec_temp[1] + vec_temp[2]
	omp_plus(vec_temp[2], current_5th, vec_temp[2], thread_num, n_threads);		// vec_temp[2] = current_5th + vec_temp[2]
	omp_multi(vec_K[3], resistance, vec_temp[2], thread_num, n_threads);		// vec_K[3] = resistance * vec_temp[2]
	omp_plus(vec_temp[2], alpha, vec_K[3], thread_num, n_threads);				// vec_temp[2] = alpha + vec_K[3]
	omp_multi(vec_K[3], inductance_inv, vec_temp[2], thread_num, n_threads);	// vec_K[3] = inductance_inv * vec_temp[2]

#pragma omp barrier
	//vec_K[4] = inductance_inv * (alpha + resistance * (current_5th + dh * ((19372. / 6561.) * vec_K[0] +
	//																	(-25360. / 2187.) * vec_K[1] +
	//																	 (64448. / 6561.) * vec_K[2] +
	//																	   (-212. / 729.) * vec_K[3])));
	omp_multi(vec_temp[0], dh * (19372. / 6561.), vec_K[0], thread_num, n_threads);	// vec_temp[0] = dh * (19372. / 6561.) * vec_K[0]
	omp_multi(vec_temp[1], dh * (-25360. / 2187.), vec_K[1], thread_num, n_threads);// vec_temp[1] = dh * (-25360. / 2187.) * vec_K[1]
	omp_multi(vec_temp[2], dh * (64448. / 6561.), vec_K[2], thread_num, n_threads);	// vec_temp[2] = dh * (64448. / 6561.) * vec_K[2]
	omp_multi(vec_temp[3], dh * (-212. / 729.), vec_K[3], thread_num, n_threads);	// vec_temp[3] = dh * (-212. / 729.) * vec_K[3]
	omp_plus(vec_temp[3], vec_temp[0], vec_temp[3], thread_num, n_threads);			// vec_temp[3] = vec_temp[0] + vec_temp[3]
	omp_plus(vec_temp[3], vec_temp[1], vec_temp[3], thread_num, n_threads);			// vec_temp[3] = vec_temp[1] + vec_temp[3]
	omp_plus(vec_temp[3], vec_temp[2], vec_temp[3], thread_num, n_threads);			// vec_temp[3] = vec_temp[2] + vec_temp[3]
	omp_plus(vec_temp[3], current_5th, vec_temp[3], thread_num, n_threads);			// vec_temp[3] = current_5th + vec_temp[3]
	omp_multi(vec_K[4], resistance, vec_temp[3], thread_num, n_threads);			// vec_K[4] = resistance * vec_temp[3]
	omp_plus(vec_temp[3], alpha, vec_K[4], thread_num, n_threads);					// vec_temp[3] = alpha + vec_K[4]
	omp_multi(vec_K[4], inductance_inv, vec_temp[3], thread_num, n_threads);		// vec_K[4] = inductance_inv * vec_temp[3]

#pragma omp barrier
	//vec_K[5] = inductance_inv * (alpha + resistance * (current_5th + dh * ((9017. / 3168.) * vec_K[0] +
	//																	   (-355. / 33.) * vec_K[1] +
	//																	(46732. / 5247.) * vec_K[2] +
	//																		(49. / 176.) * vec_K[3] +
	//																   (-5103. / 18656.) * vec_K[4])));
	omp_multi(vec_temp[0], dh * (9017. / 3168.), vec_K[0], thread_num, n_threads);	// vec_temp[0] = dh * (9017. / 3168.) * vec_K[0]
	omp_multi(vec_temp[1], dh * (-355. / 33.), vec_K[1], thread_num, n_threads);	// vec_temp[1] = dh * (-355. / 33.) * vec_K[1]
	omp_multi(vec_temp[2], dh * (46732. / 5247.), vec_K[2], thread_num, n_threads);	// vec_temp[2] = dh * (46732. / 5247.) * vec_K[2]
	omp_multi(vec_temp[3], dh * (49. / 176.), vec_K[3], thread_num, n_threads);		// vec_temp[3] = dh * (49. / 176.) * vec_K[3]
	omp_multi(vec_temp[4], dh * (-5103. / 18656.), vec_K[4], thread_num, n_threads);// vec_temp[4] = dh * (-5103. / 18656.) * vec_K[4]
	omp_plus(vec_temp[4], vec_temp[0], vec_temp[4], thread_num, n_threads);			// vec_temp[4] = vec_temp[0] + vec_temp[4]
	omp_plus(vec_temp[4], vec_temp[1], vec_temp[4], thread_num, n_threads);			// vec_temp[4] = vec_temp[1] + vec_temp[4]
	omp_plus(vec_temp[4], vec_temp[2], vec_temp[4], thread_num, n_threads);			// vec_temp[4] = vec_temp[2] + vec_temp[4]
	omp_plus(vec_temp[4], vec_temp[3], vec_temp[4], thread_num, n_threads);			// vec_temp[4] = vec_temp[3] + vec_temp[4]
	omp_plus(vec_temp[4], current_5th, vec_temp[4], thread_num, n_threads);			// vec_temp[4] = current_5th + vec_temp[4]
	omp_multi(vec_K[5], resistance, vec_temp[4], thread_num, n_threads);			// vec_K[5] = resistance * vec_temp[4]
	omp_plus(vec_temp[4], alpha, vec_K[5], thread_num, n_threads);					// vec_temp[4] = alpha + vec_K[5]
	omp_multi(vec_K[5], inductance_inv, vec_temp[4], thread_num, n_threads);		// vec_K[5] = inductance_inv * vec_temp[4]

#pragma omp barrier
	//current_4th = current_5th + dh * ((35. / 384.) * vec_K0 +
	//										//(0.) * vec_K1 +
	//								(500. / 1113.) * vec_K2 +
	//								 (125. / 192.) * vec_K3 +
	//							  (-2187. / 6784.) * vec_K4 +
	//								   (11. / 84.) * vec_K5);
	omp_multi(vec_temp[0], dh * (35. / 384.), vec_K[0], thread_num, n_threads);		// vec_temp[0] = dh * (35. / 384.) * vec_K[0]
	omp_multi(vec_temp[1], dh * (0.), vec_K[1], thread_num, n_threads);				// vec_temp[1] = dh * (0.) * vec_K[1]
	omp_multi(vec_temp[2], dh * (500. / 1113.), vec_K[2], thread_num, n_threads);	// vec_temp[2] = dh * (500. / 1113.) * vec_K[2]
	omp_multi(vec_temp[3], dh * (125. / 192.), vec_K[3], thread_num, n_threads);	// vec_temp[3] = dh * (125. / 192.) * vec_K[3]
	omp_multi(vec_temp[4], dh * (-2187. / 6784.), vec_K[4], thread_num, n_threads);	// vec_temp[4] = dh * (-2187. / 6784.) * vec_K[4]
	omp_multi(vec_temp[5], dh * (11. / 84.), vec_K[5], thread_num, n_threads);		// vec_temp[5] = dh * (11. / 84.) * vec_K[5]
	omp_plus(vec_temp[5], vec_temp[0], vec_temp[5], thread_num, n_threads);			// vec_temp[5] = vec_temp[0] + vec_temp[5]
	omp_plus(vec_temp[5], vec_temp[1], vec_temp[5], thread_num, n_threads);			// vec_temp[5] = vec_temp[1] + vec_temp[5]
	omp_plus(vec_temp[5], vec_temp[2], vec_temp[5], thread_num, n_threads);			// vec_temp[5] = vec_temp[2] + vec_temp[5]
	omp_plus(vec_temp[5], vec_temp[3], vec_temp[5], thread_num, n_threads);			// vec_temp[5] = vec_temp[3] + vec_temp[5]
	omp_plus(vec_temp[5], vec_temp[4], vec_temp[5], thread_num, n_threads);			// vec_temp[5] = vec_temp[4] + vec_temp[5]
	omp_plus(current_4th, current_5th, vec_temp[5], thread_num, n_threads);			// current_4th = current_5th + vec_temp[5]

#pragma omp barrier
	//vec_K6 = inductance_inv * (alpha + resistance * current_4th);
	omp_multi(vec_K[6], resistance, current_4th, thread_num, n_threads);			// vec_K[6] = resistance * current_4th
	omp_plus(vec_temp[6], alpha, vec_K[6], thread_num, n_threads);					// vec_temp[6] = alpha + vec_K[6]
	omp_multi(vec_K[6], inductance_inv, vec_temp[6], thread_num, n_threads);		// vec_K[6] = inductance_inv * vec_temp[6]

#pragma omp barrier
	//current_5th = current_5th + dh * ((5179. / 57600.) * vec_K0 +
	//	//(0.) * vec_K1 +
	//	(7571. / 16695.) * vec_K2 +
	//	(393. / 640.) * vec_K3 +
	//	(-92097. / 339200.) * vec_K4 +
	//	(187. / 2100.) * vec_K5 +
	//	(1. / 40.) * vec_K6);
	omp_multi(vec_temp[0], dh * (5179. / 57600.), vec_K[0], thread_num, n_threads);	// vec_temp[0] = dh * (5179. / 57600.) * vec_K[0]
	omp_multi(vec_temp[1], dh * (0.), vec_K[1], thread_num, n_threads);				// vec_temp[1] = dh * (0.) * vec_K[1]
	omp_multi(vec_temp[2], dh * (7571. / 16695.), vec_K[2], thread_num, n_threads);	// vec_temp[2] = dh * (7571. / 16695.) * vec_K[2]
	omp_multi(vec_temp[3], dh * (393. / 640.), vec_K[3], thread_num, n_threads);	// vec_temp[3] = dh * (393. / 640.) * vec_K[3]
	omp_multi(vec_temp[4], dh * (-92097. / 339200.), vec_K[4], thread_num, n_threads);// vec_temp[4] = dh * (-92097. / 339200.) * vec_K[4]
	omp_multi(vec_temp[5], dh * (187. / 2100.), vec_K[5], thread_num, n_threads);	// vec_temp[5] = dh * (187. / 2100.) * vec_K[5]
	omp_multi(vec_temp[6], dh * (1. / 40.), vec_K[6], thread_num, n_threads);		// vec_temp[6] = dh * (1. / 40.) * vec_K[6]
	omp_plus(vec_temp[6], vec_temp[0], vec_temp[6], thread_num, n_threads);			// vec_temp[6] = vec_temp[0] + vec_temp[6]
	omp_plus(vec_temp[6], vec_temp[1], vec_temp[6], thread_num, n_threads);			// vec_temp[6] = vec_temp[1] + vec_temp[6]
	omp_plus(vec_temp[6], vec_temp[2], vec_temp[6], thread_num, n_threads);			// vec_temp[6] = vec_temp[2] + vec_temp[6]
	omp_plus(vec_temp[6], vec_temp[3], vec_temp[6], thread_num, n_threads);			// vec_temp[6] = vec_temp[3] + vec_temp[6]
	omp_plus(vec_temp[6], vec_temp[4], vec_temp[6], thread_num, n_threads);			// vec_temp[6] = vec_temp[4] + vec_temp[6]
	omp_plus(vec_temp[6], vec_temp[5], vec_temp[6], thread_num, n_threads);			// vec_temp[6] = vec_temp[5] + vec_temp[6]
	omp_plus(current_5th, current_5th, vec_temp[6], thread_num, n_threads);			// current_4th = current_5th + vec_temp[6]

#pragma omp barrier

	return 0;
}
