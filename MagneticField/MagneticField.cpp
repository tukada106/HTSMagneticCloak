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
#define n_layer 2	// 2layers' condition
#define n_turn 25	// 2layers' condition
#define n_ring (n_turn * 2)
#define n_loop (n_layer - 1) * n_ring

#define step 500

using namespace std;

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);

int main() {
	// 磁場計算する点を指定(x, y, z)
	Matrix p_evaluate(3);
	*p_evaluate[0] = 0.;
	*p_evaluate[1] = 0.;
	*p_evaluate[2] = 0.05;

	// 時間計測開始
	clock_t startTime, processTime;
	startTime = clock();

	// ファイル読み取りオープン
	ostringstream oss_file_name;
	oss_file_name << "../Current/result_current.csv";
	ifstream ifs_csv_in(oss_file_name.str());

	// 電流CSVの行数を数える
	int rows_current = 0;
	string str_csv_in_line;
	while (getline(ifs_csv_in, str_csv_in_line)) rows_current++;
	ifs_csv_in.clear();
	ifs_csv_in.seekg(0, ios_base::beg);

	// 電流CSVをMatrixへ読み込み
	Matrix current(rows_current, n_loop + 1);
	string str_csv_in_comma;
	for (int row = 0; row < current.row_size(); row++) {
		if (getline(ifs_csv_in, str_csv_in_line)) {
			istringstream iss_csv_in_line(str_csv_in_line);
			for (int column = 0; column < current.column_size(); column++) {
				if (getline(iss_csv_in_line, str_csv_in_comma, ',')) {
					current[row][column] = stod(str_csv_in_comma);
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

	// Biot-Savartの法則
	Matrix MagField(current.row_size(), 3);
	Matrix temp_current(current.row_size());
	for (int column = 0; column < current.column_size(); column++) {
		// columnループ目の電流を全時間分取得
		for (int row = 0; row < current.row_size(); row++) {
			*temp_current[row] = current[row][column];
		}

		// ターン数の厚さ方向の判定
		if (column < n_ring) {

		}
		else if(column <2*n_ring){}
		while()

		// 

	}
	
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

Matrix BiotSavart(IN Matrix point, IN Matrix pos,IN Matrix tang, IN Matrix current) {
	// point	:磁場計算を行う位置ベクトル（x, y, z）
	// pos		:電流経路の位置ベクトル（PosVecを使う）
	// tang		:電流経路の接線ベクトル（TangLinVecを使う）
	// current	:1ループ分、全時間の電流
	// すべての引数はMatrixクラスを使用すること

	if (point.row_size() != 3 || point.column_size() != 1) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : BiotSavart" << endl;
		exit(1);
	}
	if (pos.row_size() != tang.row_size()) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : BiotSavart" << endl;
		exit(1);
	}
	if (pos.column_size() != 3 || tang.column_size() != 3) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : BiotSavart" << endl;
		exit(1);
	}
	if (current.column_size() != 1) {
		cerr << "Matrix size mismatch!" << endl;
		cerr << "ERR : BiotSavart" << endl;
		exit(1);
	}

	// 線素から見たpointの相対位置ベクトルを計算
	Matrix pos_rel(pos.row_size(), 4);
	for (int row = 0; row < pos.row_size(); row++) {
		pos_rel[row][0] = *point[0] - pos[row][0];
		pos_rel[row][1] = *point[1] - pos[row][1];
		pos_rel[row][2] = *point[2] - pos[row][2];
		pos_rel[row][3] = sqrt(pow(pos[row][0], 2.) + pow(pos[row][1], 2.) + pow(pos[row][2], 2.));
	}

	// Biot-Savartの法則でpointでの磁束密度ベクトルを計算
	Matrix MagField(current.row_size(), 3);
	for (int row_current = 0; row_current < current.row_size(); row_current++) {
		for (int row_tang = 0; row_tang < tang.row_size(); row_tang++) {
			MagField[row_current][0] += (tang[row_tang][1] * pos_rel[row_tang][2] - tang[row_tang][2] * pos_rel[row_tang][1]) / pow(pos_rel[row_tang][3], 3.);
			MagField[row_current][1] += (tang[row_tang][2] * pos_rel[row_tang][0] - tang[row_tang][0] * pos_rel[row_tang][2]) / pow(pos_rel[row_tang][3], 3.);
			MagField[row_current][2] += (tang[row_tang][0] * pos_rel[row_tang][1] - tang[row_tang][1] * pos_rel[row_tang][0]) / pow(pos_rel[row_tang][3], 3.);
		}
		MagField[row_current][0]*=current[]
	}
}
