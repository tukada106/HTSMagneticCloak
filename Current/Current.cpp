#include <iostream>
#include <fstream>
#include <sstream>

#include "../Inductance/matrix.h"

#define IN
#define OUT
#define Pi 3.14159265
#define mu0 (4 * Pi * 1e-7)

#define r_shield 0.015
#define w_tape 0.012
#define t_tape 0.1e-3
#define n_layer 4
#define n_turn 10
#define ring_st 0
#define ring_sp (n_turn * 2)

#define B_apply 8.48e-3
#define t_sweep 0.5e-3

using namespace std;

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	// ファイル読み取りオープン
	ifstream ifs_csv_in[n_layer - 1];
	for (int i = 0; i < n_layer - 1; i++) {
		ostringstream oss_file_name;
		oss_file_name << "../Inductance/";
		oss_file_name << "inductance_layer_base_" << i << ".csv";
		new(&ifs_csv_in[i]) ifstream(oss_file_name.str());
	}

	Matrix ind[n_layer - 1];
	for (int i = 0; i < n_layer - 1; i++) {
		new(&ind[i]) Matrix(ring_sp, n_layer - 1);
	}

	for (int layer_base = 0; layer_base < n_layer - 1; layer_base++) {
		string str_csv_in_line;
		string str_csv_in_comma;
		if (getline(ifs_csv_in[layer_base], str_csv_in_line)) {
			istringstream iss_csv_in_line(str_csv_in_line);
			for (int layer = 0; layer < n_layer - 1; layer++) {
				if (getline(iss_csv_in_line, str_csv_in_comma, ',')) {
					ind[layer_base][]
				}
				else {
					cerr << "The specified layer is missing!" << endl;
					cerr << "getline(ifs_csv_in[layer], str_csv_in_line)" << endl;
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

	// ファイル読み取りクローズ
	for (int i = 0; i < n_layer - 1; i++) {
		ifs_csv_in[i].close();
	}

	processTime = clock() - startTime;
	cout << static_cast<double>(processTime) / 1000 << " [s]" << endl;
}
