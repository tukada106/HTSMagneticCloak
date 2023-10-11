#include <iostream>
#include <fstream>
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
#define n_turn 10
#define ring_st 0
#define ring_sp (n_turn * 2)

#define step 500

using namespace std;

Matrix PosVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix PosVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix PosVec_Rad(IN Matrix init, IN Matrix end, IN long n);
Matrix TangLinVec_For(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_Rev(IN double init, IN double end, IN long n, IN double rad, IN double wid, IN double thick);
Matrix TangLinVec_Rad(IN Matrix init, IN Matrix end, IN long n);

int main() {
	clock_t startTime, processTime;
	startTime = clock();

	ofstream ofs_out("output.csv");

	Matrix ind0(ring_sp, n_layer);

	Matrix pos0_base = PosVec_For(0, Pi, step, r_shield, w_tape, 0);
	Matrix tang0_base = TangLinVec_For(0, Pi, step, r_shield, w_tape, 0);
	for (int layer = 0; layer < n_layer; layer++) {
		for (int ring = 0; ring < ring_sp; ring++) {
			double inductance = 0;
			double dist = 0;
			double dotPro = 0;
			double thickness;
			if (layer == 0 && ring == 0) {
				thickness = t_tape;
			}
			else {
				thickness = 0.;
			}

			Matrix pos(step, 3);
			Matrix tang(step, 3);
			if (layer % 2 == 0) {
				pos = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
				tang = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
			}
			else {
				pos = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
				tang = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + layer * t_tape, w_tape, thickness);
			}
			
			for (int i = 0; i < step; i++) {
				for (int j = 0; j < step; j++) {
					dist = sqrt(pow(pos0_base[i][0] - pos[j][0], 2.) +
								pow(pos0_base[i][1] - pos[j][1], 2.) +
								pow(pos0_base[i][2] - pos[j][2], 2.));
					dotPro = tang0_base[i][0] * tang[j][0] +
							 tang0_base[i][1] * tang[j][1] +
							 tang0_base[i][2] * tang[j][2];
					inductance += dotPro / dist;
				}
			}

			ind0[ring][layer] = mu0 / (4 * Pi) * inductance;
			cout << layer << ", " << ring << endl;
		}
	}

	//Matrix pos0_center = PosVec_For(0, Pi, step, r_shield, w_tape, 0);
	//Matrix pos1_center = PosVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape, 0);
	//Matrix pos_radial01_center = PosVec_Rad(pos0_center.extract_row(step).transposed(), pos1_center.extract_row(1).transposed(), step / 100);
	//Matrix pos_radial10_center = PosVec_Rad(pos1_center.extract_row(step).transposed(), pos0_center.extract_row(1).transposed(), step / 100);
	//Matrix tang0_center = TangLinVec_For(0, Pi, step, r_shield, w_tape, 0);
	//Matrix tang1_center = TangLinVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape, 0);
	//Matrix tang_radial01_center = TangLinVec_Rad(pos0_center.extract_row(step).transposed(), pos1_center.extract_row(1).transposed(), step / 100);
	//Matrix tang_radial10_center = TangLinVec_Rad(pos1_center.extract_row(step).transposed(), pos0_center.extract_row(1).transposed(), step / 100);
	//
	//Matrix pos0_edge = PosVec_For(0, Pi, step, r_shield, w_tape, t_tape);
	//Matrix pos1_edge = PosVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape, t_tape);
	//Matrix pos_radial01_edge = PosVec_Rad(pos0_edge.extract_row(step).transposed(), pos1_edge.extract_row(1).transposed(), step / 100);
	//Matrix pos_radial10_edge = PosVec_Rad(pos1_edge.extract_row(step).transposed(), pos0_edge.extract_row(1).transposed(), step / 100);
	//Matrix tang0_edge = TangLinVec_For(0, Pi, step, r_shield, w_tape, t_tape);
	//Matrix tang1_edge = TangLinVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape, t_tape);
	//Matrix tang_radial01_edge = TangLinVec_Rad(pos0_edge.extract_row(step).transposed(), pos1_edge.extract_row(1).transposed(), step / 100);
	//Matrix tang_radial10_edge = TangLinVec_Rad(pos1_edge.extract_row(step).transposed(), pos0_edge.extract_row(1).transposed(), step / 100);

	//int i, j;
	//// 自己インダクタンス計算
	//{
	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos0_center[i][0] - pos0_edge[j][0], 2.) +
	//				pow(pos0_center[i][1] - pos0_edge[j][1], 2.) +
	//				pow(pos0_center[i][2] - pos0_edge[j][2], 2.));
	//			dotPro = tang0_center[i][0] * tang0_edge[j][0] +
	//				tang0_center[i][1] * tang0_edge[j][1] +
	//				tang0_center[i][2] * tang0_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "1 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos0_center[i][0] - pos1_edge[j][0], 2.) +
	//				pow(pos0_center[i][1] - pos1_edge[j][1], 2.) +
	//				pow(pos0_center[i][2] - pos1_edge[j][2], 2.));
	//			dotPro = tang0_center[i][0] * tang1_edge[j][0] +
	//				tang0_center[i][1] * tang1_edge[j][1] +
	//				tang0_center[i][2] * tang1_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "2 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos0_center[i][0] - pos_radial01_edge[j][0], 2.) +
	//				pow(pos0_center[i][1] - pos_radial01_edge[j][1], 2.) +
	//				pow(pos0_center[i][2] - pos_radial01_edge[j][2], 2.));
	//			dotPro = tang0_center[i][0] * tang_radial01_edge[0][0] +
	//				tang0_center[i][1] * tang_radial01_edge[1][0] +
	//				tang0_center[i][2] * tang_radial01_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "3 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos0_center[i][0] - pos_radial10_edge[j][0], 2.) +
	//				pow(pos0_center[i][1] - pos_radial10_edge[j][1], 2.) +
	//				pow(pos0_center[i][2] - pos_radial10_edge[j][2], 2.));
	//			dotPro = tang0_center[i][0] * tang_radial10_edge[0][0] +
	//				tang0_center[i][1] * tang_radial10_edge[1][0] +
	//				tang0_center[i][2] * tang_radial10_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "4 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos1_center[i][0] - pos0_edge[j][0], 2.) +
	//				pow(pos1_center[i][1] - pos0_edge[j][1], 2.) +
	//				pow(pos1_center[i][2] - pos0_edge[j][2], 2.));
	//			dotPro = tang1_center[i][0] * tang0_edge[j][0] +
	//				tang1_center[i][1] * tang0_edge[j][1] +
	//				tang1_center[i][2] * tang0_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "5 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos1_center[i][0] - pos1_edge[j][0], 2.) +
	//				pow(pos1_center[i][1] - pos1_edge[j][1], 2.) +
	//				pow(pos1_center[i][2] - pos1_edge[j][2], 2.));
	//			dotPro = tang1_center[i][0] * tang1_edge[j][0] +
	//				tang1_center[i][1] * tang1_edge[j][1] +
	//				tang1_center[i][2] * tang1_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "6 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos1_center[i][0] - pos_radial01_edge[j][0], 2.) +
	//				pow(pos1_center[i][1] - pos_radial01_edge[j][1], 2.) +
	//				pow(pos1_center[i][2] - pos_radial01_edge[j][2], 2.));
	//			dotPro = tang1_center[i][0] * tang_radial01_edge[0][0] +
	//				tang1_center[i][1] * tang_radial01_edge[1][0] +
	//				tang1_center[i][2] * tang_radial01_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "7 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos1_center[i][0] - pos_radial10_edge[j][0], 2.) +
	//				pow(pos1_center[i][1] - pos_radial10_edge[j][1], 2.) +
	//				pow(pos1_center[i][2] - pos_radial10_edge[j][2], 2.));
	//			dotPro = tang1_center[i][0] * tang_radial10_edge[0][0] +
	//				tang1_center[i][1] * tang_radial10_edge[1][0] +
	//				tang1_center[i][2] * tang_radial10_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "8 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial01_center[i][0] - pos0_edge[j][0], 2.) +
	//				pow(pos_radial01_center[i][1] - pos0_edge[j][1], 2.) +
	//				pow(pos_radial01_center[i][2] - pos0_edge[j][2], 2.));
	//			dotPro = tang_radial01_center[0][0] * tang0_edge[j][0] +
	//				tang_radial01_center[1][0] * tang0_edge[j][1] +
	//				tang_radial01_center[2][0] * tang0_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "9 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial01_center[i][0] - pos1_edge[j][0], 2.) +
	//				pow(pos_radial01_center[i][1] - pos1_edge[j][1], 2.) +
	//				pow(pos_radial01_center[i][2] - pos1_edge[j][2], 2.));
	//			dotPro = tang_radial01_center[0][0] * tang1_edge[j][0] +
	//				tang_radial01_center[1][0] * tang1_edge[j][1] +
	//				tang_radial01_center[2][0] * tang1_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "10 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial01_center[i][0] - pos_radial01_edge[j][0], 2.) +
	//				pow(pos_radial01_center[i][1] - pos_radial01_edge[j][1], 2.) +
	//				pow(pos_radial01_center[i][2] - pos_radial01_edge[j][2], 2.));
	//			dotPro = tang_radial01_center[0][0] * tang_radial01_edge[0][0] +
	//				tang_radial01_center[1][0] * tang_radial01_edge[1][0] +
	//				tang_radial01_center[2][0] * tang_radial01_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "11 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial01_center[i][0] - pos_radial10_edge[j][0], 2.) +
	//				pow(pos_radial01_center[i][1] - pos_radial10_edge[j][1], 2.) +
	//				pow(pos_radial01_center[i][2] - pos_radial10_edge[j][2], 2.));
	//			dotPro = tang_radial01_center[0][0] * tang_radial10_edge[0][0] +
	//				tang_radial01_center[1][0] * tang_radial10_edge[1][0] +
	//				tang_radial01_center[2][0] * tang_radial10_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "12 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial10_center[i][0] - pos0_edge[j][0], 2.) +
	//				pow(pos_radial10_center[i][1] - pos0_edge[j][1], 2.) +
	//				pow(pos_radial10_center[i][2] - pos0_edge[j][2], 2.));
	//			dotPro = tang_radial10_center[0][0] * tang0_edge[j][0] +
	//				tang_radial10_center[1][0] * tang0_edge[j][1] +
	//				tang_radial10_center[2][0] * tang0_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "13 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial10_center[i][0] - pos1_edge[j][0], 2.) +
	//				pow(pos_radial10_center[i][1] - pos1_edge[j][1], 2.) +
	//				pow(pos_radial10_center[i][2] - pos1_edge[j][2], 2.));
	//			dotPro = tang_radial10_center[0][0] * tang1_edge[j][0] +
	//				tang_radial10_center[1][0] * tang1_edge[j][1] +
	//				tang_radial10_center[2][0] * tang1_edge[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "14 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial10_center[i][0] - pos_radial01_edge[j][0], 2.) +
	//				pow(pos_radial10_center[i][1] - pos_radial01_edge[j][1], 2.) +
	//				pow(pos_radial10_center[i][2] - pos_radial01_edge[j][2], 2.));
	//			dotPro = tang_radial10_center[0][0] * tang_radial01_edge[0][0] +
	//				tang_radial10_center[1][0] * tang_radial01_edge[1][0] +
	//				tang_radial10_center[2][0] * tang_radial01_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "15 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial10_center[i][0] - pos_radial10_edge[j][0], 2.) +
	//				pow(pos_radial10_center[i][1] - pos_radial10_edge[j][1], 2.) +
	//				pow(pos_radial10_center[i][2] - pos_radial10_edge[j][2], 2.));
	//			dotPro = tang_radial10_center[0][0] * tang_radial10_edge[0][0] +
	//				tang_radial10_center[1][0] * tang_radial10_edge[1][0] +
	//				tang_radial10_center[2][0] * tang_radial10_edge[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "16 : " << i + 1 << "/" << step << endl;
	//	}
	//}

	//cout << 0 << "/" << ring_sp << " : " << mu0 / (4 * Pi) * inductance << " [H]" << endl;
	//ofs_out << 0 << "," << mu0 / (4 * Pi) * inductance << "\n";

	//// 相互インダクタンス計算
	//for (int ring = ring_st; ring <= ring_sp; ring++) {
	//	Matrix pos0_base = PosVec_For(0, Pi, step, r_shield, w_tape, 0);
	//	Matrix pos1_base = PosVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape, 0);
	//	Matrix pos_radial01_base = PosVec_Rad(pos0_base.extract_row(step).transposed(), pos1_base.extract_row(1).transposed(), step / 100);
	//	Matrix pos_radial10_base = PosVec_Rad(pos1_base.extract_row(step).transposed(), pos0_base.extract_row(1).transposed(), step / 100);
	//	Matrix tang0_base = TangLinVec_For(0, Pi, step, r_shield, w_tape, 0);
	//	Matrix tang1_base = TangLinVec_Rev(-Pi, 0, step, r_shield + t_tape, w_tape, 0);
	//	Matrix tang_radial01_base = TangLinVec_Rad(pos0_base.extract_row(step).transposed(), pos1_base.extract_row(1).transposed(), step / 100);
	//	Matrix tang_radial10_base = TangLinVec_Rad(pos1_base.extract_row(step).transposed(), pos0_base.extract_row(1).transposed(), step / 100);
	//	
	//	Matrix pos0 = PosVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield, w_tape, 0);
	//	Matrix pos1 = PosVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + t_tape, w_tape, 0);
	//	Matrix pos_radial01 = PosVec_Rad(pos0.extract_row(step).transposed(), pos1.extract_row(1).transposed(), step / 100);
	//	Matrix pos_radial10 = PosVec_Rad(pos1.extract_row(step).transposed(), pos0.extract_row(1).transposed(), step / 100);
	//	Matrix tang0 = TangLinVec_For(ring * Pi, (ring + 1) * Pi, step, r_shield, w_tape, 0);
	//	Matrix tang1 = TangLinVec_Rev(-(ring + 1) * Pi, -ring * Pi, step, r_shield + t_tape, w_tape, 0);
	//	Matrix tang_radial01 = TangLinVec_Rad(pos0.extract_row(step).transposed(), pos1.extract_row(1).transposed(), step);
	//	Matrix tang_radial10 = TangLinVec_Rad(pos1.extract_row(step).transposed(), pos0.extract_row(1).transposed(), step);

	//	inductance = 0;
	//	dist = 0;
	//	dotPro = 0;

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos0_base[i][0] - pos0[j][0], 2.) +
	//						pow(pos0_base[i][1] - pos0[j][1], 2.) +
	//						pow(pos0_base[i][2] - pos0[j][2], 2.));
	//			dotPro = tang0_base[i][0] * tang0[j][0] +
	//					 tang0_base[i][1] * tang0[j][1] +
	//					 tang0_base[i][2] * tang0[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "1 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos0_base[i][0] - pos1[j][0], 2.) +
	//						pow(pos0_base[i][1] - pos1[j][1], 2.) +
	//						pow(pos0_base[i][2] - pos1[j][2], 2.));
	//			dotPro = tang0_base[i][0] * tang1[j][0] +
	//					 tang0_base[i][1] * tang1[j][1] +
	//					 tang0_base[i][2] * tang1[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "2 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos0_base[i][0] - pos_radial01[j][0], 2.) +
	//						pow(pos0_base[i][1] - pos_radial01[j][1], 2.) +
	//						pow(pos0_base[i][2] - pos_radial01[j][2], 2.));
	//			dotPro = tang0_base[i][0] * tang_radial01[0][0] +
	//					 tang0_base[i][1] * tang_radial01[1][0] +
	//					 tang0_base[i][2] * tang_radial01[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "3 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos0_base[i][0] - pos_radial10[j][0], 2.) +
	//						pow(pos0_base[i][1] - pos_radial10[j][1], 2.) +
	//						pow(pos0_base[i][2] - pos_radial10[j][2], 2.));
	//			dotPro = tang0_base[i][0] * tang_radial10[0][0] +
	//					 tang0_base[i][1] * tang_radial10[1][0] +
	//					 tang0_base[i][2] * tang_radial10[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "4 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos1_base[i][0] - pos0[j][0], 2.) +
	//						pow(pos1_base[i][1] - pos0[j][1], 2.) +
	//						pow(pos1_base[i][2] - pos0[j][2], 2.));
	//			dotPro = tang1_base[i][0] * tang0[j][0] +
	//					 tang1_base[i][1] * tang0[j][1] +
	//					 tang1_base[i][2] * tang0[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "5 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos1_base[i][0] - pos1[j][0], 2.) +
	//						pow(pos1_base[i][1] - pos1[j][1], 2.) +
	//						pow(pos1_base[i][2] - pos1[j][2], 2.));
	//			dotPro = tang1_base[i][0] * tang1[j][0] +
	//					 tang1_base[i][1] * tang1[j][1] +
	//					 tang1_base[i][2] * tang1[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "6 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos1_base[i][0] - pos_radial01[j][0], 2.) +
	//						pow(pos1_base[i][1] - pos_radial01[j][1], 2.) +
	//						pow(pos1_base[i][2] - pos_radial01[j][2], 2.));
	//			dotPro = tang1_base[i][0] * tang_radial01[0][0] +
	//					 tang1_base[i][1] * tang_radial01[1][0] +
	//					 tang1_base[i][2] * tang_radial01[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "7 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos1_base[i][0] - pos_radial10[j][0], 2.) +
	//						pow(pos1_base[i][1] - pos_radial10[j][1], 2.) +
	//						pow(pos1_base[i][2] - pos_radial10[j][2], 2.));
	//			dotPro = tang1_base[i][0] * tang_radial10[0][0] +
	//					 tang1_base[i][1] * tang_radial10[1][0] +
	//					 tang1_base[i][2] * tang_radial10[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "8 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial01_base[i][0] - pos0[j][0], 2.) +
	//						pow(pos_radial01_base[i][1] - pos0[j][1], 2.) +
	//						pow(pos_radial01_base[i][2] - pos0[j][2], 2.));
	//			dotPro = tang_radial01_base[0][0] * tang0[j][0] +
	//					 tang_radial01_base[1][0] * tang0[j][1] +
	//					 tang_radial01_base[2][0] * tang0[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "9 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial01_base[i][0] - pos1[j][0], 2.) +
	//						pow(pos_radial01_base[i][1] - pos1[j][1], 2.) +
	//						pow(pos_radial01_base[i][2] - pos1[j][2], 2.));
	//			dotPro = tang_radial01_base[0][0] * tang1[j][0] +
	//					 tang_radial01_base[1][0] * tang1[j][1] +
	//					 tang_radial01_base[2][0] * tang1[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "10 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial01_base[i][0] - pos_radial01[j][0], 2.) +
	//						pow(pos_radial01_base[i][1] - pos_radial01[j][1], 2.) +
	//						pow(pos_radial01_base[i][2] - pos_radial01[j][2], 2.));
	//			dotPro = tang_radial01_base[0][0] * tang_radial01[0][0] +
	//					 tang_radial01_base[1][0] * tang_radial01[1][0] +
	//					 tang_radial01_base[2][0] * tang_radial01[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "11 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial01_base[i][0] - pos_radial10[j][0], 2.) +
	//						pow(pos_radial01_base[i][1] - pos_radial10[j][1], 2.) +
	//						pow(pos_radial01_base[i][2] - pos_radial10[j][2], 2.));
	//			dotPro = tang_radial01_base[0][0] * tang_radial10[0][0] +
	//					 tang_radial01_base[1][0] * tang_radial10[1][0] +
	//					 tang_radial01_base[2][0] * tang_radial10[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "12 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial10_base[i][0] - pos0[j][0], 2.) +
	//						pow(pos_radial10_base[i][1] - pos0[j][1], 2.) +
	//						pow(pos_radial10_base[i][2] - pos0[j][2], 2.));
	//			dotPro = tang_radial10_base[0][0] * tang0[j][0] +
	//					 tang_radial10_base[1][0] * tang0[j][1] +
	//					 tang_radial10_base[2][0] * tang0[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "13 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step; j++) {
	//			dist = sqrt(pow(pos_radial10_base[i][0] - pos1[j][0], 2.) +
	//						pow(pos_radial10_base[i][1] - pos1[j][1], 2.) +
	//						pow(pos_radial10_base[i][2] - pos1[j][2], 2.));
	//			dotPro = tang_radial10_base[0][0] * tang1[j][0] +
	//					 tang_radial10_base[1][0] * tang1[j][1] +
	//					 tang_radial10_base[2][0] * tang1[j][2];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "14 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial10_base[i][0] - pos_radial01[j][0], 2.) +
	//						pow(pos_radial10_base[i][1] - pos_radial01[j][1], 2.) +
	//						pow(pos_radial10_base[i][2] - pos_radial01[j][2], 2.));
	//			dotPro = tang_radial10_base[0][0] * tang_radial01[0][0] +
	//					 tang_radial10_base[1][0] * tang_radial01[1][0] +
	//					 tang_radial10_base[2][0] * tang_radial01[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "15 : " << i + 1 << "/" << step << endl;
	//	}

	//	for (i = 0; i < step / 100; i++) {
	//		for (j = 0; j < step / 100; j++) {
	//			dist = sqrt(pow(pos_radial10_base[i][0] - pos_radial10[j][0], 2.) +
	//						pow(pos_radial10_base[i][1] - pos_radial10[j][1], 2.) +
	//						pow(pos_radial10_base[i][2] - pos_radial10[j][2], 2.));
	//			dotPro = tang_radial10_base[0][0] * tang_radial10[0][0] +
	//					 tang_radial10_base[1][0] * tang_radial10[1][0] +
	//					 tang_radial10_base[2][0] * tang_radial10[2][0];
	//			inductance += dotPro / dist;
	//		}
	//		//cout << "16 : " << i + 1 << "/" << step << endl;
	//	}


	//	cout << ring << "/" << ring_sp << " : " << mu0 / (4 * Pi) * inductance << " [H]" << endl;
	//	ofs_out << ring << "," << mu0 / (4 * Pi) * inductance << "\n";
	//}

	/*
#pragma omp parallel
{
	int i, j;
	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (i = 0; i < step; i++) {
		for (j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_For[i][0] - pos2_For[j][0], 2.) +
						pow(pos1_For[i][1] - pos2_For[j][1], 2.) +
						pow(pos1_For[i][2] - pos2_For[j][2], 2.));
			dotPro = tang1_For[i][0] * tang2_For[j][0] +
					 tang1_For[i][1] * tang2_For[j][1] +
					 tang1_For[i][2] * tang2_For[j][2];
			inductance += dotPro / dist;
		}
		cout << "1 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;

	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_For[i][0] - pos2_Rev[j][0], 2.) +
						pow(pos1_For[i][1] - pos2_Rev[j][1], 2.) +
						pow(pos1_For[i][2] - pos2_Rev[j][2], 2.));
			dotPro = tang1_For[i][0] * tang2_Rev[j][0] +
					 tang1_For[i][1] * tang2_Rev[j][1] +
					 tang1_For[i][2] * tang2_Rev[j][2];
			inductance += dotPro / dist;
		}
		cout << "2 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;

	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_Rev[i][0] - pos2_For[j][0], 2.) +
						pow(pos1_Rev[i][1] - pos2_For[j][1], 2.) +
						pow(pos1_Rev[i][2] - pos2_For[j][2], 2.));
			dotPro = tang1_Rev[i][0] * tang2_For[j][0] +
					 tang1_Rev[i][1] * tang2_For[j][1] +
					 tang1_Rev[i][2] * tang2_For[j][2];
			inductance += dotPro / dist;
		}
		cout << "3 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;

	#pragma omp for reduction(+:inductance) private(i, j, dist, dotPro)
	for (int i = 0; i < step; i++) {
		for (int j = 0; j < step; j++) {
			dist = sqrt(pow(pos1_Rev[i][0] - pos2_Rev[j][0], 2.) +
						pow(pos1_Rev[i][1] - pos2_Rev[j][1], 2.) +
						pow(pos1_Rev[i][2] - pos2_Rev[j][2], 2.));
			dotPro = tang1_Rev[i][0] * tang2_Rev[j][0] +
					 tang1_Rev[i][1] * tang2_Rev[j][1] +
					 tang1_Rev[i][2] * tang2_Rev[j][2];
			inductance += dotPro / dist;
		}
		cout << "4 : " << i + 1 << "/" << step << endl;
	}
	#pragma omp single
	cout << endl;
}	*/

	processTime = clock() - startTime;
	cout << static_cast<double>(processTime) / 1000 << " [s]" << endl;

	// 確認等用csv書き出し
	{
		for (int i = 0; i < ring_sp; i++) {
			ofs_out << ind0[i][0] << ",";
			ofs_out << ind0[i][1] << ",";
			ofs_out << ind0[i][2] << ",";
			ofs_out << ind0[i][3] << endl;
		}
	}

	// ファイルクローズ・gnuplot書き出し
	ofs_out.close();
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
