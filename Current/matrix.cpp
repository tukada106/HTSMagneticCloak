#include <iostream>
#include <math.h>

#include "matrix.h"

using namespace std;

//---------------------------------
// コンストラクタ
//---------------------------------
Matrix::Matrix() {
    // 行・列数を保存
    row = 1;
    column = 1;

    // 配列領域の動的確保
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // 値の初期化
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = 0.;
        }
    }
}
Matrix::Matrix(int size_x) {
    // 行・列数を保存
    row = size_x;
    column = 1;

    // 配列領域の動的確保
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // 値の初期化
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = 0.;
        }
    }
}
Matrix::Matrix(int size_x, int size_y) {
    // 行・列数を保存
    row = size_x;
    column = size_y;

    // 配列領域の動的確保
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // 値の初期化
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = 0.;
        }
    }
}

//---------------------------------
// コピーコンストラクタ
//---------------------------------
Matrix::Matrix(const Matrix& copy) {
    // 行・列数を保存
    row = copy.row;
    column = copy.column;

    // 配列領域の動的確保
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // 値の代入
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = copy.dpTop[i][j];
        }
    }
}

//----------------------
// デストラクタ
//----------------------
Matrix::~Matrix() {
    for (int i = 0; i < row; i++) {
        delete[] dpTop[i];
    }
    delete[] dpTop;
}

//------------------------------------
// 代入
//------------------------------------
Matrix Matrix::operator=(const Matrix& mat) {
    if (row != mat.row || column != mat.column) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : operator=" << endl;
        exit(1);
    }

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = mat.dpTop[i][j];
        }
    }

    return *this;
}

//------------------------------------
// 行列の加算
//------------------------------------
Matrix Matrix::operator+(const Matrix& mat) {
    if (row != mat.row || column != mat.column) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : operator+" << endl;
        exit(1);
    }

    Matrix ret(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            ret.dpTop[i][j] = dpTop[i][j] + mat.dpTop[i][j];
        }
    }

    return ret;
}

//------------------------------------
// 行列の減算
//------------------------------------
Matrix Matrix::operator-(const Matrix& mat) {
    if (row != mat.row || column != mat.column) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : operator-" << endl;
        exit(1);
    }

    Matrix ret(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            ret.dpTop[i][j] = dpTop[i][j] - mat.dpTop[i][j];
        }
    }

    return ret;
}

//------------------------------------
// 行列の積
//------------------------------------
Matrix Matrix::operator*(const Matrix& mat) {
    if (column != mat.row) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : operator*" << endl;
        exit(1);
    }

    Matrix ret(row, mat.column);
    for (int i = 0; i < row; i++) {
        for (int k = 0; k < column; k++) {
            for (int j = 0; j < mat.column; j++) {
                ret.dpTop[i][j] += dpTop[i][k] * mat.dpTop[k][j];
            }
        }
    }

    return ret;
}

//--------------------------------------
// 行列の定数倍
//--------------------------------------
Matrix operator*(const Matrix& mat, double con) {
    Matrix ret(mat.row, mat.column);
    for (int i = 0; i < mat.row; i++) {
        for (int j = 0; j < mat.column; j++) {
            ret[i][j] = con * mat.dpTop[i][j];
        }
    }

    return ret;
}
Matrix operator*(double con, const Matrix& mat) {
    Matrix ret(mat.row, mat.column);
    for (int i = 0; i < mat.row; i++) {
        for (int j = 0; j < mat.column; j++) {
            ret[i][j] = con * mat.dpTop[i][j];
        }
    }

    return ret;
}

//----------------------------------------
// 転置行列を返す
//----------------------------------------
Matrix Matrix::transposed() {
    Matrix ret(column, row);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            ret[j][i] = dpTop[i][j];
        }
    }

    return ret;
}

//----------------------------------------
// 逆行列を返す
//----------------------------------------
Matrix Matrix::inverse() {
    Matrix gauss_jorudan(row, column * 2);
    Matrix ret(row, column);
    double temp;

    if (row != column) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : inverse()" << endl;
        exit(1);
    }

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            gauss_jorudan[i][j] = dpTop[i][j];
        }
    }
    for (int i = 0; i < row; i++) {
        for (int j = column; j < column * 2; j++) {
            gauss_jorudan[i][j] = 0;
        }
        gauss_jorudan[i][i + column] = 1;
    }

    for (int k = 0; k < row; k++) {
        double p = 0.000000000001;
        int c = -1;
        for (int i = k; i < row; i++) {
            if (fabs(gauss_jorudan[i][k]) > p) {
                p = fabs(gauss_jorudan[i][k]);
                c = i;
            }
        }
        if (c < 0) {
            cerr << "No inverse matrix exists!" << endl;
            cerr << "ERR : inverse()" << endl;
            exit(1);
        }
        if (c != k) {
            for (int j = k; j < column * 2; j++) {
                temp = gauss_jorudan[k][j];
                gauss_jorudan[k][j] = gauss_jorudan[c][j];
                gauss_jorudan[c][j] = temp;
            }
        }

        temp = gauss_jorudan[k][k];
        for (int j = k; j < column * 2; j++) {
            gauss_jorudan[k][j] = gauss_jorudan[k][j] / temp;
        }
        for (int i = 0; i < row; i++) {
            if (i != k) {
                temp = gauss_jorudan[i][k];
                for (int j = k; j < column * 2; j++) {
                    gauss_jorudan[i][j] = gauss_jorudan[i][j] - temp * gauss_jorudan[k][j];
                }
            }
        }
    }

    for (int i = 0; i < row; i++) {
        for (int j = column; j < column * 2; j++) {
            ret[i][j - column] = gauss_jorudan[i][j];
        }
    }

    return ret;
}

//----------------------------------------
// 引数指定の行を返す
//----------------------------------------
Matrix Matrix::extract_row(const int row_num) {
    Matrix ret(1, column);

    if (row_num > row) {
        cerr << "The specified row is missing!" << endl;
        cerr << "ERR : extract_row()" << endl;
        exit(1);
    }

    for (int i = 0; i < column; i++) {
        ret[0][i] = dpTop[row_num - 1][i];
    }

    return ret;
}

//----------------------------------------
// 引数指定の列を返す
//----------------------------------------
Matrix Matrix::extract_column(const int column_num) {
    Matrix ret(column);

    if (column_num > column) {
        cerr << "The specified column is missing!" << endl;
        cerr << "ERR : extract_column()" << endl;
        exit(1);
    }

    for (int i = 0; i < row; i++) {
        *ret[i] = dpTop[i][column_num - 1];
    }

    return ret;
}

//----------------------------------------
// ノルムを返す
//----------------------------------------
double Matrix::norm() {
    double ret = 0.;

    if (column != 1) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : norm()" << endl;
        exit(1);
    }

    for (int i = 0; i < row; i++) {
        ret += pow(dpTop[i][0], 2.);
    }
    ret = sqrt(ret);

    return ret;
}

//----------------------------------------
// 外積を返す
//----------------------------------------
Matrix cross(Matrix& const matA, Matrix& const matB) {
    Matrix ret(3);

    if (matA.row_size() != 3 || matB.row_size() != 3) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : cross()" << endl;
        exit(1);
    }
    
    ret[0][0] = matA[1][0] * matB[2][0] - matA[2][0] * matB[1][0];
    ret[1][0] = matA[2][0] * matB[0][0] - matA[0][0] * matB[2][0];
    ret[2][0] = matA[0][0] * matB[1][0] - matA[1][0] * matB[0][0];

    return ret;
}

//----------------------------------------
// 内積を返す
//----------------------------------------
double dot(Matrix& const matA, Matrix& const matB) {
    double ret = 0;

    if (matA.row_size() != matB.row_size()) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : dot()" << endl;
        exit(1);
    }

    for (int i = 0; i < matA.row_size(); i++) {
        ret += matA[i][0] * matB[i][0];
    }

    return ret;
}

//----------------------------------------
// 並列処理で行列の足し算
//----------------------------------------
int omp_plus(Matrix& ret, Matrix& const matA, Matrix& const matB, int thread_num, int n_threads) {
    if (matA.row_size() != matB.row_size() || matA.column_size() != matB.column_size()) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : omp_plus" << endl;
        exit(1);
    }
    if (ret.row_size() != matA.row_size() || ret.column_size() != matA.column_size()) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : omp_plus" << endl;
        exit(1);
    }

    int loop_max = matA.row_size();
    int loop_per_thread = loop_max / n_threads;
    int loop_mod = loop_max % n_threads;
    int start = 0, end = 0;
    if (thread_num < n_threads) {
        if (thread_num < loop_mod) {
            start = thread_num * (loop_per_thread + 1);
        }
        else {
            start = thread_num * loop_per_thread + loop_mod;
        }
    }
    else {
        start = 0;
    }
    if ((thread_num + 1) * loop_per_thread <= loop_max) {
        if (thread_num < loop_mod) {
            end = (thread_num + 1) * (loop_per_thread + 1);
        }
        else {
            if ((thread_num + 1) * loop_per_thread + loop_mod <= loop_max) {
                end = (thread_num + 1) * loop_per_thread + loop_mod;
            }
            else {
                end = 0;
            }
        }
    }
    else {
        end = 0;
    }
    cout << thread_num << "/" << n_threads << "\t" << start << "->" << end << endl;

    for (int i = start; i < end; i++) {
        for (int j = 0; j < matA.column_size(); j++) {
            ret[i][j] = matA[i][j] + matB[i][j];
        }
    }

//#pragma omp barrier

    return 0;
}

//----------------------------------------
// 並列処理で行列の掛け算
//----------------------------------------
int omp_multi(Matrix& ret, Matrix& const matA, Matrix& const matB, int thread_num, int n_threads) {
    if (matA.column_size() != matB.row_size()) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : omp_multi" << endl;
        exit(1);
    }
    if (ret.row_size() != matA.row_size() || ret.column_size() != matB.column_size()) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : omp_multi" << endl;
        exit(1);
    }
    cout << "I'm " << thread_num << endl;

    int rows = matA.row_size();
    int row_mod = rows % n_threads;
    int start_step = rows * thread_num / n_threads;
    int end_step = rows * (thread_num + 1) / n_threads - 1;
    if (thread_num < row_mod) {
        start_step += thread_num;
        end_step += (thread_num + 1);
    }
    else {
        start_step += row_mod;
        if (thread_num != n_threads - 1) {
            end_step += row_mod;
        }
        else {
            end_step = rows - 1;
        }
    }

    for (int i = start_step; i <= end_step; i++) {
        for (int k = 0; k < matA.column_size(); k++) {
            for (int j = 0; j < matB.column_size(); j++) {
                ret[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }

#pragma omp barrier

    return 0;
}

//----------------------------------------
// 並列処理で実数の掛け算
//----------------------------------------
int omp_multi(Matrix& ret, double con, Matrix& const matA, int thread_num, int n_threads) {
    if (ret.row_size() != matA.row_size() || ret.column_size() != matA.column_size()) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : omp_multi" << endl;
        exit(1);
    }

    int rows = matA.row_size();
    int row_mod = rows % n_threads;
    int start_step = rows * thread_num / n_threads;
    int end_step = rows * (thread_num + 1) / n_threads - 1;
    if (thread_num < row_mod) {
        start_step += thread_num;
        end_step += (thread_num + 1);
    }
    else {
        start_step += row_mod;
        if (thread_num != n_threads - 1) {
            end_step += row_mod;
        }
        else {
            end_step = rows - 1;
        }
    }

    for (int i = start_step; i <= end_step; i++) {
        for (int j = 0; j < matA.column_size(); j++) {
            ret[i][j] = con * matA[i][j];
        }
    }

#pragma omp barrier

    return 0;
}
int omp_multi(Matrix& ret, Matrix& const matA, double con, int thread_num, int n_threads) {
    omp_multi(ret, con, matA, thread_num, n_threads);

    return 0;
}
