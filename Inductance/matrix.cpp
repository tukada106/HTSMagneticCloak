#include <iostream>
#include <math.h>

#include "matrix.h"
using namespace std;

//---------------------------------
//     通常のコンストラクタ
//---------------------------------
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
            dpTop[i][j] = 0;
        }
    }
}

//---------------------------------
//     コピーコンストラクタ
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
//   デストラクタ
//----------------------
Matrix::~Matrix() {
    for (int i = 0; i < row; i++) {
        delete[] dpTop[i];
    }
    delete[] dpTop;
}

//------------------------------------
//     代入
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
//       行列の加算
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
//       行列の減算
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
//       行列の積
//------------------------------------
Matrix Matrix::operator*(const Matrix& mat) {
    if (column != mat.row) {
        cerr << "Matrix size mismatch!" << endl;
        cerr << "ERR : operator*" << endl;
        exit(1);
    }

    Matrix ret(row, mat.column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < mat.column; j++) {
            for (int k = 0; k < column; k++) {
                ret.dpTop[i][j] += dpTop[i][k] * mat.dpTop[k][j];
            }
        }
    }

    return ret;
}

//--------------------------------------
//       行列の定数倍
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
//  転置行列を返す
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
//  逆行列を返す
//----------------------------------------
Matrix Matrix::inverse() {
    Matrix gauss_jorudan(row, column * 2);

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

    return ret;
}
