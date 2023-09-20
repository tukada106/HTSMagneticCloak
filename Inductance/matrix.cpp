#include <iostream>
#include <math.h>

#include "matrix.h"

using namespace std;

//---------------------------------
// �R���X�g���N�^
//---------------------------------
Matrix::Matrix() {
    // �s�E�񐔂�ۑ�
    row = 1;
    column = 1;

    // �z��̈�̓��I�m��
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // �l�̏�����
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = 0;
        }
    }
}
Matrix::Matrix(int size_x) {
    // �s�E�񐔂�ۑ�
    row = size_x;
    column = 1;

    // �z��̈�̓��I�m��
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // �l�̏�����
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = 0;
        }
    }
}
Matrix::Matrix(int size_x, int size_y) {
    // �s�E�񐔂�ۑ�
    row = size_x;
    column = size_y;

    // �z��̈�̓��I�m��
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // �l�̏�����
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = 0;
        }
    }
}

//---------------------------------
// �R�s�[�R���X�g���N�^
//---------------------------------
Matrix::Matrix(const Matrix& copy) {
    // �s�E�񐔂�ۑ�
    row = copy.row;
    column = copy.column;

    // �z��̈�̓��I�m��
    dpTop = new double* [row];
    for (int i = 0; i < row; i++) {
        dpTop[i] = new double[column];
    }

    // �l�̑��
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            dpTop[i][j] = copy.dpTop[i][j];
        }
    }
}

//----------------------
// �f�X�g���N�^
//----------------------
Matrix::~Matrix() {
    for (int i = 0; i < row; i++) {
        delete[] dpTop[i];
    }
    delete[] dpTop;
}

//------------------------------------
// ���
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
// �s��̉��Z
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
// �s��̌��Z
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
// �s��̐�
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
// �s��̒萔�{
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
// �]�u�s���Ԃ�
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
// �t�s���Ԃ�
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
// �O�ς�Ԃ�
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
// ���ς�Ԃ�
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
