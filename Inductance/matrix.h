#pragma once
class Matrix {

	int row;		// �s
	int column;		// ��

	double** dpTop;	// �z��̍ŏ����w���|�C���^

public:
	Matrix();						// �X�J���[�̃R���X�g���N�^
	Matrix(int size_x);				// �x�N�g���̃R���X�g���N�^
	Matrix(int size_x, int size_y);	// �s��̃R���X�g���N�^
	Matrix(const Matrix& copy);		// �R�s�[�R���X�g���N�^

	~Matrix();						// �f�X�g���N�^

	int row_size() const { return(row); }			// �s����Ԃ�
	int column_size() const { return(column); }		// �񐔂�Ԃ�

	//���Z�q�̃I�[�o�[���[�h
	double*& operator[](int i) { return(dpTop[i]); }
	Matrix operator=(const Matrix& mat);
	Matrix operator+(const Matrix& mat);
	Matrix operator-(const Matrix& mat);
	Matrix operator*(const Matrix& mat);

	friend Matrix operator*(const Matrix& mat, double con);
	friend Matrix operator*(double con, const Matrix& mat);

	//�s��̕ϊ��Ȃ�
	Matrix transposed();	// �]�u�s���Ԃ�
	Matrix inverse();		// �t�s���Ԃ�
};

Matrix cross(Matrix& const matA, Matrix& const matB);	// �O�ς�Ԃ�
double dot(Matrix& const matA, Matrix& const matB);		// ���ς�Ԃ�
