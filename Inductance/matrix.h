#pragma once
class Matrix {

	int row;		// �s
	int column;		// ��

	double** dpTop;	// �z��̍ŏ����w���|�C���^

public:
	Matrix(int i = 1, int j = 1);	// �R���X�g���N�^
	Matrix(const Matrix& copy);		// �R�s�[�R���X�g���N�^

	~Matrix();						// �f�X�g���N�^

	int row_size() { return(row); }			// �s����Ԃ�
	int column_size() { return(column); }	// �񐔂�Ԃ�

	//���Z�q�̃I�[�o�[���[�h
	double*& operator[](int i) { return(dpTop[i]); }
	Matrix operator=(const Matrix& a);
	Matrix operator+(const Matrix& a);
	Matrix operator-(const Matrix& a);
	Matrix operator*(const Matrix& a);

	friend Matrix operator*(const Matrix& a, double b);
	friend Matrix operator*(double b, const Matrix& a);

	//�s��̕ϊ��Ȃ�
	Matrix transposed();	// �]�u�s���Ԃ�
	Matrix inverse();		// �t�s���Ԃ�
};