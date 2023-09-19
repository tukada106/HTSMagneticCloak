#pragma once
class Matrix {

	int row;		// 行
	int column;		// 列

	double** dpTop;	// 配列の最初を指すポインタ

public:
	Matrix(int i = 1, int j = 1);	// コンストラクタ
	Matrix(const Matrix& copy);		// コピーコンストラクタ

	~Matrix();						// デストラクタ

	int row_size() { return(row); }			// 行数を返す
	int column_size() { return(column); }	// 列数を返す

	//演算子のオーバーロード
	double*& operator[](int i) { return(dpTop[i]); }
	Matrix operator=(const Matrix& a);
	Matrix operator+(const Matrix& a);
	Matrix operator-(const Matrix& a);
	Matrix operator*(const Matrix& a);

	friend Matrix operator*(const Matrix& a, double b);
	friend Matrix operator*(double b, const Matrix& a);

	//行列の変換など
	Matrix transposed();	// 転置行列を返す
	Matrix inverse();		// 逆行列を返す
};