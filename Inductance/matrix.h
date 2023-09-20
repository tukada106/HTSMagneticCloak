#pragma once
class Matrix {

	int row;		// 行
	int column;		// 列

	double** dpTop;	// 配列の最初を指すポインタ

public:
	Matrix();						// スカラーのコンストラクタ
	Matrix(int size_x);				// ベクトルのコンストラクタ
	Matrix(int size_x, int size_y);	// 行列のコンストラクタ
	Matrix(const Matrix& copy);		// コピーコンストラクタ

	~Matrix();						// デストラクタ

	int row_size() const { return(row); }			// 行数を返す
	int column_size() const { return(column); }		// 列数を返す

	//演算子のオーバーロード
	double*& operator[](int i) { return(dpTop[i]); }
	Matrix operator=(const Matrix& mat);
	Matrix operator+(const Matrix& mat);
	Matrix operator-(const Matrix& mat);
	Matrix operator*(const Matrix& mat);

	friend Matrix operator*(const Matrix& mat, double con);
	friend Matrix operator*(double con, const Matrix& mat);

	//行列の変換など
	Matrix transposed();	// 転置行列を返す
	Matrix inverse();		// 逆行列を返す
};

Matrix cross(Matrix& const matA, Matrix& const matB);	// 外積を返す
double dot(Matrix& const matA, Matrix& const matB);		// 内積を返す
