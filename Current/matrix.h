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
	Matrix extract_row(const int row_num);			// 引数指定の行を返す
	Matrix extract_column(const int column_num);	// 引数指定の列を返す
	double norm();
};

Matrix cross(Matrix& const matA, Matrix& const matB);	// 外積を返す
double dot(Matrix& const matA, Matrix& const matB);		// 内積を返す

int omp_plus(Matrix& ret, Matrix& const matA, Matrix& const matB, int thread_num, int n_threads);	// 並列処理で行列の足し算
int omp_multi(Matrix& ret, Matrix& const matA, Matrix& const matB, int thread_num, int n_threads);	// 並列処理で行列の掛け算
int omp_multi(Matrix& ret, Matrix& const matA, double con, int thread_num, int n_threads);	// 並列処理で行列と実数の掛け算
int omp_multi(Matrix& ret, double con, Matrix& const matA, int thread_num, int n_threads);	// 並列処理で行列と実数の掛け算
