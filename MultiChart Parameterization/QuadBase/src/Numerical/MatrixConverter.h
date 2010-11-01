#ifndef MATRIX_CONVERTER_2_H
#define MATRIX_CONVERTER_2_H


#include "MeshSparseMatrix.h"
#include "matrix.h"
#include "../hj_3rd/include/sparse/sparse.h"
using namespace zjucad::matrix;

class CMatrixConverter
{
public:
	CMatrixConverter();
	~CMatrixConverter();

public:
	static void CSparseMatrix2hjMatrix(matrix<double>& hjMatrix, CMeshSparseMatrix& csMatrix);
	static void HjMatrix2CMatrix(matrix<double>& hjMatrix, CMatrix& cMatrix);
	static void CSparseMatrix2hjCscMatrix(hj::sparse::spm_csc<double>& hjcscMatrix, CMeshSparseMatrix& cMatrix);
	static void CMatrix2hjMatrix(matrix<double>& hjMatrix, CMatrix& cMatrix);
};

#endif