#include "MatrixConverter.h"

//////////////////////////////////////////////////////////////////////
// CMatrixConverter Construction/Destruction
//////////////////////////////////////////////////////////////////////
CMatrixConverter::CMatrixConverter()
{
}
CMatrixConverter::~CMatrixConverter()
{
}

//////////////////////////////////////////////////////////////////////
// CMatrixConverter public methods
//////////////////////////////////////////////////////////////////////
void CMatrixConverter::CSparseMatrix2hjMatrix(matrix<double>& hjMatrix, CMeshSparseMatrix& csMatrix)
{
	size_t row = csMatrix.GetRowNum();
	size_t col = csMatrix.GetColNum();
	
	hjMatrix = zeros<double>(row, col);

	for (size_t i = 0; i < csMatrix.m_ColIndex.size(); i++)
	{
		IndexArray& colIndex = csMatrix.m_ColIndex[i];

		for (size_t j = 0; j < colIndex.size(); j++)
		{
			hjMatrix(colIndex[j], i) = csMatrix.m_ColData[i][j];
		}
	}

// 	for (size_t i = 0; i < row; i++)
// 	{
// 		for (size_t j = 0; j < col; j++)
// 		{
// 			csMatrix.GetElement(i, j, value);
// 			hjMatrix(i, j) = value;
// 		}
// 	}
}
void CMatrixConverter::HjMatrix2CMatrix(matrix<double>& hjMatrix, CMatrix& cMatrix)
{
	size_t row = hjMatrix.size(1);
	size_t col = hjMatrix.size(2);

	cMatrix.Init(row, col);

	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < col; j++)
		{
			cMatrix.SetElement(i, j, hjMatrix(i, j));
		}
	}
}
void CMatrixConverter::CSparseMatrix2hjCscMatrix(hj::sparse::spm_csc<double>& hjcscMatrix, CMeshSparseMatrix& cMatrix)
{
	/*
	matrix<double> hjMatrix;
	hj::sparse::spm_csc<double> cscMatrix;

	CSparseMatrix2hjMatrix(hjMatrix, cMatrix);

	hj::sparse::convert(hjMatrix, hjcscMatrix, SMALL_ZERO_EPSILON);

	return;*/

	//
	cMatrix.MakeMatrixIndexLessSeq();

	int i;
	int m_Row = cMatrix.GetRowNum();
	int m_Col = cMatrix.GetColNum();

	int nnz = cMatrix.GetNZNum();
	hjcscMatrix.resize(m_Row, m_Col, nnz);

	//
	int index = 0;
	double value = 0;

	for(i = 0; i < m_Col; ++ i)
	{
		(hjcscMatrix.ptr_)[i] = index;

		size_t n = cMatrix.m_ColIndex[i].size();

		for(size_t j = 0; j < n; ++ j)
		{
			(hjcscMatrix.idx_)[index] = cMatrix.m_ColIndex[i][j];
			(hjcscMatrix.val_)[index ++] = cMatrix.m_ColData[i][j];
		}
	}
	hjcscMatrix.ptr_[i] = nnz;

	/*
	// debug.
	matrix<double> hjMatrix;
	hj::sparse::spm_csc<double> cscMatrix;

	CSparseMatrix2hjMatrix(hjMatrix, cMatrix);

	hj::sparse::convert(hjMatrix, cscMatrix, SMALL_ZERO_EPSILON);

	for(size_t i = 0; i < hjcscMatrix.idx_.size(); i++)
	{
		ASSERT(cscMatrix.idx_[i] == hjcscMatrix.idx_[i]);
	}
	for(size_t i = 0; i < hjcscMatrix.ptr_.size(); i++)
	{
		ASSERT(cscMatrix.ptr_[i] == hjcscMatrix.ptr_[i]);
	}
	for(size_t i = 0; i < hjcscMatrix.val_.size(); i++)
	{
		ASSERT(cscMatrix.val_[i] == hjcscMatrix.val_[i]);
	}*/
}
void CMatrixConverter::CMatrix2hjMatrix(matrix<double>& hjMatrix, CMatrix& cMatrix)
{
	size_t row = cMatrix.GetNumRows();
	size_t col = cMatrix.GetNumColumns();

	hjMatrix.resize(row, col);

	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < col; j++)
		{
			hjMatrix(i, j) = cMatrix.GetElement(i, j);
		}
	}
}