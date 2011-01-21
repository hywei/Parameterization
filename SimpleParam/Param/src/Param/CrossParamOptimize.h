#ifndef CROSSPARAMOPTIMIZE_H_
#define CROSSPARAMOPTIMIZE_H_

#include "CrossParam.h"
#include "../Numerical/MeshSparseMatrix.h"
#include "../hj_3rd/include/zjucad/matrix/matrix.h"

class LinearSolver;

namespace PARAM
{	
	class CrossParamOptimizer
	{
	public:
		CrossParamOptimizer(CrossParam& _cross_param) : m_cross_param(_cross_param){}
		~CrossParamOptimizer(){}

		virtual void Optimize() = 0;

	protected:
		CrossParam& m_cross_param;
	};

	class CrossParamHarmonicOptimizer : public CrossParamOptimizer
	{
	public:
		CrossParamHarmonicOptimizer(CrossParam& _cross_param) ;
		~CrossParamHarmonicOptimizer();

		void Optimize();
		

	private:
		void ComputeParamJacobiMatrix();
		bool ComputeVertexJacobiMatrixOnSurface2(const SurfaceCoord& surface_coord, zjucad::matrix::matrix<double>& j_mat);

		/// newton method
		void SetInitVariableValue();
		void SetNewtonMethodJacobiMatrix();
		void SolveJacobiLinearEquation();
		void ComputeNewtonMethodFuncValue();

		bool CheckNewtonMethodResult(const LinearSolver& linear_solver);
		bool CheckNewtonMethodJacobiMatrix();
		double ComputeNewtonMethodError();

		void VertexRelaxtion();
		void ResignParamCoord();

	private:
		std::vector<bool> m_fixed_vtx_flag;
		std::vector<zjucad::matrix::matrix<double> > m_param_jacobi_matrix;

		std::vector<double> m_variable_value;
		std::vector<int> m_vtx_var_index_mapping;
		std::vector<int> m_var_vtx_index_mapping;
		std::vector<double> m_func_value;

		int m_vari_num;

		CMeshSparseMatrix m_nw_jacobi_mat;
	};
}

#endif