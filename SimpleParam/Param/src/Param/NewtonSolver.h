#ifndef NEWTONSOLVER_H_
#define NEWTONSOLVER_H_

#include "../hj_3rd/include/sparse/sparse.h"

namespace PARAM
{
	class NewtonSolver
	{
	public:
		NewtonSolver(int _variable_num);
		~NewtonSolver();

		void Solve();

	private:
		void SetInitVariableValue(const std::vector<double>& _variable_value);
		void SetJacobiMatrix(const hj::sparse::spm_csc<double>& _jacobi_matrix);
		void SolveLinearEquation();
		void UpdateJacobiMatrix();
		void UpdateVariableValue();

	private:
		hj::sparse::spm_csc<double> m_jacobi_matrix;
		std::vector<double> m_variable_value;

		int m_variable_num;
	};
}

#endif