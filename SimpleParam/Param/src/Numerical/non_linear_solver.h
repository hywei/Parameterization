#ifndef NON_LINEAR_SOLVER_H
#define NON_LINEAR_SOLVER_H

#include "MeshSparseMatrix.h"
#include "solver.h"
#include "../Graphite/OGF/math/symbolic/symbolic.h"
#include "../Graphite/OGF/math/symbolic/stencil.h"
#include <vector>
#include <deque>

#ifdef WIN32
#include <hj_3rd/hjlib/sparse_old/sparse.h>
#else
#include <hj_3rd/hjlib/sparse/sparse.h>
#endif

enum NonLinearSolveMethod {
	GAUSS_NEWTON, 
	QUASI_NEWTON,
	NEWTON,
	LEVENBERG_MARQUARDT,
	LAGRANGE_GAUSS_NEWTON
} ;

class NonLinearSolver : public Solver 
{

public:
	NonLinearSolver(int nb_variables) ;
	virtual ~NonLinearSolver() ;

public:
	void set_solve_method(NonLinearSolveMethod method_);
	void is_printf_info(bool is_printf){m_is_printf=is_printf;}
	// __________________ Construction _____________________

	void begin_equation() ;

	/** 
	* declares a stencil from its expression. If more performance is
	* needed, use declare_stencil(Stencil* S).
	*/
	int declare_stencil(const OGF::Symbolic::Expression& f) ;

	/** 
	* declares a specialized stencil. Client code needs to derive
	* a Stencil class, and implement evaluation, gradient and Hessian.
	* If performance is not critical, use 
	* declare_stencil(const Expression& f) instead.
	*/
	int declare_stencil(OGF::Stencil* S) ;


	void begin_stencil_instance(int stencil_id) ;
	void stencil_variable(int global_variable_id) ;
	void stencil_variable(int local_variable_id, int global_variable_id) ;
	void stencil_parameter(double value) ;
	void stencil_parameter(int local_parameter_id, double value) ;
	void end_stencil_instance() ;

	void end_equation() ;


	// add by zmy
	void clear_equation() ;
	Linear_Constraint& get_linear_constraint(){return m_linear_cons;}
	void set_init_variables_value(vector<double>& init_val_vec);
	double equations_value(vector<double>& var_val_vec);
	void set_equation_div_flag();

	OGF::Stencil* stencil(int id) { 
		ogf_assert(id >= 0 && id < int(stencil_.size())) ;
		return stencil_[id] ;  
	}

	void solve();

private:
	OGF::StencilInstance& cur_stencil_instance() {
		return *(m_equation_vec.rbegin()) ;
	}

	bool is_free(int id)   { return (id < nb_free_variables_) ;  }
	bool is_locked(int id) { return (id >= nb_free_variables_) ; }

private:
	void solve_one_iteration_gaussian_newton();
	void instanciate_gaussian_newton();

	void solve_one_iteration_newton();
	void instanciate_newton();

	void solve_one_Levenberg_marquardt(double& mu_, double& nu_);

	void solve_one_iteration_Lagrange_gaussian_newton();

	void update_variables();

	double f() ;
	double f(vector<double>& xc_);
	double norm_grad_f();
	double vec_multiply_vec(vector<double>& vec_1, vector<double>& vec_2);

private:
	// User representation
	int nb_free_variables_ ;
	int nb_locked_variables_ ;
	std::vector<OGF::Stencil*> stencil_ ;

	// State
	enum State {INITIAL, IN_EQN, IN_STENCIL_REF, CONSTRUCTED, MINIMIZED} ;
	State state_ ;
	int cur_stencil_var_ ;
	int cur_stencil_prm_ ;

	// Internal representation: problem setting
	std::deque<OGF::StencilInstance> m_equation_vec;
	Linear_Constraint m_linear_cons;

	vector<double> m_dx_ ;             // Unknown delta vector for the variables
	vector<double> m_gradient_ ;               // -gradient
	CMeshSparseMatrix m_Hessian_ ;    // Hessian 
	CMeshSparseMatrix m_Jacobi_;
	vector<double> m_xc_ ;              // Variables + constants
	double fk_ ;             // value of the function at current step
	double gk_ ;            // norm of the gradient at current step
	bool m_jacobi_ata_first_time;
	double m_mu_;
	double m_nu_;

	// Tuning
	int max_newton_iter_ ;
	double gradient_threshold_ ;

	// Solve method
	NonLinearSolveMethod solve_method_;
	bool m_stencil_use_hessian;

	//
	bool m_is_printf;
	vector<size_t> m_equ_div_flag_vec;
};

#endif
