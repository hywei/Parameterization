#include "non_linear_solver.h"
#include "../Common/stopwatch.h"
#include "../Numerical/MatrixConverter.h"

#ifdef WIN32
#include <hj_3rd/hjlib/sparse_old/sparse_multi_cl.h>
#else
#include <hj_3rd/hjlib/sparse/sparse_multi_cl.h>
#endif

#include <math.h>
#include <float.h>
#include <memory>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
static mm_rtn cache;

NonLinearSolver::NonLinearSolver(int nb_variables) 
: Solver(nb_variables) 
{
	state_ = INITIAL ;
	solve_method_ = GAUSS_NEWTON;
	max_newton_iter_ = 15 ;
	gradient_threshold_ = 1e-3 ;
	m_stencil_use_hessian = false;
	m_is_printf = true;

	::OGF::Logger::initialize();
}
NonLinearSolver::~NonLinearSolver() 
{
	for(size_t i=0; i<stencil_.size(); i++) 
	{
		delete stencil_[i] ;
	}
	::OGF::Logger::terminate();
}
void NonLinearSolver::set_solve_method(NonLinearSolveMethod method_)
{
	solve_method_ = method_;
	if (method_ == NEWTON)
	{
		m_stencil_use_hessian = true;
	}
	else if (method_ == LEVENBERG_MARQUARDT)
	{
		m_mu_ = -1.0;
		m_nu_ = 2.0;
	}
}
//////////////////////////////////////////////////////////////////////
// NonLinearSolver public methods
//////////////////////////////////////////////////////////////////////
void NonLinearSolver::begin_equation() 
{
	ogf_assert(state_ == INITIAL) ;
	state_ = IN_EQN ;

	nb_free_variables_ = 0 ;
	nb_locked_variables_ = 0 ;
	int cur_index = 0 ;

	// check if referred variable is also setted value
	int valued_var_num;
	do 
	{
		valued_var_num = 0;
		for(int i=0; i<nb_variables_; i++) {
			if (variable_[i].is_referred() && variable_[i].is_valued()) {
				int ref_var_index = variable_[i].ref_var_index();
				if (variable_[ref_var_index].is_valued()) {
					double val_ = variable_[ref_var_index].value();
					if (val_ != variable_[i].value()) {
						printf("%d and %d variable value are conflicted.\n", i, ref_var_index);
					}
				} else {
					variable_[ref_var_index].lock();
					variable_[ref_var_index].set_value(variable_[i].value());
					valued_var_num++;
				}
			}
		}
	} while(valued_var_num > 0);
	for(int i=0; i<nb_variables_; i++) {
		if (variable_[i].is_referred() && variable_[i].is_valued()) {
			variable_[i].set_ref_variable_inex(-1);
		}
	}

	for(int i=0; i<nb_variables_; i++) {
		if(!variable_[i].is_locked()) {
			set_variable_index(variable_[i], cur_index) ;
			nb_free_variables_++ ;
			cur_index++ ;
		}
	}

	//
	for(int i=0; i<nb_variables_; i++) {
		if(variable_[i].is_valued()) {
			set_variable_index(variable_[i], cur_index) ;
			cur_index++ ;
			nb_locked_variables_++ ;
		}
	}
	for(int i=0; i<nb_variables_; i++) {
		if (variable_[i].is_referred()) {
			int ref_var_index = variable_[i].ref_var_index();
			if (variable_[ref_var_index].is_valued()) {
				variable_[i].set_value(variable_[ref_var_index].value());
				set_variable_index(variable_[i], cur_index) ;
				cur_index++ ;
				nb_locked_variables_++ ;
			} else {
				set_variable_index(variable_[i], variable_[ref_var_index].index());
				variable_[i].unlock();
			}
		}
	}

	/*
	for(int i=0; i<nb_variables_; i++) {
		if(!variable_[i].is_locked()) {
			set_variable_index(variable_[i],cur_index) ;
			nb_free_variables_++ ;
			cur_index++ ;
		}
	}
	for(int i=0; i<nb_variables_; i++) {
		if(variable_[i].is_locked()) {
			set_variable_index(variable_[i],cur_index) ;
			nb_locked_variables_++ ;
			cur_index++ ;
		}
	}*/

	m_equ_div_flag_vec.push_back(0);
}
void NonLinearSolver::end_equation() 
{
	ogf_assert(state_ == IN_EQN) ;
	state_ = CONSTRUCTED ;

	//
	m_dx_.resize(nb_free_variables_) ;
	m_gradient_.resize(nb_free_variables_) ;
	m_xc_.resize(nb_variables_) ;
	for(int i=0; i<nb_variables_; i++) {
		m_xc_[variable_[i].index()] = variable_[i].value() ;
	}

	m_equ_div_flag_vec.push_back(m_equation_vec.size());

	//m_Hessian_.SetRowCol(nb_free_variables_, nb_free_variables_);
}
void NonLinearSolver::clear_equation() 
{
	m_equation_vec.clear();
	m_equ_div_flag_vec.clear();
	for(int i=0; i<nb_variables_; i++) {
		variable_[i].set_value(0.0);
	}

	state_ = INITIAL;
}
int NonLinearSolver::declare_stencil(OGF::Stencil* S) {
	stencil_.push_back(S) ;
	return (int) (stencil_.size() - 1) ;
}

int NonLinearSolver::declare_stencil(const OGF::Symbolic::Expression& f) {
	OGF::Symbolic::Stencil* S = new OGF::Symbolic::Stencil(f, m_stencil_use_hessian) ;
	
	if (m_is_printf){
		std::cout << "Stencil " << stencil_.size() << std::endl ;
		std::cout << "------------------------" << std::endl ;
		S->print(std::cout) ;
	}

	return declare_stencil(S) ;
}
void NonLinearSolver::begin_stencil_instance(int stencil_id) {
	ogf_assert(state_ == IN_EQN) ;
	state_ = IN_STENCIL_REF ;
	ogf_debug_assert(stencil_id >= 0 && stencil_id < int(stencil_.size())) ;
	cur_stencil_var_ = 0 ;
	cur_stencil_prm_ = 0 ;

	m_equation_vec.push_back(OGF::StencilInstance()) ;
	m_equation_vec.back().bind(stencil_[stencil_id]) ;
}

void NonLinearSolver::end_stencil_instance() {
	ogf_debug_assert(cur_stencil_instance().is_initialized()) ;
	ogf_assert(state_ == IN_STENCIL_REF) ;
	state_ = IN_EQN ;
}
void NonLinearSolver::stencil_variable(int global_variable_id) {
	stencil_variable(cur_stencil_var_, global_variable_id) ;
	cur_stencil_var_++ ;
}

void NonLinearSolver::stencil_variable(int local_variable_id, int global_variable_id) {
	ogf_assert(state_ == IN_STENCIL_REF) ;
	ogf_debug_assert(local_variable_id < cur_stencil_instance().nb_variables()) ;
	ogf_debug_assert(global_variable_id >= 0 && global_variable_id < nb_variables_) ;
	int internal_id = variable_[global_variable_id].index() ;
	cur_stencil_instance().global_variable_index(local_variable_id) = internal_id ;
}
void NonLinearSolver::stencil_parameter(double value) {
	stencil_parameter(cur_stencil_prm_, value) ;
	cur_stencil_prm_++ ;
}

void NonLinearSolver::stencil_parameter(int local_prm_id, double value) {
	ogf_assert(state_ == IN_STENCIL_REF) ;
	ogf_debug_assert(local_prm_id < cur_stencil_instance().nb_parameters()) ;
	cur_stencil_instance().parameter(local_prm_id) = value ;
}
void NonLinearSolver::set_init_variables_value(vector<double>& init_val_vec)
{
	for(int i=0; i<nb_variables_; i++) {
		m_xc_[variable_[i].index()] = init_val_vec[i] ;
	}
}
double NonLinearSolver::equations_value(vector<double>& var_val_vec)
{
	set_init_variables_value(var_val_vec);
	return f();
}
void NonLinearSolver::set_equation_div_flag()
{
	m_equ_div_flag_vec.push_back(m_equation_vec.size());
}
//////////////////////////////////////////////////////////////////////
// SH_Func_Rot private methods for solving
//////////////////////////////////////////////////////////////////////
void NonLinearSolver::solve()
{
	ogf_assert(state_ == CONSTRUCTED) ;

	int k = 0;
	m_jacobi_ata_first_time = true;
	for(; k<max_newton_iter_; k++) {

		//
		if (solve_method_ == GAUSS_NEWTON)
		{
			solve_one_iteration_gaussian_newton();
		}
		else if (solve_method_ == QUASI_NEWTON)
		{
		}
		else if (solve_method_ == NEWTON)
		{
			solve_one_iteration_newton();
		}
		else if (solve_method_ == LEVENBERG_MARQUARDT)
		{
			solve_one_Levenberg_marquardt(m_mu_, m_nu_);
		}
		else if (solve_method_ == LAGRANGE_GAUSS_NEWTON)
		{
			m_jacobi_ata_first_time = true;
			solve_one_iteration_Lagrange_gaussian_newton();
		}

		//
		gk_ = norm_grad_f() ;

		if (m_is_printf){
			printf("gradient value: %f\n", gk_);
		}
		
		if(gk_ < gradient_threshold_) {
			break ;
		}
	}
	printf("iteration %d times.\n", k);

	update_variables() ;
	state_ = MINIMIZED ;
}
void NonLinearSolver::solve_one_Levenberg_marquardt(double& mu_, double& nu_)
{
	//
	instanciate_gaussian_newton();

	//
	SystemStopwatch watch_;
	CMeshSparseMatrix JTJ;
	{
		CMeshSparseMatrix jacobian_matrix_AT;
		m_Jacobi_.Transpose(jacobian_matrix_AT);
		jacobian_matrix_AT.Multiply(m_Jacobi_, JTJ);
	}

	// if first time
	double avg_ii;
	if (mu_ < 0)
	{
		mu_ = JTJ.MaxAii() * 1e-3;
		avg_ii = JTJ.AvgAii();
	}

	bool is_step_ok = false;
	while (!is_step_ok)
	{
		// add \mu*I here.
		for (int i = 0; i < JTJ.GetRowNum(); i++)
		{
			JTJ.AddElement(i, i, mu_);
		}
		hj::sparse::spm_csc<double> spm_ATA;
		CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_ATA, JTJ);

		std::auto_ptr<hj::sparse::solver> m_solver_;
		m_solver_.reset(hj::sparse::solver::create(spm_ATA, "cholmod"));

		bool su = m_solver_->solve(&m_gradient_[0], &m_dx_[0]);
		if (!su)
		{
			printf("solve dx failed.\n");
		}

		// find \mu here.
		double F_xold = f();
		vector<double> tmp_mx_(m_xc_);
		for(int i=0; i<nb_free_variables_; i++) {
			tmp_mx_[i] -= m_dx_[i] ;
		}

		// update \mu and \nu
		double F_xnew = f(tmp_mx_);
		double F_diff = F_xold - F_xnew;

		double L_diff = 0; 
		vector<double> L_val_vec(m_gradient_.size());
		for (size_t i = 0; i < L_val_vec.size(); i++)
		{
			L_val_vec[i] = -mu_*m_dx_[i] - m_gradient_[i];
		}
		for (size_t i = 0; i < L_val_vec.size(); i++)
		{
			L_diff += (-m_dx_[i] * L_val_vec[i]);
		}
		L_diff /= 2.0;

		double rho_ = F_diff / L_diff;

		if (rho_ > 0 )
		{
			mu_ *= std::max(1.0/ 3.0, (1.0 - pow((2*rho_ - 1.0), 3.0)));
			nu_ = 2.0;
			
			// update
			for(int i=0; i<nb_free_variables_; i++) {
				m_xc_[i] -= m_dx_[i] ;
			}
			is_step_ok = true;
		}
		else
		{
			// remove \mu*I here.
			for (int i = 0; i < JTJ.GetRowNum(); i++)
			{
				JTJ.AddElement(i, i, -mu_);
			}

			mu_ *= nu_;
			nu_ *= 2.0;
		}
	}
	watch_.print_elapsed_time();
}
void NonLinearSolver::solve_one_iteration_Lagrange_gaussian_newton()
{
	//
	instanciate_gaussian_newton();

	// H = JT * J
	SystemStopwatch watch_;
	hj::sparse::spm_csc<double> spm_ATA;
	{
		hj::sparse::spm_csc<double> spm_A;
		CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_A, m_Jacobi_);

		if(m_jacobi_ata_first_time) 
		{
			cache = spm_dmm(true, spm_A, false, spm_A);
			m_jacobi_ata_first_time = false;
		}
		spm_dmm(true, spm_A, false, spm_A, spm_ATA, &cache);
	}

	// solve 
	std::auto_ptr<hj::sparse::solver> m_solver_;
	m_solver_.reset(hj::sparse::solver::create(spm_ATA, "cholmod"));

	//
	vector<double> Winv_Gfk(nb_free_variables_);
	bool su = m_solver_->solve(&m_gradient_[0], &Winv_Gfk[0]);
	if (!su) {
		printf("solve W-1_Gfk failed.\n");
	}

	//
	vector<double> Ak_(nb_free_variables_);
	fill(Ak_.begin(), Ak_.end(), 0.0);
	double ck = m_linear_cons.right_b();
	vector<pair<int, double> > cons_coef_vec = m_linear_cons.get_constraint_coef();
	for (size_t i = 0; i < cons_coef_vec.size(); i++)
	{
		pair<int, double>& coef_ = cons_coef_vec[i];
		int gi = variable_[coef_.first].index() ;
		if(is_free(gi)) {
			Ak_[gi] = coef_.second;
		}
		ck += m_xc_[gi]*coef_.second;
	}
	vector<double> Winv_AkT(nb_free_variables_);
	su = m_solver_->solve(&Ak_[0], &Winv_AkT[0]);
	if (!su) {
		printf("solve W-1_AkT failed.\n");
	}

	//
	double AkWinv_Gfk = vec_multiply_vec(Ak_, Winv_Gfk);
	double AkWinv_AkT = vec_multiply_vec(Ak_, Winv_AkT);
	double lambda_k_1 = (AkWinv_Gfk - ck) / AkWinv_AkT;

	//
	vector<double> AkT_lambda_Gfk(Ak_.size());
	for (size_t i = 0; i < Ak_.size(); i++)
	{
		AkT_lambda_Gfk[i] = Ak_[i]*lambda_k_1 - m_gradient_[i];
	}

	su = m_solver_->solve(&AkT_lambda_Gfk[0], &m_dx_[0]);
	if (!su) {
		printf("solve W-1_(AkT*lambda-Gfk) failed.\n");
	}

	//
	watch_.print_elapsed_time();

	// update
	for(int i=0; i<nb_free_variables_; i++) {
		m_xc_[i] += m_dx_[i] ;
	}
}
void NonLinearSolver::solve_one_iteration_gaussian_newton()
{
	//
	instanciate_gaussian_newton();

	
	// H = JT * J
	SystemStopwatch watch_;
	hj::sparse::spm_csc<double> spm_ATA;
	{
		hj::sparse::spm_csc<double> spm_A;
		CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_A, m_Jacobi_);

		if(m_jacobi_ata_first_time) 
		{
			cache = spm_dmm(true, spm_A, false, spm_A);
			m_jacobi_ata_first_time = false;
		}
		spm_dmm(true, spm_A, false, spm_A, spm_ATA, &cache);
	}

	// solve 
	std::auto_ptr<hj::sparse::solver> m_solver_;
	m_solver_.reset(hj::sparse::solver::create(spm_ATA, "cholmod"));

	bool su = m_solver_->solve(&m_gradient_[0], &m_dx_[0]);
	watch_.print_elapsed_time();
	if (!su)
	{
		printf("solve dx failed.\n");
	}

	// update
	for(int i=0; i<nb_free_variables_; i++) {
		m_xc_[i] -= m_dx_[i] ;
	}
}
void NonLinearSolver::instanciate_gaussian_newton()
{
	// if empty, init the jacobi matrix.
	if (m_Jacobi_.GetRowNum() == 0)
	{
		int row_size = (int) m_equation_vec.size();
		int col_size = nb_free_variables_;
		m_Jacobi_.SetRowCol(row_size, col_size);
	}
	vector<double> function_vector(m_Jacobi_.GetRowNum());

	// fill the jacobi matrix
	for (size_t j = 0; j < m_equation_vec.size(); j++)
	{
		OGF::StencilInstance& RS = m_equation_vec[j];
		int row_ = (int) j;
		OGF::Stencil* S = RS.stencil() ;
		int N = RS.nb_variables() ;
		OGF::Symbolic::Context args ;

		for(int i=0; i<N; i++) {
			args.variables.push_back(m_xc_[RS.global_variable_index(i)]) ;
		}
		for(int i=0; i<S->nb_parameters(); i++) {
			args.parameters.push_back(RS.parameter(i)) ;
		}

		for(int i=0; i<N; i++) {
			int gi = RS.global_variable_index(i) ;
			if(is_free(gi)) {
				m_Jacobi_.SetElement(row_, gi, S->g(i,args));
			}
		}
		function_vector[j] = S->f(args);
	}

	//  g = JT*f
	CMeshSparseMatrix jacobian_matrix_AT;
	m_Jacobi_.Transpose(jacobian_matrix_AT);
	jacobian_matrix_AT.MultiplyVector(function_vector, m_gradient_);

	//
	double f2_sum_sum = 0;
	for (size_t i = 1; i < m_equ_div_flag_vec.size(); i++)
	{
		size_t start_equ = m_equ_div_flag_vec[i-1];
		size_t end_equ = m_equ_div_flag_vec[i];

		double f2_sum = 0;
		for (size_t j = start_equ; j < end_equ; j++)
		{
			f2_sum += function_vector[j] * function_vector[j];
		}
		f2_sum_sum += f2_sum;
		printf("%d part: %f  ", i-1, f2_sum);
	}
	printf("\nfunction value: %f \n", f2_sum_sum);
}
void NonLinearSolver::solve_one_iteration_newton()
{
	//
	instanciate_newton();

	//
	SystemStopwatch watch_;
	
	hj::sparse::spm_csc<double> spm_ATA;
	CMatrixConverter::CSparseMatrix2hjCscMatrix(spm_ATA, m_Hessian_);

	// solve 
	std::auto_ptr<hj::sparse::solver> m_solver_;
	m_solver_.reset(hj::sparse::solver::create(spm_ATA, "cholmod"));

	bool su = m_solver_->solve(&m_gradient_[0], &m_dx_[0]);
	watch_.print_elapsed_time();
	if (!su)
	{
		printf("solve dx failed.\n");
	}

	double f_xold = f();
	// update
	for(int i=0; i<nb_free_variables_; i++) {
		m_xc_[i] -= m_dx_[i] ;
	}

	//double f_xnew = f();
	//double f_diff = f_xold - f_xnew;
}
void NonLinearSolver::instanciate_newton()
{
	// if empty, init the hessian matrix.
	if (m_Hessian_.GetRowNum() == 0)
	{
		int row_size = nb_free_variables_;
		int col_size = nb_free_variables_;
		m_Hessian_.SetRowCol(row_size, col_size);
	}

	// fill the hessian matrix
	for (size_t k = 0; k < m_equation_vec.size(); k++)
	{
		OGF::StencilInstance& RS = m_equation_vec[k];
		int row_ = (int) k;
		OGF::Stencil* S = RS.stencil() ;
		int N = RS.nb_variables() ;
		OGF::Symbolic::Context args ;

		for(int i=0; i<N; i++) {
			args.variables.push_back(m_xc_[RS.global_variable_index(i)]) ;
		}
		for(int i=0; i<S->nb_parameters(); i++) {
			args.parameters.push_back(RS.parameter(i)) ;
		}

		//
		for(int i=0; i<N; i++) {
			int gi = RS.global_variable_index(i) ;
			if(is_free(gi)) {
				//
				m_gradient_[gi] += S->g(i,args) ;

				//
				for(int j=0; j<=i; j++) {
					int gj = RS.global_variable_index(j) ;
					if(is_free(gj)) {
						double gij = S->G(i,j,args) ;
						if(gij != 0.0) {
							m_Hessian_.AddElement(gi, gj, gij) ;

							// G (Hessian) is symmetric, but watch out:
							// do not add diagonal elements twice !!!
							if(gi != gj) {
								m_Hessian_.AddElement(gj, gi, gij) ;
							}
						}
					}
				}
			}
		}
	}
}
void NonLinearSolver::update_variables()
{
	for(int i=0; i<nb_variables_; i++) {
		variable_[i].set_value(m_xc_[variable_[i].index()]) ;
	}
}
double NonLinearSolver::norm_grad_f()
{
	double result = 0 ;
	for(int i=0; i<nb_free_variables_; i++) {
		result += (m_gradient_[i] * m_gradient_[i]) ;
	}
	return ::sqrt(result) ;
}
double NonLinearSolver::f()
{
	return f(m_xc_);
}
double NonLinearSolver::f(vector<double>& xc_)
{
	double rst_ = 0;
	for (size_t k = 0; k < m_equation_vec.size(); k++)
	{
		OGF::StencilInstance& RS = m_equation_vec[k];
		int row_ = (int) k;
		OGF::Stencil* S = RS.stencil() ;
		int N = RS.nb_variables() ;
		OGF::Symbolic::Context args ;

		for(int i=0; i<N; i++) {
			args.variables.push_back(xc_[RS.global_variable_index(i)]) ;
		}
		for(int i=0; i<S->nb_parameters(); i++) {
			args.parameters.push_back(RS.parameter(i)) ;
		}
		double val_ = S->f(args) ;

		//
		if (solve_method_ == GAUSS_NEWTON ||
			solve_method_ == LEVENBERG_MARQUARDT)
		{
			rst_ += val_ * val_;
		}
		else
		{
			rst_ += val_;
		}
	}

	return rst_;
}
double NonLinearSolver::vec_multiply_vec(vector<double>& vec_1, vector<double>& vec_2)
{
	assert(vec_1.size()==vec_2.size());

	double sum_ = 0;
	for (size_t i = 0; i < vec_1.size(); i++)
	{
		sum_ += vec_1[i]*vec_2[i];
	}
	return sum_;
}
