//
// This class is the implementation of a row-dominant sparse matrix.
// Primarily, the class is designed for sparse mesh structure.
//
//////////////////////////////////////////////////////////////////////

#ifndef LINEAR_SOLVER_2_H
#define LINEAR_SOLVER_2_H

#include "MeshSparseMatrix.h"
#include "solver.h"
#include <hj_3rd/hjlib/sparse/sparse.h>

class LinearSolver: public Solver
{
public:
	LinearSolver(int nb_variables);
    virtual ~LinearSolver();

public:
	// __________________ Construction _____________________
	void begin_equation() ;

	void begin_row() ;
	void set_right_hand_side(double b_) ;
	void add_coefficient(int index_, double a_) ;
	void end_row() ;
	void end_equation() ;

	//
	void solve();
	void solve_from_file();

	//
	void factorize();
	void renew_right_b(vector<double>& right_b_vec);
	void set_equation_div_flag();

	void equations_value(vector<double>& var_val_vec);
	size_t get_equation_size(){return m_equation_vec.size();}

	void is_printf_info(bool is_){m_is_printf_info=is_;}
	void print_to_file(vector<double>& var_val_vec, string filename);
	void write_to_file(string ata_filename, string atb_filename);

private:

	void set_solve_matrix();
	void set_solve_b();
	void update_variables();

	bool is_free(int id)   { return (id < nb_free_variables_) ;  }
	bool is_locked(int id) { return (id >= nb_free_variables_) ; }

	void print_f(vector<double>& xc_);
	void print_equation_value(vector<double>& function_vec);

private:
	// User representation
	int nb_free_variables_ ;
	int nb_locked_variables_ ;

	vector<vector<pair<int, double> > > m_equation_vec;
	vector<pair<int, double> > m_current_equ;
	vector<double> m_right_b_vec;
	vector<double> m_xc_;

private:
	//std::auto_ptr<hj::sparse::solver> m_solver_;
	CMeshSparseMatrix m_solve_matrix_AT_;
	vector<double> m_solve_b_vec;

private:
	bool factorize_state;

private:
	vector<size_t> m_equ_div_flag_vec;
	bool m_is_printf_info;

	vector<pair<int, double> > m_tmp_current_equ;
	vector<vector<pair<int, double> > > m_tmp_equation_vec;
};

#endif
