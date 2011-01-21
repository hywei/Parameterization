#include "CrossParamOptimize.h"
#include "ParamDistortion.h"
#include "TriangleTransFunctor.h"
#include "../Numerical/linear_solver.h"
#include "../hj_3rd/include/zjucad/matrix/io.h"

namespace PARAM
{
	CrossParamHarmonicOptimizer::CrossParamHarmonicOptimizer(CrossParam& _cross_param)
		: CrossParamOptimizer(_cross_param){}

	CrossParamHarmonicOptimizer::~CrossParamHarmonicOptimizer(){}

	void CrossParamHarmonicOptimizer::Optimize()
	{
		printf("Harmonic maps optimizer!\n");	

	    const QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const vector<ParamCoord>& param_coord_array = param_1.GetParamCoord();
		const boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();
		int vert_num = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();
		m_vari_num =0;		
		m_vtx_var_index_mapping.clear();
		m_vtx_var_index_mapping.resize(vert_num, -1);
		
		m_fixed_vtx_flag.clear();
		m_fixed_vtx_flag.resize(vert_num, false);

		for(int vid=0; vid<vert_num; ++vid)
		{
			/// the conners and boundary vertices are fixed
			if( param_1.IsConnerVertex(vid) || p_mesh_1->m_BasicOp.IsBoundaryVertex(vid))    
			{				
				continue;
			}			
			m_vtx_var_index_mapping[vid] = m_vari_num;
		
			m_vari_num ++; 
		}
		m_vari_num *= 2; /// each vertex has two variables (u_i, v_i)

		m_func_value.clear();
		m_func_value.resize(vert_num*2);

		m_variable_value.clear();
		m_variable_value.resize(vert_num*2);
		
		SetInitVariableValue();
		int iter_time = 4;
		for(int k=0; k<iter_time; ++k)
		{
			SolveJacobiLinearEquation();
		}

		ComputeNewtonMethodFuncValue();
		ResignParamCoord();
	}

	void CrossParamHarmonicOptimizer::SolveJacobiLinearEquation()
	{
		SetNewtonMethodJacobiMatrix();
		ComputeNewtonMethodFuncValue();		

		printf("Solve Jacobi Linear Equation!\n");

		QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();
		const PolyIndexArray& vtx_adj_vtx_array = p_mesh_1->m_Kernel.GetVertexInfo().GetAdjVertices();
		int vert_num = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		const CMeshSparseMatrix& lap_coef_matrix = param_1.GetLaplaceMatrix(); 

		LinearSolver linear_solver(m_vari_num );

		linear_solver.begin_equation();

		for(int vid = 0; vid < vert_num; ++vid)
		{
			if(p_mesh_1->m_BasicOp.IsBoundaryVertex(vid)) continue; /// there are no laplance equation in boundary
			int var_index = m_vtx_var_index_mapping[vid];		   
            
//			if(var_index == -1) continue;

			int to_chart_id = param_1.GetVertexChartID(vid);

			for(size_t uv=0; uv<2; ++uv)
			{
				int row_idx = vid*2 + uv;

				const vector<int>& row_index = m_nw_jacobi_mat.m_RowIndex[row_idx];
				const vector<double>& row_data = m_nw_jacobi_mat.m_RowData[row_idx];

				double right_b_(0);

				linear_solver.begin_row();
				for(size_t k=0; k<row_index.size(); ++k)
				{				
					int adj_vtx = row_index[k] / 2;
					int adj_var_index = m_vtx_var_index_mapping[adj_vtx];
					if(adj_var_index == -1) continue;
					int _mod = row_index[k]%2;

					int from_chart_id = param_1.GetVertexChartID(adj_vtx);

					linear_solver.add_coefficient(adj_var_index*2 + _mod, row_data[k]);

// 					if(from_chart_id == to_chart_id)
// 					{
// 						linear_solver.add_coefficient(adj_var_index*2 + _mod, row_data[k]);
// 					}else
// 					{
// 						ChartTransFun trans_func;
// 						TransMode trans_mode = (_mod == 0) ? TRANS_S_MODE : TRANS_T_MODE;
// 
// 						int from_vid = adj_vtx, to_vid = vid;
// 						if(param_1.GetTransFuncBetweenTwoVertics(from_vid, to_vid, trans_mode, trans_func) != 0)
// 						{
// 							printf("Can't find translation function between %d and %d!\n", from_vid, to_vid);
// 						}
// 						if(trans_func.m_mode == S_MODE)
// 						{
// 							adj_var_index = adj_var_index*2 ;
// 						}else
// 						{
// 							adj_var_index = adj_var_index*2 + 1;
// 						}					
// 						linear_solver.add_coefficient(adj_var_index, row_data[k]*trans_func.m_coef);
// 
// 						//right_b_ -= (row_data[k]*trans_func.m_div);
// 					}
				} 
				linear_solver.set_right_hand_side(-m_func_value[row_idx]);
				linear_solver.end_row();
			}
		}

		linear_solver.end_equation();

		linear_solver.solve();	

		for(int vid =0; vid<vert_num; ++vid)
		{
			int var_index = m_vtx_var_index_mapping[vid];
			if( var_index == -1) continue;

			double vari_u_value = linear_solver.variable(var_index*2).value();
			double vari_v_value = linear_solver.variable(var_index*2+1).value();
			m_variable_value[vid*2] += vari_u_value;
			m_variable_value[vid*2+1] += vari_v_value;

		  // printf("%d %.6lf %.6lf %.6lf %.6lf\n", vid, vari_u_value, vari_v_value, m_variable_value[vid*2], m_variable_value[vid*2 + 1]);			
		}

		CheckNewtonMethodResult(linear_solver);
	   
	}

	void CrossParamHarmonicOptimizer::SetInitVariableValue()
	{		
		const QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();
		int vert_num_1 = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		for(int vid=0; vid<vert_num_1; ++vid)
		{
			m_variable_value[vid*2] = param_1.GetParamCoord(vid).s_coord;
			m_variable_value[vid*2+1] = param_1.GetParamCoord(vid).t_coord;
		}
	}

	void CrossParamHarmonicOptimizer::ComputeNewtonMethodFuncValue()
	{
		printf("Compute Function Value!\n");                         
		
		/// compute F(x_i)
		const QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();
		const PolyIndexArray& vtx_adj_vtx_array = p_mesh_1->m_Kernel.GetVertexInfo().GetAdjVertices();
		int vert_num = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();
		
		const CMeshSparseMatrix& lap_coef_matrix = param_1.GetLaplaceMatrix(); 

		const QuadParam& param_2 = m_cross_param.GetQuadParam2();
		const boost::shared_ptr<MeshModel> p_mesh_2 = param_2.GetMeshModel();
		const PolyIndexArray& face_index_array = p_mesh_2->m_Kernel.GetFaceInfo().GetIndex();
		const CoordArray& vtx_coord_array = p_mesh_2->m_Kernel.GetVertexInfo().GetCoord();
		ParamDistortion param_distortion_2(param_2);

		TriTransFunctor tri_trans_functor(p_mesh_2);

		for(int vid=0; vid<vert_num; ++vid)
		{
			//if(m_vtx_var_index_mapping[vid] == -1) continue;
			if(p_mesh_1->m_BasicOp.IsBoundaryVertex(vid)) 
			{
				m_func_value[vid*2] = m_func_value[vid*2 + 1] = 0;
				continue;
			}

			int vtx_chart_id = param_1.GetVertexChartID(vid);
			ParamCoord vtx_param_coord(m_variable_value[vid*2], m_variable_value[vid*2+1]);
			SurfaceCoord vtx_surface_coord;
			if(!param_distortion_2.ComputeSurfaceCoord(ChartParamCoord(vtx_param_coord, vtx_chart_id), vtx_surface_coord))
			{
				printf("@@@@@@Warning: Can't find valid barycentric for vertex %d : (%d %lf %lf)!\n", 
					vid, vtx_chart_id, vtx_param_coord.s_coord, vtx_param_coord.t_coord);
				printf("Surface Coord : (%d %lf %lf %lf)\n", vtx_surface_coord.face_index, vtx_surface_coord.barycentric[0],
					vtx_surface_coord.barycentric[1], vtx_surface_coord.barycentric[2]);

				
			}
			int face_id_surface_2 = vtx_surface_coord.face_index;

			const vector<int>& row_index = lap_coef_matrix.m_RowIndex[vid];
			const vector<double>& row_data = lap_coef_matrix.m_RowData[vid];

			double u_func_value(0.0), v_func_value(0.0);
			for(size_t k=0; k<row_index.size(); ++k)
			{
				int adj_vtx = row_index[k];
				int chart_id = param_1.GetVertexChartID(adj_vtx);
				ParamCoord param_coord(m_variable_value[adj_vtx*2], m_variable_value[adj_vtx*2 + 1]);
// 				if(chart_id != vtx_chart_id)
// 				{
// 					param_1.TransParamCoordBetweenTwoChart(chart_id, vtx_chart_id, param_1.GetParamCoord(adj_vtx), param_coord);
// 				}else
// 				{
// 					param_coord = param_1.GetParamCoord(adj_vtx);
// 				}
				SurfaceCoord surface_coord;
				if(!param_distortion_2.ComputeSurfaceCoord(ChartParamCoord(param_coord, chart_id), surface_coord))
				{
				//	printf("@@@@@@Warning: Can't find valid barycentric for vertex %d!\n", adj_vtx);
				}
				int face_id = surface_coord.face_index;
				Coord by_coord = surface_coord.barycentric;

				vector<Coord> tri_vtx_coord(3);
				const IndexArray& face_index = face_index_array[face_id];
				for(size_t i=0; i<3; ++i) tri_vtx_coord[i] = vtx_coord_array[face_index[i]];
				vector<Coord2D> local_coord_array = ParamDistortion::ComputeLocal2DCoord(tri_vtx_coord);
				Coord2D vtx_local_coord = local_coord_array[0]*by_coord[0] + local_coord_array[1]*by_coord[1] + local_coord_array[2]*by_coord[2];
				Coord2D trans_vtx_local_coord(vtx_local_coord);
				tri_trans_functor.TransLocalCoordBetweenCharts(face_id, vtx_local_coord, face_id_surface_2, trans_vtx_local_coord);

				u_func_value += (row_data[k]) * trans_vtx_local_coord[0];
				v_func_value += (row_data[k]) * trans_vtx_local_coord[1];

			}
			m_func_value[vid*2] = u_func_value;
			m_func_value[vid*2+1] = v_func_value; 

			//printf("Function value : %d %lf %lf\n", vid, u_func_value, v_func_value);
		}

		ComputeNewtonMethodError();
	}

	void CrossParamHarmonicOptimizer::ComputeParamJacobiMatrix()
	{
		const QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();
		int vert_num_1 = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		ParamDistortion para_distortion_2(m_cross_param.GetQuadParam2());
		m_param_jacobi_matrix.clear(); 
		m_param_jacobi_matrix.resize(vert_num_1);

		for(int vid=0; vid<vert_num_1; ++vid)
		{
			int chart_id = param_1.GetVertexChartID(vid);
			ParamCoord param_coord( m_variable_value[vid*2], m_variable_value[vid*2 + 1]);
			SurfaceCoord surface_coord;
			if(!para_distortion_2.ComputeSurfaceCoord(ChartParamCoord(param_coord, chart_id), surface_coord))
			{
				printf("@@@@@@Can't find valid barycentric , vid = %d!@@@@@@\n", vid);
			}

			if(!ComputeVertexJacobiMatrixOnSurface2(surface_coord, m_param_jacobi_matrix[vid]))
			{
				printf("@@@@@@Can't compute vertex %d's Jacobi matrix!@@@@@@\n",vid);
			}
		}
	}

	bool CrossParamHarmonicOptimizer::ComputeVertexJacobiMatrixOnSurface2(const SurfaceCoord& surface_coord, 
		zjucad::matrix::matrix<double>& j_mat)
	{
		boost::shared_ptr<MeshModel> p_mesh = m_cross_param.GetQuadParam2().GetMeshModel();
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const PolyIndexArray& vtx_adj_face_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();

	    ParamDistortion param_distortion(m_cross_param.GetQuadParam2());
 		j_mat = param_distortion.ComputeTriJacobiMatrix(surface_coord.face_index);
 		return true;

// 		/// another algorithm
		j_mat.resize(2, 2); 
		j_mat(0, 0) = j_mat(0, 1) = j_mat(1, 0) = j_mat(1, 1) = 0;

		const IndexArray& fv_array = face_index_array[surface_coord.face_index];

		const Coord& by_coord = surface_coord.barycentric;
		int zero_num=0;
		for(size_t k=0; k<3; ++k)
		{
			if(fabs(by_coord[k]) < ParamDistortion::EPSILON) zero_num ++;
		}
		if(zero_num == 0) /// the position is in a triangle's interior
		{
			j_mat = param_distortion.ComputeLocalJacobiMatrix(surface_coord.face_index);
		}else if(zero_num == 1) /// the position is in a triangle's edge
		{
			int vid1, vid2;
			for(size_t k=0; k<3; ++k)
			{
				if(fabs(by_coord[k]) < ParamDistortion::EPSILON) 
				{
					vid1 = fv_array[(k+1)%3];
					vid2 = fv_array[(k+2)%3];
					break;
				}
			}
		
			const IndexArray& adj_face1 = vtx_adj_face_array[vid1];
			const IndexArray& adj_face2 = vtx_adj_face_array[vid2];

			vector<int> edge_adj_face_array;
			for(size_t i=0; i<adj_face1.size(); ++i)
			{
			   int fid = adj_face1[i];
			   if(find(adj_face2.begin(), adj_face2.end(), fid) != adj_face2.end())
			   {
				   edge_adj_face_array.push_back(fid);
			   }
			}
			assert(edge_adj_face_array.size() == 1 || edge_adj_face_array.size() == 2);
			if(edge_adj_face_array.size() == 1)
			{
				j_mat = param_distortion.ComputeTriJacobiMatrix(edge_adj_face_array[0]);
			}else if(edge_adj_face_array.size() == 2)
			{
				j_mat = (param_distortion.ComputeTriJacobiMatrix(edge_adj_face_array[0]) + 
					param_distortion.ComputeTriJacobiMatrix(edge_adj_face_array[1])) / 2;
			}
		}else if(zero_num == 2)
		{
			int vid(-1);
			for(size_t k= 0; k<3; ++k)
			{
				if(fabs(by_coord[k]) < ParamDistortion::EPSILON) continue;
				vid = fv_array[k];
			}
			assert(vid != -1);
			
			const IndexArray& adj_faces = vtx_adj_face_array[vid];
			for(size_t k=0; k<adj_faces.size(); ++k)
			{
				j_mat += param_distortion.ComputeTriJacobiMatrix(adj_faces[k]);
			}
			j_mat = j_mat/(adj_faces.size());
		}
		return true;
	}

	double CrossParamHarmonicOptimizer::ComputeNewtonMethodError()
	{
		double total_error = 0;
		for(size_t k=0; k<m_func_value.size(); ++k)
		{
			total_error += m_func_value[k]*m_func_value[k];
		}
		total_error = sqrt(total_error);
		printf("Newton Method Total Error is %lf\n", total_error);

		return total_error;
	}
	
	void CrossParamHarmonicOptimizer::ResignParamCoord()
	{
		QuadParam& quad_param = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh = quad_param.GetMeshModel();
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		for(int vid = 0; vid < vert_num; ++vid)
		{
			int var_index = m_vtx_var_index_mapping[vid];
			if(var_index != -1)
			{
				quad_param.SetVtxParamCoord(vid, ParamCoord(m_variable_value[vid*2],
					m_variable_value[vid*2 + 1]));
			}
		}

		quad_param.AdjustVertexGroup();
		quad_param.SetAllChartTexture();
	}

	bool CrossParamHarmonicOptimizer::CheckNewtonMethodResult(const LinearSolver& linear_solver)
	{
		const QuadParam& quad_param = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh = quad_param.GetMeshModel();

		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		double func_value_error=0;

		for(int vid = 0; vid < vert_num; ++vid)
		{
			//if(p_mesh->m_BasicOp.IsBoundaryVertex(vid)) continue;
			if(m_vtx_var_index_mapping[vid] == -1) continue;

			for(int uv=0; uv<2; ++uv)
			{
				int row_idx = vid*2 + uv;

				const vector<int>& row_index = m_nw_jacobi_mat.m_RowIndex[row_idx];
				const vector<double>& row_data = m_nw_jacobi_mat.m_RowData[row_idx];
				
				double func_value = 0;
				for(size_t k=0; k<row_index.size(); ++k)
				{				
					int adj_vtx = row_index[k] / 2;
					int adj_var_index = m_vtx_var_index_mapping[adj_vtx];
					if(adj_var_index == -1) continue;
					int _mod = row_index[k]%2;

					double cur_value = linear_solver.variable(adj_var_index*2 + _mod).value();

					func_value += cur_value*row_data[k];				
				}
				double cur_var_error = fabs(func_value + m_func_value[row_idx]);
			//	printf("current variable error: %lf\n", cur_var_error);

				func_value_error += cur_var_error*cur_var_error;
			}
		}
		func_value_error = sqrt(func_value_error);
		printf("Total function error: %lf\n", func_value_error);

		return true;
	}
	

	bool CrossParamHarmonicOptimizer::CheckNewtonMethodJacobiMatrix()
	{
		double delta_x = 0.0001;

		vector<double> _func_value(m_func_value.size());

		const QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();

		int vert_num = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		for(int vid=0; vid<vert_num; ++vid)
		{
			int var_index = m_vtx_var_index_mapping[vid];
			if(var_index != -1)
			{
				const zjucad::matrix::matrix<double>& jacobi_mat = m_param_jacobi_matrix[vid]; 
				
			}
		}
		return true;
	}

	void CrossParamHarmonicOptimizer::SetNewtonMethodJacobiMatrix()
	{
		printf("Set Newton Method Jacobi Matrix...\n");
		const QuadParam& param_1 = m_cross_param.GetQuadParam1();
		const QuadParam& param_2 = m_cross_param.GetQuadParam2();

		boost::shared_ptr<MeshModel> p_mesh_1 = param_1.GetMeshModel();
		boost::shared_ptr<MeshModel> p_mesh_2 = param_2.GetMeshModel();

		ParamDistortion param_distortion_1(param_1);
		ParamDistortion param_distortion_2(param_2);	
		
		const PolyIndexArray& face_list_2 = p_mesh_2->m_Kernel.GetFaceInfo().GetIndex();
		const CoordArray& vtx_coord_array_2 = p_mesh_2->m_Kernel.GetVertexInfo().GetCoord();

		int vert_num_1 = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		TriTransFunctor tri_trans_functor(p_mesh_2);
		
		const CMeshSparseMatrix& lap_mat = param_1.GetLaplaceMatrix();

		m_nw_jacobi_mat.ClearData();
		m_nw_jacobi_mat.SetRowCol(2*vert_num_1, 2*vert_num_1);

		for(int vid = 0; vid < vert_num_1; ++vid)
		{
			
			const vector<int>& lap_row_index = lap_mat.m_RowIndex[vid];
			const vector<double>& lap_row_data = lap_mat.m_RowData[vid];

			int cur_vtx_chart_id = param_1.GetVertexChartID(vid);
			ParamCoord vtx_param_coord(m_variable_value[vid*2], m_variable_value[vid*2+1]);
			SurfaceCoord vtx_surface_coord;
			if(!param_distortion_2.ComputeSurfaceCoord(ChartParamCoord(vtx_param_coord, cur_vtx_chart_id), vtx_surface_coord))
			{
				printf("@@@@@@Warning: Can't find valid barycentric for vertex %d : (%d %lf %lf)!\n", 
					vid, cur_vtx_chart_id, vtx_param_coord.s_coord, vtx_param_coord.t_coord);
				printf("Surface Coord : (%d %lf %lf %lf)\n", vtx_surface_coord.face_index, vtx_surface_coord.barycentric[0],
					vtx_surface_coord.barycentric[1], vtx_surface_coord.barycentric[2]);

			}

			int std_frame_id = vtx_surface_coord.face_index;
			const IndexArray& cur_face = face_list_2[vtx_surface_coord.face_index];
			int std_chart_id = param_2.GetVertexChartID(cur_face[0]);

			std_chart_id = cur_vtx_chart_id;

			for(size_t k=0; k<lap_row_index.size(); ++k)
			{
				int adj_vtx = lap_row_index[k];
				int chart_id = param_1.GetVertexChartID(adj_vtx);
				ParamCoord param_coord(param_1.GetParamCoord(adj_vtx));				

				if(chart_id != std_chart_id) 
				{
					param_1.TransParamCoordBetweenTwoChart(chart_id, std_chart_id, param_1.GetParamCoord(adj_vtx), param_coord);
				}

				SurfaceCoord surface_coord;
				param_distortion_2.ComputeSurfaceCoord(ChartParamCoord(param_coord, std_chart_id), surface_coord);

				int cur_frame_id = surface_coord.face_index;
				const IndexArray& faces = face_list_2[cur_frame_id];
				vector<Coord> vtx_pos_array_3d(3);
				for(int i=0; i<3; ++i) vtx_pos_array_3d[i] = vtx_coord_array_2[faces[i]];
				vector<Coord2D> vtx_pos_array_2d = ComputeTriVertexLocalCoord(vtx_pos_array_3d);				
				for(int i=0; i<3; ++i)
				{
					Coord2D tmp_vtx_pos_2d = vtx_pos_array_2d[i];
					tri_trans_functor.TransLocalCoordBetweenCharts(cur_frame_id, tmp_vtx_pos_2d, std_frame_id, vtx_pos_array_2d[i]);
				}

				vector<ParamCoord> vtx_pc_array(3);
				for(int i=0; i<3; ++i)
				{
					vtx_pc_array[i] = param_2.GetParamCoord(faces[i]);
					int vtx_chart_id = param_2.GetVertexChartID(faces[i]);
					if(vtx_chart_id != std_chart_id)
					{
						ParamCoord tmp_vtx_pc = vtx_pc_array[i];
						param_2.TransParamCoordBetweenTwoChart(vtx_chart_id, std_chart_id, tmp_vtx_pc, vtx_pc_array[i]);
					}
				}

				zjucad::matrix::matrix<double> t_mat = param_1.GetTransMatrixBetweenTwoCharts(chart_id, cur_vtx_chart_id);			 
				zjucad::matrix::matrix<double> j_mat = param_distortion_2.ComputeTriJacobiMatrix(vtx_pos_array_2d, vtx_pc_array);

				zjucad::matrix::matrix<double> sub_t_mat(2, 2);
				for(int i=0; i<2; ++i) for(int j=0; j<2; ++j) sub_t_mat(i, j) = t_mat(i, j);

				zjucad::matrix::matrix<double> res_mat = j_mat * sub_t_mat;				

				double value1 = lap_row_data[k]*(res_mat(0, 0)), value2 = lap_row_data[k]*(res_mat(0, 1));
				double value3 = lap_row_data[k]*(res_mat(1, 0)), value4 = lap_row_data[k]*(res_mat(1, 1));
				
				m_nw_jacobi_mat.SetElement(vid*2, adj_vtx*2, value1);
				m_nw_jacobi_mat.SetElement(vid*2, adj_vtx*2+1, value2);
				m_nw_jacobi_mat.SetElement(vid*2+1, adj_vtx*2, value3);
				m_nw_jacobi_mat.SetElement(vid*2+1, adj_vtx*2+1, value4);				
				
			}
	   
		}
 	}
}