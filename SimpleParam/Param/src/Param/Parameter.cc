#include "ChartCreator.h"
#include "Parameter.h"
#include "TransFunctor.h"
#include "Barycentric.h"
#include "TriDistortion.h"

#include "../ModelMesh/MeshModel.h"
#include "../Numerical/linear_solver.h"
#include "../Numerical/MeshSparseMatrix.h"
#include <hj_3rd/zjucad/matrix/matrix.h>
#include <hj_3rd/zjucad/matrix/io.h>

#include <iostream>
#include <queue>
#include <set>
#include <limits>
#include <fstream>
#include <cmath>

namespace PARAM
{
	Parameter::Parameter(boost::shared_ptr<MeshModel> _p_mesh) : p_mesh(_p_mesh){}
	Parameter::~Parameter(){}

	bool Parameter::LoadPatchFile(const std::string& file_name)
	{
		p_chart_creator = boost::shared_ptr<ChartCreator> ( new ChartCreator(p_mesh));
		p_chart_creator->LoadPatchFile(file_name);
		if(!p_chart_creator->FormParamCharts())
		{
			std::cout<<"Error: Cannot compute parameteriztion!\n";
			return false;
		}
		return true;
	}

	bool Parameter::ComputeParamCoord()
	{
		if(p_mesh == NULL) return false;
		if(p_chart_creator == NULL)
		{
			std::cout<<"Error : Please load quad file first!\n";
			return false;
		}		
        
		SetInitFaceChartLayout();
		SetInitVertChartLayout();        

// 		TransFunctor tri_functor(p_chart_creator);
// 		zjucad::matrix::matrix<double> tran_mat = tri_functor.GetTransMatrix(0, 1);
// 		std::cout<< tran_mat <<std::endl;
// 
// 		ParamCoord to_param_coord(0, 0);	
// 		
// 		TransParamCoordBetweenCharts(1, 3, ParamCoord(1, 1), to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 
// 		TransParamCoordBetweenCharts(1, 0, ParamCoord(1, 1), to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 		
// 		ParamCoord pre_coord = to_param_coord;
// 		TransParamCoordBetweenCharts(0, 2, pre_coord, to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 
// 		pre_coord = to_param_coord;
// 		TransParamCoordBetweenCharts(2, 3, pre_coord, to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 		
// 
// 		TransParamCoordBetweenCharts(2, 3, ParamCoord(0, 0), to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 
// 		TransParamCoordBetweenCharts(2, 3, ParamCoord(1, 1), to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 
// 		TransParamCoordBetweenCharts(3, 2, ParamCoord(1, 0), to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;
// 
// 		TransParamCoordBetweenCharts(3, 2, ParamCoord(0.633975, 1.36603), to_param_coord);
// 		std::cout << to_param_coord.s_coord <<" " << to_param_coord.t_coord << std::endl;



		CMeshSparseMatrix lap_mat;
		SetLapMatrixCoef(p_mesh, lap_mat);
		//SetLapMatrixCoefWithMeanValueCoord(p_mesh, lap_mat);
		
		int loop_num = 1;
		for(int k=0; k<loop_num; ++k)
		{
			SolveParameter(lap_mat);
			if(k < loop_num-1)
			{
                GetOutRangeVertices(m_out_range_vert_array);
				AdjustPatchBoundary();
			}
		}	   		

		AdjustPatchBoundary();
        GetOutRangeVertices(m_out_range_vert_array);
        
		ResetFaceChartLayout();
		SetMeshFaceTextureCoord();

		//GetOutRangeVertices(m_out_range_vert_array);
		SetChartVerticesArray();

		ComputeDistortion();

		CheckFlipedTriangle();

		return true;
	}

	void Parameter::SetInitFaceChartLayout()
	{        		
        int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();
        m_face_chart_array.clear(); m_face_chart_array.resize(face_num);

		const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
		for(size_t k=0; k<patch_array.size(); ++k)
		{
			const ParamPatch& param_patch = patch_array[k];
			const std::vector<int>& faces_in_patch = param_patch.m_face_index_array;
			for(size_t i=0; i<faces_in_patch.size(); ++i)
			{
				int fid = faces_in_patch[i];
				m_face_chart_array[fid] = k;
			}
		}

        ofstream fout("temp.txt");
        for(size_t k=0; k<(size_t)face_num; ++k){
            fout << m_face_chart_array[k] <<" ";
        }
        fout << std::endl;
        fout.close();
	}

	void Parameter::SetInitVertChartLayout()
	{
        
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		m_vert_chart_array.clear(); m_vert_chart_array.resize(vert_num);

		const PolyIndexArray& adj_face_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();

		for(int vid = 0; vid < vert_num; ++vid)
		{
			const IndexArray& adj_faces = adj_face_array[vid];
			std::map<int, int> chart_count_num;
			for(size_t i=0; i<adj_faces.size(); ++i)
			{
				int face_chart_id = m_face_chart_array[adj_faces[i]];
				chart_count_num[face_chart_id] ++;
			}
			int chart_id(-1), max_num(-1);
			for(std::map<int, int>::const_iterator im = chart_count_num.begin(); im != chart_count_num.end(); ++im)
			{
				if(im->second > max_num) { max_num = im->second; chart_id = im->first; }
			}
			m_vert_chart_array[vid] = chart_id;
		}

        ofstream fout("temp.txt");
        for(int k=0; k<vert_num; ++k){
            fout << m_vert_chart_array[k] << ' ';
        }
        fout << std::endl;
        fout.close();
	}
   
	void Parameter::SetBoundaryVertexParamValue(LinearSolver* p_linear_solver /* = NULL */)
	{
		const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
		const std::vector<ParamChart>& chart_array = p_chart_creator->GetChartArray();
		const std::vector<PatchEdge>& patch_edge_array = p_chart_creator->GetPatchEdgeArray();		
		const std::vector<PatchConner>& patch_conner_array = p_chart_creator->GetPatchConnerArray();

		for(size_t k=0; k<patch_array.size(); ++k)
		{
			int chart_id = k;
			const ParamPatch& param_patch = patch_array[k];
			const ParamChart& param_chart = chart_array[k];

			const std::vector<int>& patch_conners = param_patch.m_conner_index_array;
			const std::vector<int>& patch_edges = param_patch.m_edge_index_array;

			/// fix conner vertex
			for(size_t i=0; i<patch_conners.size(); ++i)
			{
				int conner_idx = patch_conners[i];
                const PatchConner& conner = patch_conner_array[conner_idx];
				int conner_vid = conner.m_mesh_index;
				if(m_vert_chart_array[conner_vid] == chart_id) 
				{
					m_vert_param_coord_array[conner_vid].s_coord =
						param_chart.m_conner_param_coord_array[i].s_coord;
					m_vert_param_coord_array[conner_vid].t_coord =
						param_chart.m_conner_param_coord_array[i].t_coord;	

					if(p_linear_solver){
						int var_index = conner_vid*2;
						p_linear_solver->variable(var_index).lock();
						p_linear_solver->variable(var_index).set_value(param_chart.m_conner_param_coord_array[i].s_coord);
						p_linear_solver->variable(var_index + 1).lock();
						p_linear_solver->variable(var_index + 1).set_value(param_chart.m_conner_param_coord_array[i].t_coord);
					}
				}
			}

			/// fix boundary vertex
			for(size_t i=0; i<patch_edges.size(); ++i)
			{
				int edge_idx = patch_edges[i];
				const PatchEdge& patch_edge = patch_edge_array[edge_idx];

				if(patch_edge.m_nb_patch_index_array.size() == 1) /// this is a boundary patch edge
				{
					std::pair<int, int> conner_pair = patch_edge.m_conner_pair_index;
					int conner_vid1 = patch_conner_array[conner_pair.first].m_mesh_index;
					int conner_vid2 = patch_conner_array[conner_pair.second].m_mesh_index;

					int conner_idx_1 = GetConnerIndexInPatch(conner_vid1, k);
					int conner_idx_2 = GetConnerIndexInPatch(conner_vid2, k);
                    
					double start_s_coord = param_chart.m_conner_param_coord_array[conner_idx_1].s_coord;
					double start_t_coord = param_chart.m_conner_param_coord_array[conner_idx_1].t_coord;
					double end_s_coord = param_chart.m_conner_param_coord_array[conner_idx_2].s_coord;
					double end_t_coord = param_chart.m_conner_param_coord_array[conner_idx_2].t_coord;

					const std::vector<int>& mesh_path = patch_edge.m_mesh_path;
					double path_len = ComputeMeshPathLength(mesh_path, 0, mesh_path.size());

					/// arc length parameterization
					for(size_t j=1; j<mesh_path.size()-1; ++j)
					{
						int mesh_vert = mesh_path[j];
						if(m_vert_chart_array[mesh_vert] != chart_id)
						{
							std::cout<<"Adjust boundary vertex" << mesh_vert <<" from chart "
								<< m_vert_chart_array[mesh_vert] << " to " << chart_id << std::endl;
							m_vert_chart_array[mesh_vert] = chart_id;
						}

						double arc_len = ComputeMeshPathLength(mesh_path, 0, j+1);
						double lambda = arc_len / path_len;
						double s_coord = (1-lambda)*start_s_coord + lambda*end_s_coord;
						double t_coord = (1-lambda)*start_t_coord + lambda*end_t_coord;

						m_vert_param_coord_array[mesh_vert].s_coord = s_coord;						
						m_vert_param_coord_array[mesh_vert].t_coord = t_coord;

						if(p_linear_solver){
							int var_index = mesh_vert*2;
							p_linear_solver->variable(var_index).lock();
							p_linear_solver->variable(var_index).set_value(s_coord);
							p_linear_solver->variable(var_index + 1).lock();
							p_linear_solver->variable(var_index + 1).set_value(t_coord);
						}
					}
				}
			}

		}
	}

	void Parameter::SolveParameter(const CMeshSparseMatrix& lap_mat)
	{
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

        m_vert_param_coord_array.clear();
		m_vert_param_coord_array.resize(vert_num);
        
		vector<int> vari_index_mapping;
		SetVariIndexMapping(vari_index_mapping);

        // for(size_t k=0; k<m_vert_param_coord_array.size(); ++k){
        //     if(vari_index_mapping[k] == -1)
        //         std::cout << m_vert_chart_array[k] <<" " <<  m_vert_param_coord_array[k].s_coord << ' ' <<
        //             m_vert_param_coord_array[k].t_coord << std::endl;
        // }

		int vari_num = (int)vari_index_mapping.size()*2;
				
		std::cout << "Begin solve parameterization: variable num "<< vari_num << std::endl;

		LinearSolver linear_solver(vari_num);

		SetBoundaryVertexParamValue(&linear_solver);        
		
		linear_solver.begin_equation();
		for(int vid = 0; vid < vert_num; ++vid)
		{		
			/// there are no laplance equation on boundary vertex						
			if(vari_index_mapping[vid] == -1) continue;

		    int to_chart_id = m_vert_chart_array[vid];
			
			const std::vector<int>& row_index = lap_mat.m_RowIndex[vid];
			const std::vector<double>& row_data = lap_mat.m_RowData[vid];

			for(int st=0; st<2; ++st)
			{
				linear_solver.begin_row();

				bool have_boundary_vert(false);
				for(size_t k=0; k<row_index.size(); ++k)
				{
					int col_vert = row_index[k];
					int val_index = vari_index_mapping[col_vert];
					if(val_index == -1) 
					{
						have_boundary_vert = true;
						break;
					}
				}

				double right_b = 0;
				for(size_t k=0; k<row_index.size(); ++k)
				{
					int col_vert = row_index[k];
					int from_chart_id = m_vert_chart_array[col_vert];
					int var_index = vari_index_mapping[col_vert];
                    
					double lap_weight = row_data[k];
									   
					if(from_chart_id == to_chart_id)
					{
						if( var_index == -1)
						{
							double st_value = (st==0) ? m_vert_param_coord_array[col_vert].s_coord : 
								m_vert_param_coord_array[col_vert].t_coord;
							
							right_b -= lap_weight*st_value;
						}else
						{
							linear_solver.add_coefficient(var_index*2 + st, lap_weight);
						}
					}else
					{
						if(var_index == -1) 
						{
							ParamCoord param_coord;
							TransParamCoordBetweenCharts(from_chart_id, to_chart_id, col_vert, 
								m_vert_param_coord_array[col_vert], param_coord);
							double st_value = (st==0) ? param_coord.s_coord : param_coord.t_coord;
							right_b -= lap_weight*st_value;
						}else
						{
                            
                            zjucad::matrix::matrix<double> trans_mat = GetTransMatrix(col_vert, vid, from_chart_id, to_chart_id);
							double a = (st == 0) ? trans_mat(0, 0) : trans_mat(1, 0);
							double b = (st == 0) ? trans_mat(0, 1) : trans_mat(1, 1);
							double c = (st == 0) ? trans_mat(0, 2) : trans_mat(1, 2);
                            
							///! a*u + b*v + c
							if(!(fabs(a) < LARGE_ZERO_EPSILON))
							{
								linear_solver.add_coefficient(var_index*2, lap_weight*a);
							}
							if(!(fabs(b) < LARGE_ZERO_EPSILON))
							{
								linear_solver.add_coefficient(var_index*2+1, lap_weight*b);
							}
							right_b -= c*lap_weight;
						}
					}					
				}
				linear_solver.set_right_hand_side(right_b);
				linear_solver.end_row();

			}
		}
		linear_solver.end_equation();

		linear_solver.solve();

		for(int vid=0; vid < vert_num; ++vid)
		{
			int vari_index = vari_index_mapping[vid];
			if(vari_index != -1)
			{
				m_vert_param_coord_array[vid].s_coord = linear_solver.variable(vari_index*2).value();
				m_vert_param_coord_array[vid].t_coord = linear_solver.variable(vari_index*2+1).value();
			}
		}

        ofstream fout ("parame.txt");
        for(size_t k=0; k<m_vert_param_coord_array.size(); ++k){
            fout <<  m_vert_param_coord_array[k].s_coord << ' ' <<
                m_vert_param_coord_array[k].t_coord << std::endl;
        }
        fout.close();

	}

	void Parameter::SetVariIndexMapping(std::vector<int>& vari_index_mapping)
	{
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		vari_index_mapping.clear();
		vari_index_mapping.resize(vert_num, -1);

        std::vector<int> conner_vertices;
        const std::vector<PatchConner>& patch_conners = p_chart_creator->GetPatchConnerArray();

        for(size_t k=0; k<patch_conners.size(); ++k){
            int vid = patch_conners[k].m_mesh_index;
            conner_vertices.push_back(vid);
        }
            
        
		int vari_num = 0;
		for(int vid=0; vid < vert_num; ++vid){
			if(p_mesh->m_BasicOp.IsBoundaryVertex(vid) ){
//				|| find(conner_vertices.begin(), conner_vertices.end(), vid) != conner_vertices.end()){
				continue;
			}
			vari_index_mapping[vid] = vari_num++;
		}
	}

	void Parameter::AdjustPatchBoundary()
	{
		const std::vector<ParamChart>& param_chart_array = p_chart_creator->GetChartArray();
		const PolyIndexArray& vert_adjvertices_array =p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();

		bool swap_able = false;
		do{
			swap_able = false;
			int adjust_num(0);
			for(int vid = 0; vid < (int)vert_adjvertices_array.size(); ++vid)
			{
				int chart_id = m_vert_chart_array[vid];
				ParamCoord param_coord = m_vert_param_coord_array[vid];
				const ParamChart& param_chart = param_chart_array[chart_id];
				if(param_chart.InValidRangle(param_coord)) continue;

				bool flag = false;
				double out_range_error = ComputeOutRangeError4EqualTriangle(param_coord);
				double min_out_range_error = out_range_error;
				int min_error_chart_id = chart_id;
				ParamCoord min_error_param_coord;

				const IndexArray& adj_vertices = vert_adjvertices_array[vid];
				for(size_t k=0; k<adj_vertices.size(); ++k)
				{
					int adj_chart_id = m_vert_chart_array[adj_vertices[k]];
					if(chart_id == adj_chart_id) continue;

					ParamCoord adj_param_coord = m_vert_param_coord_array[adj_vertices[k]];
					const ParamChart& adj_chart = param_chart_array[adj_chart_id];
					if(!adj_chart.InValidRangle(adj_param_coord)) continue;

					if(chart_id != adj_chart_id)
					{
// 						TransParamCoordBetweenCharts(chart_id, adj_chart_id, 
// 							m_vert_param_coord_array[vid], param_coord);
						TransParamCoordBetweenCharts(chart_id, adj_chart_id, vid, 
							m_vert_param_coord_array[vid], param_coord);
					}

					if(adj_chart.InValidRangle(param_coord))
					{
						m_vert_chart_array[vid] = adj_chart_id;
						m_vert_param_coord_array[vid] = param_coord;
						flag = true; 
						swap_able = true;
						break;
					}else
					{
						double cur_out_range_error = ComputeOutRangeError4EqualTriangle(param_coord);
						if(cur_out_range_error < min_out_range_error)
						{
							min_out_range_error = cur_out_range_error;
							min_error_chart_id = adj_chart_id;
							min_error_param_coord = param_coord;
						}
					}

				}
				if(flag == false)
				{
					// 	if(chart_id == min_error_chart_id) 
					// 	 std::cout<<"Can't adjust this vertex " << vid << std::endl;
					// 				m_vert_chart_array[vid] = min_error_chart_id;
					// 				m_vert_param_coord_array[vid] = min_error_param_coord;
				}else
				{
					adjust_num ++;
				}
			}
			std::cout << "Adjust " << adjust_num << "vertices." << std::endl;
		}while(swap_able);
		
	}

	void Parameter::GetOutRangeVertices(std::vector<int>& out_range_vert_array) const
	{
		const std::vector<ParamChart>& chart_array = p_chart_creator->GetChartArray();

		out_range_vert_array.clear();
		for(size_t vid=0; vid<m_vert_param_coord_array.size(); ++vid)
		{
			int vert_chart_id = m_vert_chart_array[vid];
			const ParamChart& param_chart = chart_array[vert_chart_id];
			if(!param_chart.InValidRangle(m_vert_param_coord_array[vid]))
			{
				out_range_vert_array.push_back(vid);
			}
		}

		std::cout<<"There are " << out_range_vert_array.size() << " out range vertices.\n";
	}

	bool Parameter::FindValidChartForOutRangeVertex(int out_range_vert, int max_ringe_num /* = 5 */)
	{
		int init_chart_id = m_vert_chart_array[out_range_vert];
		ParamCoord init_param_coord = m_vert_param_coord_array[out_range_vert];

		const std::vector<ParamChart>& param_chart_array = p_chart_creator->GetChartArray();
		const std::vector<ParamPatch>& param_patch_array = p_chart_creator->GetPatchArray();

		int valid_chart_id(-1);
		ParamCoord valid_param_coord;

		std::set<int> visited_chart;

		int steps=0;
		std::queue<int> q;
		q.push(init_chart_id); 
		q.push(steps);
		visited_chart.insert(init_chart_id);

		double min_out_range_error = numeric_limits<double>::infinity();
		int min_error_valid_chart_id(-1);
		ParamCoord min_error_valid_param_coord;

		while(!q.empty())
		{
			int cur_chart_id = q.front(); q.pop();
			steps = q.front(); q.pop();

			ParamCoord cur_param_coord;
			TransParamCoordBetweenCharts(init_chart_id, cur_chart_id, 
				out_range_vert,init_param_coord, cur_param_coord);
			const ParamChart& cur_param_chart = param_chart_array[cur_chart_id];
			if(cur_param_chart.InValidRangle(cur_param_coord))
			{
				valid_chart_id = cur_chart_id; 
				valid_param_coord = cur_param_coord;
				break;
			}else
			{
				double cur_out_range_error = ComputeOutRangeError4EqualTriangle(cur_param_coord);
				if(cur_out_range_error < min_out_range_error)
				{
					min_out_range_error = cur_out_range_error;
					min_error_valid_param_coord = cur_param_coord;
					min_error_valid_chart_id = cur_chart_id;
				}
			}

			if(steps > max_ringe_num) break;

			const ParamPatch& cur_param_patch = param_patch_array[cur_chart_id];
			const std::vector<int>& chart_neighbor = cur_param_patch.m_nb_patch_index_array;

			for(size_t k=0; k<chart_neighbor.size(); ++k)
			{
				int nb_chart_id = chart_neighbor[k];
				if(visited_chart.find(nb_chart_id) == visited_chart.end())
				{
					q.push(nb_chart_id);
					q.push(steps+1);
					visited_chart.insert(nb_chart_id);
				}
			}

		}

		if(valid_chart_id != -1)
		{
			m_vert_chart_array[out_range_vert] = valid_chart_id;
			m_vert_param_coord_array[out_range_vert] = valid_param_coord;
			return true;
		}else
		{
			//: TODO: we should find a min-unvalid-error chart as this vertex's chart
			m_vert_chart_array[out_range_vert] = min_error_valid_chart_id;
			m_vert_param_coord_array[out_range_vert] = min_error_valid_param_coord;
			return false;
		}
		return false;
	}

	void Parameter::TransParamCoordBetweenCharts(int from_chart_id, int to_chart_id, 
		const ParamCoord& from_param_coord, ParamCoord& to_param_coord) const
	{	
		TransFunctor tran_functor(p_chart_creator);
		zjucad::matrix::matrix<double> tran_mat = 
			tran_functor.GetTransMatrix(from_chart_id, to_chart_id);

		to_param_coord.s_coord = tran_mat(0, 0)*from_param_coord.s_coord + 
			tran_mat(0, 1)*from_param_coord.t_coord + tran_mat(0, 2);
		to_param_coord.t_coord = tran_mat(1, 0)*from_param_coord.s_coord +
			tran_mat(1, 1)*from_param_coord.t_coord + tran_mat(1, 2);
	}

	void Parameter::TransParamCoordBetweenCharts(int from_chart_id, int to_chart_id, int vid, 
		const ParamCoord& from_param_coord, ParamCoord& to_param_coord) const
	{
		TransFunctor tran_functor(p_chart_creator);
		
		bool is_ambiguity = IsAmbiguityChartPair(from_chart_id, to_chart_id);
		
		if(is_ambiguity) {
			tran_functor.TransParamCoordBetweenAmbiguityCharts(
			vid, from_chart_id, to_chart_id, from_param_coord, to_param_coord);
		}else{
			zjucad::matrix::matrix<double> tran_mat = 
				tran_functor.GetTransMatrix(from_chart_id, to_chart_id);

			to_param_coord.s_coord = tran_mat(0, 0)*from_param_coord.s_coord + 
				tran_mat(0, 1)*from_param_coord.t_coord + tran_mat(0, 2);
			to_param_coord.t_coord = tran_mat(1, 0)*from_param_coord.s_coord +
				tran_mat(1, 1)*from_param_coord.t_coord + tran_mat(1, 2);
		}
	}

	zjucad::matrix::matrix<double> Parameter::GetTransMatrix(int from_vid, int to_vid, int from_chart_id, int to_chart_id) const
	{
		TransFunctor tran_functor(p_chart_creator);
		bool is_ambiguity = IsAmbiguityChartPair(from_chart_id, to_chart_id);
		if(is_ambiguity){
			return tran_functor.GetTransMatrixBetweenAmbiguityCharts(from_vid, from_chart_id, to_vid, to_chart_id);
		}else{
			return tran_functor.GetTransMatrix(from_chart_id, to_chart_id);
		}
	}


	double Parameter::ComputeMeshPathLength(const std::vector<int>& mesh_path, int start_idx, int end_idx) const
	{
		assert(start_idx < end_idx && start_idx >=0 && start_idx < (int)mesh_path.size() -1&&
			end_idx >0 && end_idx <=(int) mesh_path.size());

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		double len = 0.0;
		for(int k=start_idx+1; k!=end_idx; ++k)
		{
			int vtx1 = mesh_path[k-1];
			int vtx2 = mesh_path[k];
			len += (vCoord[vtx2] - vCoord[vtx1]).abs(); 
		}
		return len;
	}

	int Parameter::GetConnerIndexInPatch(int conner_id, int patch_id)
	{
		const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
		const std::vector<PatchConner>& conner_array = p_chart_creator->GetPatchConnerArray();

		assert(patch_id < (int)patch_array.size());
		const ParamPatch& param_patch = patch_array[patch_id];
		const std::vector<int>& patch_conner_array = param_patch.m_conner_index_array;

		for(size_t k=0; k<patch_conner_array.size(); ++k)
		{
			int conner_idx = patch_conner_array[k];
			if(conner_array[conner_idx].m_mesh_index == conner_id) return (int)k;
		}

		return -1;
	}

	void Parameter::ResetFaceChartLayout()
	{
		const std::vector<ParamChart>& param_chart_array = p_chart_creator->GetChartArray();

		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		size_t face_num = face_list_array.size();
		for (size_t i = 0; i < face_num; ++i)
		{
			const IndexArray& faces= face_list_array[i];
			std::vector<int> chart_id_vec;
			for (size_t j = 0; j < 3; j++)
			{
				chart_id_vec.push_back(m_vert_chart_array[faces[j]]);

			}

			int std_chart_id(-1);
			int min_error_chart_id(-1);
			double min_out_range_error = numeric_limits<double>::infinity();
			for(size_t j=0; j<3; ++j)
			{
				int chart_id = m_vert_chart_array[faces[j]];
				bool flag = true;
				double out_range_error = 0;
				for(size_t k=0; k<3; ++k)
				{
					int cur_chart_id = m_vert_chart_array[faces[k]];
					ParamCoord cur_param_coord = m_vert_param_coord_array[faces[k]];
					if(cur_chart_id != chart_id)
					{
						TransParamCoordBetweenCharts(cur_chart_id,chart_id, faces[k],
							m_vert_param_coord_array[faces[k]], cur_param_coord);
					}
					const ParamChart& param_chart = param_chart_array[cur_chart_id];
					if(!param_chart.InValidRangle(cur_param_coord))
					{
						flag = false;		
						out_range_error += ComputeOutRangeError4EqualTriangle(cur_param_coord);
					}
				}
				if(flag == true)
				{
					std_chart_id = chart_id;
					break;
				}else{
					if(min_out_range_error > out_range_error) 
					{
						min_out_range_error = out_range_error;
						min_error_chart_id = chart_id;
					}
				}
			}
			if(std_chart_id == -1)
			{
				m_unset_layout_face_array.push_back(i);
				if(min_error_chart_id != -1) std_chart_id = min_error_chart_id;
				else{
				//std::cout<<"Can't find valid chart for this face's three vertices!\n";
					int max_times(-1);
					for(size_t k=0; k<chart_id_vec.size(); ++k){
						int c = chart_id_vec[k];
						int cur_times = count(chart_id_vec.begin(), chart_id_vec.end(), c);
						if(cur_times > max_times)
						{
							max_times = cur_times;
							std_chart_id = c;
						}
					}
				}
			}

			m_face_chart_array[i] = std_chart_id;
		}

		ofstream fout("unset_face.txt");
		for(size_t k=0; k<m_unset_layout_face_array.size(); ++k)
		{
			fout << m_unset_layout_face_array[k] << " ";
		}
		fout << std::endl;

		int colors[48][3] = 
		{
			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};

		
		ColorArray& colorArray = p_mesh->m_Kernel.GetFaceInfo().GetColor();
		colorArray.clear();
		colorArray.resize(face_num);

		for(size_t i = 0; i < m_face_chart_array.size(); ++ i)
		{
			int index = m_face_chart_array[i] % 48;
			colorArray[i] = Color(colors[index][0], colors[index][1], colors[index][2]);
		}

	}

	void Parameter::SetMeshFaceTextureCoord()
	{
		PolyTexCoordArray& face_texcoord_array = p_mesh->m_Kernel.GetFaceInfo().GetTexCoord();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		size_t face_num = face_list_array.size();
		face_texcoord_array.clear(); face_texcoord_array.resize(face_num);

		for(size_t fid = 0; fid < face_list_array.size(); ++fid)
		{
			TexCoordArray& face_tex = face_texcoord_array[fid];
			face_tex.clear(); face_tex.resize(3);

			const IndexArray& faces = face_list_array[fid];			
			int face_chart_id = m_face_chart_array[fid];
            //		face_chart_id = FindBestChartIDForTriShape(fid);

			for(int k=0; k<3; ++k)
			{
				int vid = faces[k];
				int vert_chart_id = m_vert_chart_array[vid];
				ParamCoord vert_param_coord = m_vert_param_coord_array[vid];
				if(vert_chart_id != face_chart_id)
				{
					TransParamCoordBetweenCharts(vert_chart_id, face_chart_id, vid,
						m_vert_param_coord_array[vid], vert_param_coord);
				}
				face_tex[k] = TexCoord(vert_param_coord.s_coord, vert_param_coord.t_coord);
			}
		}

	}

	void Parameter::SetChartVerticesArray()
	{
		int colors[16][3] = 
		{
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};

		int chart_num = p_chart_creator->GetChartNumber();
		m_chart_vertices_array.clear();
		m_chart_vertices_array.resize(chart_num);

		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		ColorArray& vtx_color_array = p_mesh->m_Kernel.GetVertexInfo().GetColor();
		vtx_color_array.clear(); vtx_color_array.resize(vert_num);

		for(int vid = 0; vid<vert_num; ++vid)
		{
			int chart_id = m_vert_chart_array[vid];
			m_chart_vertices_array[chart_id].push_back(vid);

			vtx_color_array[vid] = Color(colors[chart_id%16][0], colors[chart_id%16][1], colors[chart_id%16][2]);
		}
		
	}

    double Parameter::ComputeOutRangeError4Square(ParamCoord param_coord) const
    {
        double error=0;
        if(GreaterEqual(param_coord.s_coord, 1) || LessEqual(param_coord.s_coord, 0)){
            if(LessEqual(param_coord.s_coord, 0)) error += param_coord.s_coord*param_coord.s_coord;
            else if(GreaterEqual(param_coord.s_coord, 1)) error += (param_coord.s_coord-1)*(param_coord.s_coord-1);
        }if(GreaterEqual(param_coord.t_coord, 1) || LessEqual(param_coord.t_coord, 0))	{
            if(LessEqual(param_coord.t_coord, 0)) error += param_coord.t_coord*param_coord.t_coord;
            else if(GreaterEqual(param_coord.t_coord, 1)) error += (param_coord.t_coord-1)*(param_coord.t_coord-1);
        }
        return sqrt(error);
    }

	double Parameter::ComputeOutRangeError4EqualTriangle(ParamCoord param_coord) const
	{
		/// you must be sure that this parameter coordinate is outrange
		double sqrt_3 = sqrt(3.0)/2;
		double dist1 = dis_ptoline(0, 0, 0.5, sqrt_3, param_coord.s_coord, param_coord.t_coord);
		double dist2 = dis_ptoline(0, 0, 1.0, 0, param_coord.s_coord, param_coord.t_coord);
		double dist3 = dis_ptoline(1.0, 0, 0.5, sqrt_3, param_coord.s_coord, param_coord.t_coord);

		double min_dist = min(dist1, min(dist2, dist3));
		return min_dist;
	}

	int Parameter::FindBestChartIDForTriShape(int fid) const
	{
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& faces = face_list_array[fid];

		int c_0 = m_vert_chart_array[faces[0]];
		int c_1 = m_vert_chart_array[faces[1]];
		int c_2 = m_vert_chart_array[faces[2]];

		if(c_0 == c_1 && c_0 == c_2) return c_0;

		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& vtx_adjacent_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
		std::set<int> candidate_chart_set;

		for(int k=0; k<3; ++k)
		{
			int vid = faces[k];
			const IndexArray& adj_vertices = vtx_adjacent_array[vid];
			for(size_t i=0; i<adj_vertices.size(); ++i)
			{
				int chart_id = m_vert_chart_array[adj_vertices[i]];
				candidate_chart_set.insert(chart_id);
			}
		}

		double min_angle_error = numeric_limits<double>::infinity();
		int min_error_chart_id = -1;
		for(std::set<int>::const_iterator is = candidate_chart_set.begin(); is!=candidate_chart_set.end(); ++is)
		{
			int chart_id = *is;
			std::vector<Coord2D> vtx_param_coord_array(3);
			for(int k=0; k<3; ++k)
			{
				int cur_vid = faces[k];
				int cur_chart_id = m_vert_chart_array[cur_vid];
				ParamCoord cur_param_coord = m_vert_param_coord_array[cur_vid];
				if(cur_chart_id != chart_id)
				{
					TransParamCoordBetweenCharts(cur_chart_id, chart_id, cur_vid, 
						m_vert_param_coord_array[cur_vid], cur_param_coord);
				}
				vtx_param_coord_array[k] =Coord2D(cur_param_coord.s_coord, cur_param_coord.t_coord);
			}
			double angle_error = 0;
			for(int k=0; k<3; ++k)
			{
				int vid1 = faces[k];
				int vid2 = faces[(k+1)%3];
				int vid3 = faces[(k+2)%3];
				Coord vec_1 = vtx_coord_array[vid2] - vtx_coord_array[vid1];
				Coord vec_2 = vtx_coord_array[vid3] - vtx_coord_array[vid1];
				double angle_3d = angle(vec_1, vec_2);

				Coord2D vec_3 = vtx_param_coord_array[(k+1)%3] - vtx_param_coord_array[k];
				Coord2D vec_4 = vtx_param_coord_array[(k+2)%3] - vtx_param_coord_array[k];
				double angle_2d = angle(vec_3, vec_4);

				angle_error += (angle_3d - angle_2d) * (angle_3d - angle_2d);
			}
			angle_error = sqrt(angle_error);

			if(angle_error < min_angle_error)
			{
				min_angle_error = angle_error;
				min_error_chart_id = chart_id;
			}
		}

		return min_error_chart_id;
	}

	bool Parameter::FindCorrespondingOnSurface(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const
	{
		int chart_id = chart_param_coord.chart_id;
		ParamCoord param_coord = chart_param_coord.param_coord;
		
		if(!FindCorrespondingInChart(chart_param_coord, chart_id, surface_coord))
		{
			std::cout<<"Cann't find corresponding vertex" << std::endl;			
			return false;
		}

		return true;
	}

	bool Parameter::FindCorrespondingInChart(const ChartParamCoord& chart_param_coord, 
		int chart_id, SurfaceCoord& surface_coord) const
	{
		int origin_chart_id =chart_param_coord.chart_id;
		ParamCoord origin_param_coord = chart_param_coord.param_coord;
		if(origin_chart_id != chart_id) 
		{
			TransParamCoordBetweenCharts(origin_chart_id, chart_id, chart_param_coord.param_coord, origin_param_coord);
		}

		const PolyIndexArray& vert_adj_face_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		const std::vector<int>& vertics_array = m_chart_vertices_array[chart_id];
		
		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();
		std::vector<bool> face_visited_flag(face_num, false);

		double min_out_bary_error = numeric_limits<double>::infinity();
	
		for(size_t k=0; k<vertics_array.size(); ++k)
		{
			int vid = vertics_array[k];		
			const IndexArray& adj_face_array = vert_adj_face_array[vid];
			for(size_t i=0; i<adj_face_array.size(); ++i)
			{
				int fid = adj_face_array[i];
				if(!face_visited_flag[fid])
				{
					std::vector<ParamCoord> vert_param_coord(3);
					const IndexArray& faces = face_list_array[fid];
					for(size_t j=0; j<3; ++j)
					{
						vert_param_coord[j] = m_vert_param_coord_array[faces[j]];
						if(m_vert_chart_array[faces[j]] != chart_id)
						{
							TransParamCoordBetweenCharts(m_vert_chart_array[faces[j]], chart_id,
								m_vert_param_coord_array[faces[j]], vert_param_coord[j]);
						}
					}
					Barycentrc baryc = ComputeVertexBarycentric(vert_param_coord, origin_param_coord);
					if(IsValidBarycentic(baryc))
					{
						surface_coord = SurfaceCoord(fid, baryc);
						return true;
					}else
					{
						double out_by_error = ComputeErrorOutValidBarycentric(baryc);
						if(out_by_error < min_out_bary_error)
						{
							min_out_bary_error = out_by_error;
							surface_coord = SurfaceCoord(fid, baryc);
						}
					}
				}
			}
		}
		return false;
	}

	void Parameter::ComputeDistortion()
	{

		TriDistortion tri_distortion(*this);
		tri_distortion.ComputeDistortion();

		const std::vector<double>& face_harmonic_distortion = tri_distortion.GetFaceHarmonicDistortion();
		const std::vector<double>& face_isometric_distortion = tri_distortion.GetFaceIsometricDistortion();

		FaceValue2VtxColor(p_mesh, face_harmonic_distortion);
		
	}

	void Parameter::CheckFlipedTriangle()
	{
		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();

		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		m_fliped_face_array.clear();

		for(int fid =0; fid < face_num; ++fid)
		{
			const IndexArray& faces = face_list_array[fid];
			int face_chart_id = m_face_chart_array[fid];
			std::vector<double> u(3), v(3);
			for(int i=0; i<3; ++i)
			{
				int vid = faces[i];
				int chart_id = m_vert_chart_array[vid];
				ParamCoord param_coord = m_vert_param_coord_array[vid];
				if(chart_id != face_chart_id)
				{
					TransParamCoordBetweenCharts(chart_id, face_chart_id, vid,
						m_vert_param_coord_array[vid], param_coord);
				}
				u[i] = param_coord.s_coord;
				v[i] = param_coord.t_coord;
			}
			/// (v2-v0)*(u1-u0) - (u2-u0)*(v1-v0)
			double area = (u[1] - u[0]) * (v[2] - v[0]) - (u[2] - u[0]) * (v[1]-v[0]);
			if(fabs(area) < LARGE_ZERO_EPSILON) continue;
			if(area < 0 ) m_fliped_face_array.push_back(fid);
		}
		std::cout << "There are " << m_fliped_face_array.size() <<" fliped faces." << std::endl;
	}

	bool Parameter::IsAmbiguityChartPair(int chart_id_1, int chart_id_2) const
	{
		const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
		const ParamPatch& patch_1 = patch_array[chart_id_1];
		const ParamPatch& patch_2 = patch_array[chart_id_2];

		const std::vector<PatchEdge>& patch_edge_array = p_chart_creator->GetPatchEdgeArray();

		/// find the common edge of these two charts
		const std::vector<int>& patch_edge_array_1 = patch_1.m_edge_index_array;
		const std::vector<int>& patch_edge_array_2 = patch_2.m_edge_index_array;

		std::vector<int> common_edges;
		for(size_t k=0; k<patch_edge_array_1.size(); ++k){
			if(find(patch_edge_array_2.begin(), patch_edge_array_2.end(), patch_edge_array_1[k]) 
				!= patch_edge_array_2.end()) common_edges.push_back(patch_edge_array_1[k]);
		}

		if(common_edges.size() >=2 ) return true;
		return false;
	}
}
