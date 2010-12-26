#include "QuadChartCreator.h"
#include "QuadParameter.h"
#include "QuadTransFunctor.h"
#include "QuadDistortion.h"
#include "Barycentric.h"

#include "../ModelMesh/MeshModel.h"
#include "../Numerical/linear_solver.h"
#include "../Numerical/MeshSparseMatrix.h"
#include "../Common/HSVColor.h"

#include <hj_3rd/zjucad/matrix/matrix.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <queue>
#include <limits>


namespace PARAM
{
	QuadParameter::QuadParameter(boost::shared_ptr<MeshModel> _p_mesh) : p_mesh(_p_mesh){}
	QuadParameter::~QuadParameter(){}

	bool QuadParameter::LoadQuadFile(const std::string& file_name)
	{
		p_quad_chart_creator = boost::shared_ptr<QuadChartCreator> ( new QuadChartCreator(p_mesh));
		p_quad_chart_creator->LoadQuadFile(file_name);
		if(!p_quad_chart_creator->FormParamQuadCharts())
		{
			std::cout<<"Error: Cannot compute parameteriztion!\n";
			return false;
		}
		return true;
	}

	bool QuadParameter::ComputeParamCoord()
	{
		if(p_mesh == NULL) return false;
		if(p_quad_chart_creator == NULL)
		{
			std::cout<<"Error : Please load quad file first!\n";
			return false;
		}

		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		m_face_chart_array.clear(); m_face_chart_array.resize(face_num);
		m_vert_chart_array.clear(); m_vert_chart_array.resize(vert_num);

		SetInitFaceChartLayout();
		SetInitVertChartLayout();
	

		CMeshSparseMatrix lap_mat;
		SetLapMatrixCoef(p_mesh, lap_mat);

		m_vert_param_coord_array.clear();
		m_vert_param_coord_array.resize(vert_num);

		int loop_num = 5;
		for(int k=0; k<loop_num; ++k)
		{
			SolveQuadParameter(lap_mat);
			if(k < loop_num-1)
			{
				AdjustPatchBoundary();
				//VertexRelatextion();
			}
		}	   		

		ResetFaceChartLayout();
		SetMeshFaceTextureCoord();

		GetOutRangeVertices(m_out_range_vert_array);

		ComputeDistortion();

		SetChartVerticesArray();

		return true;
	}

	void QuadParameter::SetInitFaceChartLayout()
	{
		const std::vector<QuadPatch>& patch_array = p_quad_chart_creator->GetQuadPatchArray();
		for(size_t k=0; k<patch_array.size(); ++k)
		{
			const QuadPatch& quad_patch = patch_array[k];
			const std::vector<int>& faces_in_patch = quad_patch.m_face_index_array;
			for(size_t i=0; i<faces_in_patch.size(); ++i)
			{
				int fid = faces_in_patch[i];
				m_face_chart_array[fid] = k;
			}
		}

		if(p_mesh == NULL)
			return ;

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

		size_t face_num = p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();
		ColorArray& colorArray = p_mesh->m_Kernel.GetFaceInfo().GetColor();
		colorArray.clear();
		colorArray.resize(face_num);

		for(size_t i = 0; i < m_face_chart_array.size(); ++ i)
		{
			int index = m_face_chart_array[i] % 48;
			colorArray[i] = Color(colors[index][0], colors[index][1], colors[index][2]);
		}
	}

	void QuadParameter::SetInitVertChartLayout()
	{
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
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
	}
   
	void QuadParameter::SetBoundaryVertexParamValue()
	{
		const std::vector<QuadPatch>& patch_array = p_quad_chart_creator->GetQuadPatchArray();
		const std::vector<QuadChart>& chart_array = p_quad_chart_creator->GetQuadChartArray();
		const std::vector<PatchEdge>& patch_edge_array = p_quad_chart_creator->GetPatchEdgeArray();		
		const std::vector<PatchConner>& patch_conner_array = p_quad_chart_creator->GetPatchConnerArray();



		for(size_t k=0; k<patch_array.size(); ++k)
		{
			int chart_id = k;
			const QuadPatch& quad_patch = patch_array[k];
			const QuadChart& quad_chart = chart_array[k];

			const std::vector<int>& quad_conners = quad_patch.m_conner_index_array;
			const std::vector<int>& quad_edges = quad_patch.m_edge_index_array;

			/// fix conner vertex
			for(size_t i=0; i<quad_conners.size(); ++i)
			{
				int conner_idx = quad_conners[i];
				int conner_vid = patch_conner_array[conner_idx];
				if(m_vert_chart_array[conner_vid] == chart_id) 
				{
					m_vert_param_coord_array[conner_vid].s_coord =
						quad_chart.m_conner_param_coord_array[i].s_coord;
					m_vert_param_coord_array[conner_vid].t_coord =
						quad_chart.m_conner_param_coord_array[i].t_coord;				
				}
			}

			/// fix boundary vertex
			for(size_t i=0; i<quad_edges.size(); ++i)
			{
				int edge_idx = quad_edges[i];
				const PatchEdge& patch_edge = patch_edge_array[edge_idx];

				if(patch_edge.m_neighbor_patch_array.size() == 1) /// this is a boundary patch edge
				{
					std::pair<int, int> conner_pair = patch_edge.GetConnerPair();

					int conner_idx_1 = GetConnerIndexInPatch(conner_pair.first, k);
					int conner_idx_2 = GetConnerIndexInPatch(conner_pair.second, k);


					double start_s_coord = quad_chart.m_conner_param_coord_array[conner_idx_1].s_coord;
					double start_t_coord = quad_chart.m_conner_param_coord_array[conner_idx_1].t_coord;
					double end_s_coord = quad_chart.m_conner_param_coord_array[conner_idx_2].s_coord;
					double end_t_coord = quad_chart.m_conner_param_coord_array[conner_idx_2].t_coord;

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
					}
				}
			}

		}
	}

	void QuadParameter::SolveQuadParameter(const CMeshSparseMatrix& lap_mat)
	{
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		
		vector<int> vari_index_mapping;
		SetVariIndexMapping(vari_index_mapping);
		SetBoundaryVertexParamValue();
		
		int vari_num = (int)vari_index_mapping.size()*2;

		LinearSolver linear_solver(vari_num);
		
		QuadTransFunctor quad_trans_functor(p_quad_chart_creator);

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

				double right_b = 0;
				for(size_t k=0; k<row_index.size(); ++k)
				{
					int col_vert = row_index[k];
					int from_chart_id = m_vert_chart_array[col_vert];								    

					int var_index = vari_index_mapping[col_vert];

					if(from_chart_id == to_chart_id)
					{
						if( var_index == -1)
						{
							double st_value = (st==0) ? m_vert_param_coord_array[col_vert].s_coord : 
								m_vert_param_coord_array[col_vert].t_coord;
							right_b -= row_data[k]*st_value;
						}else
						{
							linear_solver.add_coefficient(var_index*2 + st, row_data[k]);
						}
					}else
					{
						if(var_index == -1)
						{
							ParamCoord param_coord;
							TransParamCoordBetweenCharts(from_chart_id, to_chart_id, 
								m_vert_param_coord_array[col_vert], param_coord);
							double st_value = (st==0) ? param_coord.s_coord : param_coord.t_coord;
							right_b -= row_data[k]*st_value;
						}else
						{
							zjucad::matrix::matrix<double> trans_mat = quad_trans_functor.GetTransMatrix(from_chart_id, to_chart_id);
							double a = (st == 0) ? trans_mat(0, 0) : trans_mat(1, 0);
							double b = (st == 0) ? trans_mat(0, 1) : trans_mat(1, 1);
							double c = (st == 0) ? trans_mat(0, 2) : trans_mat(1, 2);
							assert( fabs(a) < LARGE_ZERO_EPSILON || fabs(b) < LARGE_ZERO_EPSILON);
							///! a*u + b*v + c
							if(!(fabs(a) < LARGE_ZERO_EPSILON))
							{
								linear_solver.add_coefficient(var_index*2, row_data[k]*a);
							}
							if(!(fabs(b) < LARGE_ZERO_EPSILON))
							{
								linear_solver.add_coefficient(var_index*2+1, row_data[k]*b);
							}
							right_b -= c*row_data[k];
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
		
	}

	void QuadParameter::SetVariIndexMapping(std::vector<int>& vari_index_mapping)
	{
		const std::vector<PatchConner>& patch_conner_array = p_quad_chart_creator->GetPatchConnerArray();

		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		vari_index_mapping.clear();
		vari_index_mapping.resize(vert_num, -1);

		int vari_num = 0;
		for(int vid=0; vid < vert_num; ++vid)
		{
			if(p_mesh->m_BasicOp.IsBoundaryVertex(vid) ||
				find(patch_conner_array.begin(), patch_conner_array.end(), vid) != patch_conner_array.end())
			{
				continue;
			}
			vari_index_mapping[vid] = vari_num++;
		}
	}

	void QuadParameter::AdjustPatchBoundary()
	{
		const std::vector<QuadChart>& quad_chart_array = p_quad_chart_creator->GetQuadChartArray();
		const PolyIndexArray& vert_adjvertices_array =p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();

		int adjust_num(0);
		for(int vid = 0; vid < (int)vert_adjvertices_array.size(); ++vid)
		{
			int chart_id = m_vert_chart_array[vid];
			ParamCoord param_coord = m_vert_param_coord_array[vid];
			const QuadChart& quad_chart = quad_chart_array[chart_id];
			if(quad_chart.InValidRangle(param_coord)) continue;

			bool flag = false;
			double out_range_error = ComputeOutRangeError(param_coord);
			double min_out_range_error = out_range_error;
			int min_error_chart_id = chart_id;
			ParamCoord min_error_param_coord;

			const IndexArray& adj_vertices = vert_adjvertices_array[vid];
			for(size_t k=0; k<adj_vertices.size(); ++k)
			{
				int adj_chart_id = m_vert_chart_array[adj_vertices[k]];
				if(chart_id == adj_chart_id) continue;

				ParamCoord adj_param_coord = m_vert_param_coord_array[adj_vertices[k]];
				const QuadChart& adj_quad_chart = quad_chart_array[adj_chart_id];
				if(!adj_quad_chart.InValidRangle(adj_param_coord)) continue;

				if(chart_id != adj_chart_id)
				{
					TransParamCoordBetweenCharts(chart_id, adj_chart_id, 
						m_vert_param_coord_array[vid], param_coord);
				}

				if(adj_quad_chart.InValidRangle(param_coord))
				{
					m_vert_chart_array[vid] = adj_chart_id;
					m_vert_param_coord_array[vid] = param_coord;
					flag = true; 
					break;
				}else
				{
					double cur_out_range_error = ComputeOutRangeError(param_coord);
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
// 				if(chart_id == min_error_chart_id) 
// 					std::cout<<"Can't adjust this vertex " << vid << std::endl;
				m_vert_chart_array[vid] = min_error_chart_id;
				m_vert_param_coord_array[vid] = min_error_param_coord;
			}else
			{
				adjust_num ++;
			}
		}
		std::cout << "Adjust " << adjust_num << "vertices." << std::endl;
	}

	void QuadParameter::VertexRelatextion()
	{
		vector<int> out_range_vert_array;
		GetOutRangeVertices(out_range_vert_array);		

		int adjust_vert_num = 0;
		for(size_t k=0; k<out_range_vert_array.size(); ++k)
		{
			int out_range_vert = out_range_vert_array[k];
			if(!FindValidChartForOutRangeVertex(out_range_vert))
			{
				std::cout<<"Can't find valid chart for this our range vertex " << out_range_vert << std::endl;
			}else
			{
				adjust_vert_num ++;
			}
		}

		std::cout<<"Adjust " << adjust_vert_num <<" vertices.\n";

		GetOutRangeVertices(out_range_vert_array);
		ResetFaceChartLayout();
	}

	void QuadParameter::GetOutRangeVertices(std::vector<int>& out_range_vert_array) const
	{
		const std::vector<QuadChart>& quad_chart_array = p_quad_chart_creator->GetQuadChartArray();

		out_range_vert_array.clear();
		for(size_t vid=0; vid<m_vert_param_coord_array.size(); ++vid)
		{
			int vert_chart_id = m_vert_chart_array[vid];
			const QuadChart& quad_chart = quad_chart_array[vert_chart_id];
			if(!quad_chart.InValidRangle(m_vert_param_coord_array[vid]))
			{
				out_range_vert_array.push_back(vid);
			}
		}

		std::cout<<"There are " << out_range_vert_array.size() << " out range vertices.\n";
	}

	bool QuadParameter::FindValidChartForOutRangeVertex(int out_range_vert, int max_ringe_num /* = 5 */)
	{
		int init_chart_id = m_vert_chart_array[out_range_vert];
		ParamCoord init_param_coord = m_vert_param_coord_array[out_range_vert];

		const std::vector<QuadChart>& quad_chart_array = p_quad_chart_creator->GetQuadChartArray();
		const std::vector<QuadPatch>& quad_patch_array = p_quad_chart_creator->GetQuadPatchArray();

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
			TransParamCoordBetweenCharts(init_chart_id, cur_chart_id, init_param_coord, cur_param_coord);
			const QuadChart& cur_quad_chart = quad_chart_array[cur_chart_id];
			if(cur_quad_chart.InValidRangle(cur_param_coord))
			{
				valid_chart_id = cur_chart_id; 
				valid_param_coord = cur_param_coord;
				break;
			}else
			{
				double cur_out_range_error = ComputeOutRangeError(cur_param_coord);
				if(cur_out_range_error < min_out_range_error)
				{
					min_out_range_error = cur_out_range_error;
					min_error_valid_param_coord = cur_param_coord;
					min_error_valid_chart_id = cur_chart_id;
				}
			}

			if(steps > max_ringe_num) break;

			const QuadPatch& cur_quad_patch = quad_patch_array[cur_chart_id];
			const std::vector<int>& chart_neighbor = cur_quad_patch.m_neighbor_patch_array;

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

	void QuadParameter::TransParamCoordBetweenCharts(int from_chart_id, int to_chart_id, 
		const ParamCoord& from_param_coord, ParamCoord& to_param_coord) const
	{
		QuadTransFunctor tran_functor(p_quad_chart_creator);
		
		zjucad::matrix::matrix<double> tran_mat = tran_functor.GetTransMatrix(from_chart_id, to_chart_id);

		to_param_coord.s_coord = tran_mat(0, 0)*from_param_coord.s_coord + 
			tran_mat(0, 1)*from_param_coord.t_coord + tran_mat(0, 2);
		to_param_coord.t_coord = tran_mat(1, 0)*from_param_coord.s_coord +
			tran_mat(1, 1)*from_param_coord.t_coord + tran_mat(1, 2);
	}



	double QuadParameter::ComputeMeshPathLength(const std::vector<int>& mesh_path, int start_idx, int end_idx) const
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

	int QuadParameter::GetConnerIndexInPatch(int conner_id, int patch_id)
	{
		const std::vector<QuadPatch>& patch_array = p_quad_chart_creator->GetQuadPatchArray();
		const std::vector<int>& conner_array = p_quad_chart_creator->GetPatchConnerArray();

		assert(patch_id < (int)patch_array.size());
		const QuadPatch& quad_patch = patch_array[patch_id];
		const std::vector<int>& patch_conner_array = quad_patch.m_conner_index_array;

		for(size_t k=0; k<patch_conner_array.size(); ++k)
		{
			int conner_idx = patch_conner_array[k];
			if(conner_array[conner_idx] == conner_id) return (int)k;
		}

		return -1;
	}

	void QuadParameter::ResetFaceChartLayout()
	{
		const std::vector<QuadChart>& quad_chart_array = p_quad_chart_creator->GetQuadChartArray();

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
			for(size_t j=0; j<3; ++j)
			{
				int chart_id = m_vert_chart_array[faces[j]];
				bool flag = true;
				for(size_t k=0; k<3; ++k)
				{
					int cur_chart_id = m_vert_chart_array[faces[k]];
					ParamCoord cur_param_coord = m_vert_param_coord_array[faces[k]];
					if(cur_chart_id != chart_id)
					{
						TransParamCoordBetweenCharts(cur_chart_id,chart_id, 
							m_vert_param_coord_array[faces[k]], cur_param_coord);
					}
					const QuadChart& quad_chart = quad_chart_array[cur_chart_id];
					if(!quad_chart.InValidRangle(cur_param_coord))
					{
						flag = false;
						break;
					}
				}
				if(flag == true)
				{
					std_chart_id = chart_id;
					break;
				}
			}
			if(std_chart_id == -1)
			{
				m_unset_layout_face_array.push_back(i);
				//std::cout<<"Can't find valid chart for this face's three vertices!\n";
				int max_times(-1);
				for(size_t k=0; k<chart_id_vec.size(); ++k)
				{
					int c = chart_id_vec[k];
					int cur_times = count(chart_id_vec.begin(), chart_id_vec.end(), c);
					if(cur_times > max_times)
					{
						max_times = cur_times;
						std_chart_id = c;
					}
				}
			}

			m_face_chart_array[i] = std_chart_id;
		}

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

	void QuadParameter::SetMeshFaceTextureCoord()
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
			face_chart_id = FindBestChartIDForTriShape(fid);

			for(int k=0; k<3; ++k)
			{
				int vid = faces[k];
				int vert_chart_id = m_vert_chart_array[vid];
				ParamCoord vert_param_coord = m_vert_param_coord_array[vid];
				if(vert_chart_id != face_chart_id)
				{
					TransParamCoordBetweenCharts(vert_chart_id, face_chart_id, 
						m_vert_param_coord_array[vid], vert_param_coord);
				}
				face_tex[k] = TexCoord(vert_param_coord.s_coord, vert_param_coord.t_coord);
			}
		}

	}

	void QuadParameter::SetChartVerticesArray()
	{
		int colors[8][3] = 
		{
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0}
		};

		int chart_num = p_quad_chart_creator->GetChartNumber();
		m_chart_vertices_array.clear();
		m_chart_vertices_array.resize(chart_num);

		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

// 		ColorArray& vtx_color_array = p_mesh->m_Kernel.GetVertexInfo().GetColor();
// 		vtx_color_array.clear(); vtx_color_array.resize(vert_num);

		for(int vid = 0; vid<vert_num; ++vid)
		{
			int chart_id = m_vert_chart_array[vid];
			m_chart_vertices_array[chart_id].push_back(vid);
//			vtx_color_array[vid] = Color(colors[chart_id%8][0], colors[chart_id%8][1], colors[chart_id%8][2]);
		}

		
	}

	int QuadParameter::FindBestChartIDForTriShape(int fid) const
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
					TransParamCoordBetweenCharts(cur_chart_id, chart_id, 
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

	void QuadParameter::ComputeDistortion() 
	{
		QuadDistortion quad_distortion(*this);
		quad_distortion.ComputeDistortion();

		const std::vector<double>& face_harmonic_distortion = quad_distortion.GetFaceHarmonicDistortion();
		const std::vector<double>& face_isometric_distortion = quad_distortion.GetFaceIsometricDistortion();

		FaceValue2VtxColor(face_harmonic_distortion);
	//	FaceValue2VtxColor(face_isometric_distortion);
	}

	bool QuadParameter::FindCorrespondingOnSurface(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const
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

	bool QuadParameter::FindCorrespondingInChart(const ChartParamCoord& chart_param_coord, 
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

	void QuadParameter::FaceValue2VtxColor(const std::vector<double>& face_value)
	{
		ColorArray& face_color_array = p_mesh->m_Kernel.GetFaceInfo().GetColor();
		face_color_array.clear();

		CHSVColor color;
		color.RGBtoHSV(1, 0, 0);
		color.m_S = 0.9f;
		color.m_V = 0.9f;

		size_t face_num = p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();

		double min = 1e20, max = -1e20;
		for(size_t i = 0; i < face_num; ++i)
		{
			if(face_value[i] < min) min = face_value[i];
			if(face_value[i] > max) max = face_value[i];
		}

		double range = (max - min) * 1.1;

		bool eq = ALMOST_EQUAL_LARGE(range, 0.0);

		for(size_t i = 0; i < face_value.size();++i)
		{
			float R, G, B;
			if (eq)
			{
				color.m_H = (float) 0.5 * 255;
			}
			else
			{
				double prop = (face_value[i] - min) / range;

				//color.m_S = prop;
				color.m_H = (float)(1 - prop) * 255;
			}
			color.HSVtoRGB(&R, &G, &B);
			Color c(R, G, B);
			face_color_array.push_back(c);
		}
	}
}
