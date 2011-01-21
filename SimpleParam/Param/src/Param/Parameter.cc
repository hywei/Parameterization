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
#include <hj_3rd/hjlib/math/blas_lapack.h>
#include <hj_3rd/zjucad/matrix/lapack.h>

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

	void Parameter::OptimizeAmbiguityPatch()
	{
		if(!p_chart_creator) return;
		//p_chart_creator->OptimizeAmbiguityPatchShape();
		p_chart_creator->OptimizePatchShape();
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


		CMeshSparseMatrix lap_mat;
		SetLapMatrixCoef(p_mesh, lap_mat);		

		CMeshSparseMatrix lap_mat_with_mean_value;
		SetLapMatrixCoefWithMeanValueCoord(p_mesh, lap_mat_with_mean_value);

		size_t face_num = (size_t)p_mesh->m_Kernel.GetModelInfo().GetFaceNum();
		m_stiffen_weight.clear();
		m_stiffen_weight.resize(face_num, 1.0);

		m_flippd_face.clear();
		m_flippd_face.resize(face_num, false);

		int loop_num = 6;
		for(int k=0; k<loop_num; ++k)
		{						
			CMeshSparseMatrix lap_mat_with_stiffen;
//			SetLapMatrixWithStiffeningWeight(lap_mat_with_stiffen);

			CMeshSparseMatrix meanvalue_lap_mat_with_stiffen;
//			SetMeanValueLapMatrixWithStiffeningWeight(meanvalue_lap_mat_with_stiffen);
//			SolveParameter(meanvalue_lap_mat_with_stiffen);
//			SolveParameter(lap_mat_with_stiffen);
//			SolveParameter(lap_mat_with_mean_value);
			SolveParameter(lap_mat_with_mean_value);
			if(k < loop_num)
			{
                GetOutRangeVertices(m_out_range_vert_array);
				AdjustPatchBoundary();
				ConnerRelocating();
// 				LocalStiffening();		
// 				CheckFlipedTriangle();
// 				if(m_fliped_face_array.size() == 0) break;
			}
			
		}	   		

		GetOutRangeVertices(m_out_range_vert_array);
		AdjustPatchBoundary();
//		VertexRelalaxation();
        GetOutRangeVertices(m_out_range_vert_array);
        
		ResetFaceChartLayout();
//		SetMeshFaceTextureCoord();
		
		SetChartVerticesArray();

 		ComputeDistortion();

// 
 		CheckFlipedTriangle();

		return true;
	}

	void Parameter::FixAdjustedVertex(bool with_conner /* = false */)
	{
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		int fix_num=0;
		int total_fix_num = 0;
		do{
			fix_num = 0;
			for(int vid=0; vid < vert_num; ++vid)
			{
				if(!with_conner && IsConnerVertex(vid)) continue;
				if(FixAdjustedVertex(vid)) fix_num++; 
			}
			total_fix_num += fix_num;
		}while(fix_num !=0);

		std::cout << "Fix " << total_fix_num << " adjusted vertices" <<std::endl;
	}

	bool Parameter::FixAdjustedVertex(int vid)
	{
		PolyIndexArray& vtxAdjVtxArray = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
		IntArray& vAdjVtx = vtxAdjVtxArray[vid];
		int chart_id = m_vert_chart_array[vid];

		// check this vertex's neighbor vertex.
		std::set<int> vtxNeighChartSet;
		std::vector<int> vtxNeighChartVector;

		for (size_t j = 0 ; j < vAdjVtx.size(); j++)
		{
			int adjVtxGroup = m_vert_chart_array[vAdjVtx[j]];
			vtxNeighChartSet.insert(adjVtxGroup);
			vtxNeighChartVector.push_back(adjVtxGroup);
		}

		// if only 1 neighbor chart, this vertex is island vertex, or insider vertex.
		if (vtxNeighChartSet.size() == 1)
		{
			int adj_chart_id = vtxNeighChartVector[0];

			// if this vertex is island vertex.
			if (adj_chart_id != chart_id)
			{
				ParamCoord old_param_coord = m_vert_param_coord_array[vid];
				TransParamCoordBetweenCharts(chart_id, adj_chart_id, vid, old_param_coord, m_vert_param_coord_array[vid]);
				m_vert_chart_array[vid] = adj_chart_id;
				return true;
			}
		}else{
			// if all neighbor chart are others, this vertex island vertex
			// make it to one of them.
			if (vtxNeighChartSet.find(chart_id) == vtxNeighChartSet.end())
			{
				int max_appears_chart_id(-1), max_num(-1);
				for(size_t i=0; i<vtxNeighChartVector.size(); ++i)
				{
					int cid = vtxNeighChartVector[i];
					int appear_num = count(vtxNeighChartVector.begin(), vtxNeighChartVector.end(), cid);
					if(cid > max_num){ max_num = appear_num; max_appears_chart_id = cid;}
				}

				ParamCoord old_param_coord = m_vert_param_coord_array[vid];
				TransParamCoordBetweenCharts(chart_id, max_appears_chart_id, vid,
					old_param_coord, m_vert_param_coord_array[vid]);
				m_vert_chart_array[vid] = max_appears_chart_id;
				return true;
			}
		}

		return false;
	}

	void Parameter::ConnerRelocating()
	{
		std::vector<PatchConner>& conner_vec = p_chart_creator->GetPatchConnerArray();

		for(size_t i=0; i<conner_vec.size(); ++i)
		{
			PatchConner& conner = conner_vec[i];
			
			int conner_vid = conner.m_mesh_index;
 
			ParamCoord new_conner_coord; 
			ComputeConnerVertexNewParamCoord(conner_vid, new_conner_coord);

			std::vector<int> nb_vec;
			p_mesh->m_BasicOp.GetNeighborhoodVertex(conner_vid, 40, false, nb_vec);

			std::vector< std::pair<int, double> > node_candidate_vec;
			for (size_t j = 0; j < nb_vec.size(); j++)
			{
				int vid = nb_vec[j];
				Coord2D vert_st(m_vert_param_coord_array[vid].s_coord, m_vert_param_coord_array[vid].t_coord);
				if(vid == conner_vid){
					vert_st = Coord2D(new_conner_coord.s_coord, new_conner_coord.t_coord);
				}
				ParamCoord conner_pc;
				if(GetConnerParamCoord(m_vert_chart_array[vid], i, conner_pc))
				{
					Coord2D conner_st(conner_pc.s_coord, conner_pc.t_coord);
					double dis = (conner_st - vert_st).abs();
					node_candidate_vec.push_back(std::make_pair(vid, dis));
				}
			}

			int best_candidate_node(-1);
			double min_dist = std::numeric_limits<double>::max();
			for(size_t i=0; i<node_candidate_vec.size(); ++i)
			{
				if(min_dist > node_candidate_vec[i].second)
				{
					min_dist = node_candidate_vec[i].second;
					best_candidate_node = node_candidate_vec[i].first;
				}
			}
			std::cout <<conner.m_mesh_index << " " << best_candidate_node << std::endl;
			conner.m_mesh_index = best_candidate_node;

		}

	}

	void Parameter::ComputeConnerVertexNewParamCoord(int conner_vid, ParamCoord& new_pc) const
	{
		const PolyIndexArray& vAdjVertices = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
		const IndexArray& adj_vertics = vAdjVertices[conner_vid];

		std::vector<ParamCoord> adj_pc_vec(adj_vertics.size());
		int chart_id = m_vert_chart_array[conner_vid];

		new_pc.s_coord = 0.0; new_pc.t_coord = 0.0;
		for(size_t k=0; k<adj_vertics.size(); ++k)
		{
			int vid = adj_vertics[k];
			int cur_chart_id = m_vert_chart_array[vid];
			adj_pc_vec[k] = m_vert_param_coord_array[vid];
			if(cur_chart_id != chart_id)
			{
				TransParamCoordBetweenCharts(cur_chart_id, chart_id, vid, 
					m_vert_param_coord_array[vid], adj_pc_vec[k]);
			}
			new_pc.s_coord += adj_pc_vec[k].s_coord;
			new_pc.t_coord += adj_pc_vec[k].t_coord;
		}
		new_pc.s_coord /= adj_vertics.size();
		new_pc.t_coord /= adj_vertics.size();

	}

	bool Parameter::GetConnerParamCoord(int chart_id, int conner_idx, ParamCoord& conner_pc) const
	{
		const std::vector<ParamPatch>& patch_vec = p_chart_creator->GetPatchArray();
		const std::vector<ParamChart>& chart_vec = p_chart_creator->GetChartArray();			
		const ParamChart& chart = chart_vec[chart_id];
		const std::vector<int>& patch_conners = patch_vec[chart_id].m_conner_index_array;
		for(size_t i=0; i<patch_conners.size(); ++i)
		{
			if(conner_idx == patch_conners[i])
			{
				conner_pc = chart.m_conner_param_coord_array[i];
				return true;
			}
		}

		return false;
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

        m_face_patch_array = m_face_chart_array;
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

        m_vert_patch_array = m_vert_chart_array;
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
		int vari_num = SetVariIndexMapping(vari_index_mapping);
		vari_num *=2;
		
				
		std::cout << "Begin solve parameterization: variable num "<< vari_num << std::endl;

		LinearSolver linear_solver(vari_num);

		SetBoundaryVertexParamValue();
	 
		
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

		linear_solver.set_equation_div_flag();




		linear_solver.solve();



		for(int vid=0; vid < vert_num; ++vid)
		{
			int vari_index = vari_index_mapping[vid];
			if(vari_index != -1)
			{
				m_vert_param_coord_array[vid].s_coord = linear_solver.variable(vari_index*2).value();
				m_vert_param_coord_array[vid].t_coord = linear_solver.variable(vari_index*2+1).value();

				if(fabs(m_vert_param_coord_array[vid].s_coord) < LARGE_ZERO_EPSILON) m_vert_param_coord_array[vid].s_coord = 0.0;
				if(fabs(m_vert_param_coord_array[vid].s_coord - 1) < LARGE_ZERO_EPSILON) m_vert_param_coord_array[vid].s_coord = 1.0;
				if(fabs(m_vert_param_coord_array[vid].t_coord) < LARGE_ZERO_EPSILON) m_vert_param_coord_array[vid].t_coord = 0.0;
				if(fabs(m_vert_param_coord_array[vid].t_coord - 1) < LARGE_ZERO_EPSILON) m_vert_param_coord_array[vid].t_coord = 1.0;
			}
		}

        ofstream fout ("parame.txt");
        for(size_t k=0; k<m_vert_param_coord_array.size(); ++k){
            fout <<  m_vert_param_coord_array[k].s_coord << ' ' <<
                m_vert_param_coord_array[k].t_coord << std::endl;
        }
        fout.close();

	}

	int Parameter::SetVariIndexMapping(std::vector<int>& vari_index_mapping)
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
			if(p_mesh->m_BasicOp.IsBoundaryVertex(vid)
				|| find(conner_vertices.begin(), conner_vertices.end(), vid) != conner_vertices.end()){
				continue;
			}
			vari_index_mapping[vid] = vari_num++;
		}

		return vari_num;
	}

	void Parameter::AdjustPatchBoundary()
	{
		const std::vector<ParamChart>& param_chart_array = p_chart_creator->GetChartArray();
		const PolyIndexArray& vert_adjvertices_array =p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();

		bool swap_able = false, tag=false;
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
				double out_range_error ;
				if(param_chart.m_conner_param_coord_array.size() == 3){				
					out_range_error = ComputeOutRangeError4TriangleChart(param_coord, chart_id);
				}else if(param_chart.m_conner_param_coord_array.size() == 4)
					out_range_error = ComputeOutRangeError4Square(param_coord);

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
						double cur_out_range_error ;
						if(adj_chart.m_conner_param_coord_array.size() == 3)
							cur_out_range_error = ComputeOutRangeError4TriangleChart(param_coord, adj_chart_id);						  
						else if(adj_chart.m_conner_param_coord_array.size() == 4)
							cur_out_range_error = ComputeOutRangeError4Square(param_coord);
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
					
// 					m_vert_chart_array[vid] = min_error_chart_id;
// 					m_vert_param_coord_array[vid] = min_error_param_coord;
				}else
				{
					adjust_num ++;
				}
			}
			std::cout << "Adjust " << adjust_num << "vertices." << std::endl;
			
			
		}while(swap_able);

		std::vector<int> out_range_vertices;
		GetOutRangeVertices(out_range_vertices);

		FixAdjustedVertex(true);
		
	}

	void Parameter::VertexRelalaxation()
	{
		vector<int> out_range_vert_array;
		GetOutRangeVertices(out_range_vert_array);		

		int adjust_vert_num = 0;
		for(size_t k=0; k<out_range_vert_array.size(); ++k)
		{
			int out_range_vert = out_range_vert_array[k];
			if(!FindValidChartForOutRangeVertex(out_range_vert, 1))
			{
				std::cout<<"Can't find valid chart for this our range vertex " << out_range_vert << std::endl;
			}else
			{
				adjust_vert_num ++;
			}
		}

		std::cout<<"Adjust " << adjust_vert_num <<" vertices.\n";
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
				double cur_out_range_error;
				if(cur_param_chart.m_conner_param_coord_array.size() == 4)
					cur_out_range_error = ComputeOutRangeError4Square(cur_param_coord);
				else if(cur_param_chart.m_conner_param_coord_array.size() == 3)
					cur_out_range_error = ComputeOutRangeError4TriangleChart(cur_param_coord, cur_chart_id);
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
			//m_vert_chart_array[out_range_vert] = min_error_valid_chart_id;
			//m_vert_param_coord_array[out_range_vert] = min_error_valid_param_coord;
			return false;
		}
		return false;
	}

// 	void Parameter::TransParamCoordBetweenCharts(int from_chart_id, int to_chart_id, 
// 		const ParamCoord& from_param_coord, ParamCoord& to_param_coord) const
// 	{	
// 		TransFunctor tran_functor(p_chart_creator);
// 		zjucad::matrix::matrix<double> tran_mat = 
// 			tran_functor.GetTransMatrix(from_chart_id, to_chart_id);
// 
// 		to_param_coord.s_coord = tran_mat(0, 0)*from_param_coord.s_coord + 
// 			tran_mat(0, 1)*from_param_coord.t_coord + tran_mat(0, 2);
// 		to_param_coord.t_coord = tran_mat(1, 0)*from_param_coord.s_coord +
// 			tran_mat(1, 1)*from_param_coord.t_coord + tran_mat(1, 2);
// 	}

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
				tran_functor.GetTransMatrix(from_chart_id, to_chart_id, vid);

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
			return tran_functor.GetTransMatrix(from_vid, from_chart_id, to_vid, to_chart_id);
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
 						if(param_chart.m_conner_param_coord_array.size() == 3)
 							out_range_error += ComputeOutRangeError4TriangleChart(cur_param_coord, cur_chart_id);
 						else if(param_chart.m_conner_param_coord_array.size() == 4)
 							out_range_error += ComputeOutRangeError4Square(cur_param_coord);

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

		//std::cout << "There are " << m_unset_layout_face_array.size() << "unset faces." << std::endl;

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
           // 		face_chart_id = FindBestChartIDForTriShape(fid);

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
		int colors[56][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

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

			vtx_color_array[vid] = Color(colors[chart_id%56][0], colors[chart_id%56][1], colors[chart_id%56][2]);
		}
		
	}

    double Parameter::ComputeOutRangeError4Square(ParamCoord param_coord) const
    {
        double error=0;
        if(GreaterEqual(param_coord.s_coord, 1) || LessEqual(param_coord.s_coord, 0))
		{
            if(LessEqual(param_coord.s_coord, 0)) error += param_coord.s_coord*param_coord.s_coord;
            else if(GreaterEqual(param_coord.s_coord, 1)) error += (param_coord.s_coord-1)*(param_coord.s_coord-1);
        }
		if(GreaterEqual(param_coord.t_coord, 1) || LessEqual(param_coord.t_coord, 0))	
		{
            if(LessEqual(param_coord.t_coord, 0)) error += param_coord.t_coord*param_coord.t_coord;
            else if(GreaterEqual(param_coord.t_coord, 1)) error += (param_coord.t_coord-1)*(param_coord.t_coord-1);
        }
        return sqrt(error);
    }

	double Parameter::ComputeOutRangeError4TriangleChart(ParamCoord param_coord, int chart_id)
	{
		const ParamChart& chart = (p_chart_creator->GetChartArray())[chart_id];
		Coord2D p(param_coord.s_coord, param_coord.t_coord);
		Coord2D a(chart.m_conner_param_coord_array[0].s_coord, chart.m_conner_param_coord_array[0].t_coord);
		Coord2D b(chart.m_conner_param_coord_array[1].s_coord, chart.m_conner_param_coord_array[1].t_coord);
		Coord2D c(chart.m_conner_param_coord_array[2].s_coord, chart.m_conner_param_coord_array[2].t_coord);
		return DistanceToTriangle(p, a, b, c);
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
			const std::vector<ParamPatch>& patch_vec = p_chart_creator->GetPatchArray();
			const ParamPatch& param_patch = patch_vec[chart_id];
			const std::vector<int>& nb_patchs = param_patch.m_nb_patch_index_array;
			for(size_t i=0; i<nb_patchs.size(); ++i)
			{
				int adj_chart_id = nb_patchs[i];
				if(FindCorrespondingInChart(chart_param_coord, adj_chart_id, surface_coord)) break;
			}
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
							TransParamCoordBetweenCharts(m_vert_chart_array[faces[j]], chart_id, faces[j], 
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

		ofstream fout("distortion.txt");
		for(size_t k=0; k<face_isometric_distortion.size(); ++k){
			fout << face_isometric_distortion[k] << std::endl;
		}
		fout.close();

		//FaceValue2VtxColor(p_mesh, face_harmonic_distortion);
		std::vector<double> face_value = face_isometric_distortion;
// 		double sum_value = 0.0;
// 		for(size_t k=0; k<face_value.size(); ++k)
// 		{
// 			sum_value += face_value[k]*face_value[k];
// 		}
// 		sum_value = sqrt(sum_value);
// 		for(size_t k=0; k<face_value.size(); ++k)
// 		{
// 			face_value[k] /= sum_value;
// 		}
		FaceValue2VtxColor(p_mesh, face_value);
		
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
			if(area < 0 ){
				m_fliped_face_array.push_back(fid);
//				std::cout << fid << " : " << m_stiffen_weight[fid] << std::endl;
			}
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
	
	bool Parameter::IsConnerVertex(int vert_vid) const
	{
		if(p_chart_creator == NULL) return true;
		const std::vector<PatchConner>& conners = p_chart_creator->GetPatchConnerArray();

		for(size_t k=0; k<conners.size(); ++k)
		{
			if(conners[k].m_mesh_index == vert_vid) return true;
		}
		return false;
	}

	void Parameter::ComputeFaceSignFuncValue(std::vector<int>& face_sign_func_value)// const
	{
		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();
		face_sign_func_value.clear(); 
		face_sign_func_value.resize(face_num, 0);

		const PolyIndexArray& face_list = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		int flip_tri_num = 0;
		for(int k=0; k<face_num; ++k)
		{
			const IndexArray& face = face_list[k];
			int face_chart_id = m_face_chart_array[k];
			std::vector<double> u(3), v(3);
			for(int i=0; i<3; ++i)
			{
				int vid = face[i];
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
			
			double f_value = (u[1]-u[0]) * (v[2]-v[0]) - (u[2]-u[0]) * (v[1]-v[0]);
			face_sign_func_value[k] = SignFunc(f_value);
			if(face_sign_func_value[k] == -1){
				flip_tri_num ++;
				
			}
		}
		std::cout << "Fliped triangles number: " << flip_tri_num << std::endl;
	}

	void Parameter::ComputeFaceLocalDistortion(const std::vector<int>& face_sign_func_value, 
		std::vector<double>& face_distortion) const
	{
		TriDistortion tri_distortion(*this);
		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();
		face_distortion.clear(); face_distortion.resize(face_num);
		for(int k=0; k<face_num; ++k)
		{
			zjucad::matrix::matrix<double> jaobi_mat = tri_distortion.ComputeParamJacobiMatrix(k);			
			zjucad::matrix::matrix<double> U, S, VT;
			svd(jaobi_mat, U, S, VT);
			double s1 = S(0, 0), s2 = S(1, 1);
			face_distortion[k] = fabs( face_sign_func_value[k]*s1 - 1) 
				+ fabs( face_sign_func_value[k]*s2 - 1);
		}
	}


	void Parameter::UpdateStiffeningWeight(const std::vector<double>& face_distortion)
	{
		const PolyIndexArray& face_list = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		size_t face_num = face_list.size();
		
		for(size_t fid=0; fid<face_num; ++fid)
		{
			const IndexArray& face = face_list[fid];
			double delta_d = face_distortion[fid];
			delta_d = 0.0;
			for(int i=0; i<3; ++i)
			{
				int vid1 = face[i], vid2 = face[(i+1)%3];
				std::vector<int> adj_faces = GetMeshEdgeAdjFaces(p_mesh, vid1, vid2);
				assert(adj_faces.size() == 2);
				int _fid = (adj_faces[0] == fid) ? adj_faces[1] : adj_faces[0];
				delta_d += face_distortion[_fid];
			}
			delta_d /= 3.0;

			//delta_d = face_distortion[fid];

			m_stiffen_weight[fid] += std::min(delta_d, 15.0);
		}

		ofstream fout("Stiffen-Weight.txt");
		for(size_t fid=0; fid<face_num; ++fid)
		{
			fout << m_stiffen_weight[fid] << std::endl;
		}
		fout.close();

	}

	void Parameter::LocalStiffening()
	{						
		std::vector<int> face_sign_func_value;
		ComputeFaceSignFuncValue(face_sign_func_value);
		
 		std::vector<double> face_distortion;
 		ComputeFaceLocalDistortion(face_sign_func_value, face_distortion);

		const PolyIndexArray& face_list = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const PolyIndexArray& adjface_list = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();

		UpdateStiffeningWeight(face_distortion);

		size_t face_num = face_list.size();

		for(size_t k=0; k<face_num; ++k)
		{
			const IndexArray& face = face_list[k];
			
		}

		for(size_t k=0; k < face_num; ++k)
		{
			const IndexArray& face = face_list[k];
			if(face_sign_func_value[k] < 0 || m_flippd_face[k] == true)
			{
				m_flippd_face[k] = true;
				m_stiffen_weight[k] = (m_stiffen_weight[k] + 20000)/2;

				std::vector<int> nb_1_ring;
				p_mesh->m_BasicOp.GetFaceNeighborhood(k, nb_1_ring);

				for(size_t i=0; i<nb_1_ring.size(); ++i)
				{
					int nb1_fid = nb_1_ring[i];
					if(face_sign_func_value[nb1_fid] >= 0 && m_flippd_face[nb1_fid] == false){
						m_stiffen_weight[nb1_fid] = (m_stiffen_weight[nb1_fid] + 2000)/2;
					}

					std::vector<int> nb_2_ring;
					p_mesh->m_BasicOp.GetFaceNeighborhood(nb1_fid, nb_2_ring);
					for(size_t j=0; j<nb_2_ring.size(); ++j)
					{
						int nb2_fid = nb_2_ring[j];
						if(find(nb_1_ring.begin(), nb_1_ring.end(), nb2_fid) != nb_1_ring.end())
						{
							if(face_sign_func_value[nb2_fid] >=0 && m_flippd_face[nb2_fid] == false){
								m_stiffen_weight[nb2_fid] = (m_stiffen_weight[nb2_fid] + 200) /2;
							}
						}

						std::vector<int> nb_3_ring;
						p_mesh->m_BasicOp.GetFaceNeighborhood(nb2_fid, nb_3_ring);
						for(size_t k=0; k<nb_3_ring.size(); ++k)
						{
							if(find(nb_1_ring.begin(), nb_1_ring.end(), nb1_fid) != nb_1_ring.end() &&
								find(nb_2_ring.begin(), nb_2_ring.end(), nb2_fid) != nb_2_ring.end())
							{
								int nb3_fid = nb_3_ring[k];
								if(face_sign_func_value[nb3_fid] >= 0 && m_flippd_face[nb3_fid] == false){
									m_stiffen_weight[nb3_fid] = (m_stiffen_weight[nb3_fid] + 20) / 2;
								}
							}
						}
					}
				}


// 				for(int i=0; i<3; ++i)
// 				{
// 					int vid = face[i];
// 					const IndexArray& adjfaces = adjface_list[vid];
// 
// 					for(size_t j=0; j<adjfaces.size(); ++j)
// 					{
// 						int adj_fid = adjfaces[j];
// 						if(face_sign_func_value[adj_fid] >=0)
// 							m_stiffen_weight[adjfaces[j]] = 500;
// 					}
// 				}
			}else{
				m_stiffen_weight[k] = 1.0;
			}
		} 		

	}

	void Parameter::SetLapMatrixWithStiffeningWeight(CMeshSparseMatrix& stiffen_lap_mat)
	{
		std::vector<Coord> cot_coef_vec;
		SetCotCoef(p_mesh, cot_coef_vec);

		const PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		int i, j, k, n;
		int fID, vID;
		int row, col;		

		size_t vert_num = p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();
		//
		int nb_variables_ = (int) vert_num;
		stiffen_lap_mat.SetRowCol(nb_variables_, nb_variables_);

		for(i = 0; i < nb_variables_; ++ i)
		{
			const IndexArray& adjFaces = vAdjFaces[i];
			row = i;
			n = (int) adjFaces.size();

			for (j = 0; j < n; j++)
			{
				fID = adjFaces[j];
				const IndexArray& f = fIndex[fID];

				// Find the position of vertex i in face fID
				for(k = 0; k < 3; ++ k)
				{
					if(f[k] == i) break;
				}
				assert(k!=3);

				vID = f[(k+1)%3];
				col = vID;

				double value;
				double stiffen_weight = m_stiffen_weight[fID];
				double coef1 = cot_coef_vec[fID][(k+2)%3] * stiffen_weight;
				stiffen_lap_mat.GetElement(row, col, value);
				value -= coef1;
				stiffen_lap_mat.SetElement(row, col, value);

				vID = f[(k+2)%3];
				col = vID;

				double coef2 = cot_coef_vec[fID][(k+1)%3] * stiffen_weight;
				stiffen_lap_mat.GetElement(row, col, value);
				value -= coef2;
				stiffen_lap_mat.SetElement(row, col, value);

				stiffen_lap_mat.GetElement(row, row, value);
				value += coef1 + coef2; 
				stiffen_lap_mat.SetElement(row, row, value);
			}
		}

	}

	void Parameter::SetMeanValueLapMatrixWithStiffeningWeight(CMeshSparseMatrix& lap_matrix)
	{
		if(p_mesh == NULL) return;
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		std::vector< std::vector<double> > vert_tan_coef(vert_num);
		const PolyIndexArray& vert_adj_faces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const PolyIndexArray& vert_adj_vertices = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_list = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		for(int vid=0; vid < vert_num; ++vid){
			if(p_mesh->m_BasicOp.IsBoundaryVertex(vid)) continue;
			const IndexArray& adj_vert_array = vert_adj_vertices[vid];
			for(size_t i=0; i<adj_vert_array.size(); ++i)
			{
				int cur_vtx = adj_vert_array[i];
				int nxt_vtx = adj_vert_array[(i+1)%adj_vert_array.size()];
				Coord e1 = vCoord[cur_vtx] - vCoord[vid];
				Coord e2 = vCoord[nxt_vtx] - vCoord[vid];

				double a = angle(e1, e2)/2;
				vert_tan_coef[vid].push_back(tan(a));
			}
		}	   

		//
		int nb_variables_ = vert_num;
		lap_matrix.SetRowCol(nb_variables_, nb_variables_);
		

		for(int i = 0; i < nb_variables_; ++ i)
		{
			if(p_mesh->m_BasicOp.IsBoundaryVertex(i)) continue;
			const IndexArray& adjVertices = vert_adj_vertices[i];
			const IndexArray& adjFaces = vert_adj_faces[i];
			int row = i;

			int adj_num = adjVertices.size();
			for(int j=0; j<adj_num; ++j)
			{
				int fid=-1;
				{
					int vid1 = adjVertices[j];
					int vid2 = adjVertices[(j+adj_num-1)%adj_num];
					int vid3 = i;
					if(vid1 > vid2) swap(vid1, vid2);
					if(vid2 > vid3) swap(vid2, vid3);
					if(vid1 > vid2) swap(vid1, vid2);

					for(size_t k=0; k<adjFaces.size(); ++k)
					{
						const IndexArray& face = face_list[adjFaces[k]];
						int _vid1 = face[0], _vid2 = face[1], _vid3 = face[2];
						if(_vid1 > _vid2) swap(_vid1, _vid2);
						if(_vid2 > _vid3) swap(_vid2, _vid3);
						if(_vid1 > _vid2) swap(_vid1, _vid2);

						if(vid1 == _vid1 && vid2 == _vid2 && vid3 == _vid3)
						{
							fid = adjFaces[k]; break;
						}
					}
				}
				if(fid == -1){
					std::cerr <<"error, fid == -1" << std::endl;
				}
				
				double stiffen_weight = m_stiffen_weight[fid];			   
				int col = adjVertices[j];
				double edge_len = (vCoord[row]-vCoord[col]).abs();
				double coef = (vert_tan_coef[i][j] + vert_tan_coef[i][(j+adj_num-1)%adj_num])/edge_len;
				coef *= stiffen_weight;
				double cur_value;
				lap_matrix.GetElement(row, col, cur_value);
				cur_value -= coef;
				lap_matrix.SetElement(row, col, cur_value);

				lap_matrix.GetElement(row, row, cur_value);
				cur_value += coef;
				lap_matrix.SetElement(row, row, cur_value);

			}

		}
	}

	std::vector<ParamCoord> Parameter::GetFaceVertParamCoord(int fid) const
	{
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& face= fIndex[fid];

		std::vector<ParamCoord> pc_vec(3);
		int chart_id = m_face_chart_array[fid];
		for(int i=0; i<3; ++i)
		{
			int vid = face[i];
			int cur_chart_id = m_vert_chart_array[vid];
			pc_vec[i] = m_vert_param_coord_array[vid];
			if(cur_chart_id != chart_id)
			{
				TransParamCoordBetweenCharts(cur_chart_id, chart_id, vid, 
					m_vert_param_coord_array[vid], pc_vec[i]);
			}
		}
		return pc_vec;
	}

	void ChartCreator::AlginQuadPatchOrient()
	{
		int origin_chart_id = 0;
		std::set<int> aligned_chart;
		std::queue<int> q;
		q.push(origin_chart_id);
		aligned_chart.insert(origin_chart_id);

		map< int, int> orient_map;

		while(!q.empty())
		{
			int pid = q.front(); q.pop();
			const ParamPatch& patch = m_patch_array[pid];
			const std::vector<int>& nb_charts = patch.m_nb_patch_index_array;
			const std::vector<int>& edges = patch.m_edge_index_array;
			for(size_t i=0; i<4; ++i)
			{
				int e_id = edges[i];
				int nb_pid = nb_charts[i];
				if(aligned_chart.find(nb_pid) != aligned_chart.end()) continue;

				ParamPatch& nb_patch = m_patch_array[nb_pid];
				std::vector<int>& nb_edges = nb_patch.m_edge_index_array;
				int ps = find(nb_edges.begin(), nb_edges.end(), e_id) - nb_edges.begin();
			}
		}
	}

}
