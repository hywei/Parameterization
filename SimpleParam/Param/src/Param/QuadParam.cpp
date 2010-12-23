#include "QuadParam.h"
#include "../ModelMesh/MeshModel.h"
#include "../Numerical/linear_solver.h"
#include "../Common/HSVColor.h"
#include <fstream>
#include <queue>
#include <map>
#include <set>
#include <limits>
#include <gl/GLAux.h>
#include "TriangleTransFunctor.h"
#include "Barycentric.h"

#include "../hj_3rd/include/math/blas_lapack.h"
#include "../hj_3rd/include/zjucad/matrix/lapack.h"
using namespace std;

namespace PARAM
{
	QuadParam::QuadParam(){}

	QuadParam::~QuadParam(){}

	int QuadParam::AttchMesh(const boost::shared_ptr<MeshModel> mesh)
	{
		if( mesh == NULL)
			return -1;

		p_mesh = mesh;

		size_t vertex_num = p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();
		size_t face_num =  p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();

		m_vertex_group.clear();
		m_face_group.clear();

		m_vertex_group.resize(vertex_num, -1);
		m_face_group.resize(face_num, -1);

		return 0;
	}

	int QuadParam::LoadQuadFile(const std::string& quad_file_name)
	{

		ifstream fin(quad_file_name.c_str());
		if(fin.fail())
		{
			printf("Cannot load %s\n", quad_file_name.c_str());
			return -1;
		}

		m_chart_node_array.clear();
		m_chart_path_array.clear();

		int node_number, mesh_index;
		fin >> node_number;
		m_chart_node_array.resize(node_number);
		for(int k = 0; k < node_number; ++k)
		{
		   fin >> mesh_index;
			m_chart_node_array[k] = QuadNode(k, mesh_index);
		}

		int path_number, vertex_number;
		fin >> path_number;
		m_chart_path_array.resize(path_number);
		for(int k = 0; k < path_number; ++k)
		{
			m_chart_path_array[k] = QuadPath(k);
			std::vector<int>& mesh_path = m_chart_path_array[k].m_path;
			fin >> vertex_number;
			mesh_path.resize(vertex_number);
			for(int i = 0; i < vertex_number; ++i)
			{
				fin >> mesh_path[i];
			}
		}

		int chart_number;
		fin >> chart_number;
		m_chart_array.resize(chart_number);
		for(int k=0; k<chart_number; ++k)
		{
			m_chart_array[k].m_id = k;
			vector<int>& node_array = m_chart_array[k].m_node_index_array;
			vector<int>& path_array = m_chart_array[k].m_path_index_array;
			node_array.resize(4);
			path_array.resize(4);
			for(int i=0; i<4; ++i) 
			{
				fin >> node_array[i];
			}		
			for(int i=0; i<4; ++i)
			{
				fin >> path_array[i];
			}
		}


		fin.close();
		return 0;
	}
	
	bool QuadParam::IsConnerVertex(int vid) const
	{
		for(size_t k=0; k<m_chart_node_array.size(); ++k)
		{
			if(vid == m_chart_node_array[k].m_mesh_index) return true;
		}
		return false;
	}

	int QuadParam::FormQuadChart()
	{
		map< pair<int, int>, vector<int> > m_edge_pid_mapping;
		for(size_t k=0; k<m_chart_path_array.size(); ++k)
		{
			int path_id = m_chart_path_array[k].m_id;
			const vector<int>& q_path = m_chart_path_array[k].m_path;
			for(size_t i=1; i<q_path.size(); ++i)
			{
				int vtx1 = q_path[i-1];
				int vtx2 = q_path[i];

				vector<int>& pid_array1 = m_edge_pid_mapping[make_pair(vtx1, vtx2)];
				if(find(pid_array1.begin(), pid_array1.end(), path_id) == pid_array1.end())
					pid_array1.push_back(path_id);
				vector<int>& pid_array2 = m_edge_pid_mapping[make_pair(vtx2, vtx1)];
				if(find(pid_array2.begin(), pid_array2.end(), path_id) == pid_array2.end())
					pid_array2.push_back(path_id);

			}
		}
		
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		size_t face_num = face_index_array.size(); 

		vector<bool> face_visited_flag(face_num, false);
		int chart_id = 0;
		for(size_t k=0; k<face_num; ++k)
		{
			set<int> quad_path_set;
			if(face_visited_flag[k] == false)
			{
				std::set<int> face_set;
				if(FloodFillChart((int)k, face_set, quad_path_set, 
					face_visited_flag, m_edge_pid_mapping) == 0)
				{
					if(quad_path_set.size() != 4)
					{
						printf("quad path array size is %d, can't form chart!\n", quad_path_set.size());
						continue;
					}

					vector<int> quad_path_array;
					quad_path_array.assign(quad_path_set.begin(), quad_path_set.end());
					int chart_id = GetChartIdFromQuadPath(quad_path_array);

					for(set<int>::iterator is = face_set.begin(); is!=face_set.end(); ++is)
					{
						m_face_group[*is] = chart_id;
					}

				}
			}
		}

		printf("Form %d quad charts!\n", m_chart_array.size());

		SetVertexGroup();
		SetFaceGroupColor(m_face_group);

		/// set each chart's four conner's local coordinate
		m_chart_nodes_local_coord.clear(); 
		m_chart_nodes_local_coord.resize(m_chart_array.size());
		for(size_t cid = 0; cid < m_chart_array.size(); ++cid)
		{
			vector<Coord2D>& node_local_coord = m_chart_nodes_local_coord[cid];
		    node_local_coord.resize(4);
			node_local_coord[0] = Coord2D(0, 0); node_local_coord[1] = Coord2D(1, 0);
			node_local_coord[2] = Coord2D(1, 1); node_local_coord[3] = Coord2D(0, 1);
		}


		return 0;
	}

	void QuadParam::SetVertexGroup()
	{
		const PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();

		for(size_t k=0; k<m_vertex_group.size(); ++k)
		{
			const IndexArray& faces = vAdjFaces[k];
			map<int, int> chart_count;
			for(size_t i=0; i<faces.size(); ++i)
			{
				int fid = faces[i];
				int chart_id = m_face_group[fid];
				chart_count[chart_id] ++;
			}

			int chart_id(-1), num(-1);
			for(map<int, int>::iterator im = chart_count.begin(); im!=chart_count.end(); ++im)
			{
				if(num < im->second)
				{
					chart_id = im->first;
					num = im->second;
				}
			}

			m_vertex_group[k] = chart_id;
		}
	}

	int QuadParam::FormQuadHalfEdge()
	{
		m_quad_he_array.clear();
		/// form each chart's four half edges
		for(size_t k = 0; k<m_chart_array.size(); ++k)
		{
			FormChartHalfEdge(m_chart_array[k].m_id);
		}

		/// set pair for each half edge

		map<int, vector<int> > qpath_qhe_mapping; /// the mapping between a quad-edge and its half-edge
		for(size_t k=0; k<m_quad_he_array.size(); ++k)
		{
			QuadHalfEdge& q_he = m_quad_he_array[k];
			int q_path_idx = q_he.m_path_id;
			qpath_qhe_mapping[q_path_idx].push_back(k);
		}

		for(size_t k=0; k<m_quad_he_array.size(); ++k)
		{
			QuadHalfEdge& q_he = m_quad_he_array[k];
			if(q_he.m_pair_id == -1)
			{
				int q_path_idx = q_he.m_path_id;
				const vector<int>& he_index_array = qpath_qhe_mapping[q_path_idx];
				assert(he_index_array.size() == 1 || he_index_array.size() == 2);

				if(he_index_array.size() == 2)
				{
					int another_he_idx = (he_index_array[0] == k) ? he_index_array[1] :
						he_index_array[0];
					q_he.m_pair_id = another_he_idx;
					QuadHalfEdge& q_he_2 = m_quad_he_array[another_he_idx];
					q_he_2.m_pair_id = k;
				}
			}
		}

		return 0;
	}

	int QuadParam::FormChartNeighFunc()
	{
		m_chart_neigh_func.clear();
		m_chart_neigh_func.resize(m_chart_array.size()*4);

		for(size_t i=0; i<m_chart_array.size(); ++i)
		{
			const QuadChart& q_chart = m_chart_array[i];
			const vector<int>& c_he_array = q_chart.m_he_index_array;

			for(size_t j=0; j<c_he_array.size(); ++j)
			{
				ChartNeighFun fun;
				fun.m_from_chart_id = q_chart.m_id;

				int he_idx = c_he_array[j];
				const QuadHalfEdge& he = m_quad_he_array[he_idx];
				if(he.m_pair_id != -1)
				{
					const QuadHalfEdge& pair_he = m_quad_he_array[he.m_pair_id];
					fun.m_to_chart_id = pair_he.m_chart_id;
					fun.SetTransType(he.m_type, pair_he.m_type);
				}else
				{
					fun.m_to_chart_id = -1;
					fun.m_type = UN_DEFINEED;
				}
				m_chart_neigh_func[i*4+j] = fun;
			}
			
		}

		return 0;
	}


	int QuadParam::FormChartHalfEdge(int chart_id)
	{
		assert(chart_id < (int) m_chart_array.size());
		QuadChart& q_chart = m_chart_array[chart_id];
	
		assert(q_chart.m_node_index_array.size() == 4);		
		q_chart.m_he_index_array.clear();
		q_chart.m_he_index_array.resize(4);


		for(size_t k=0; k<4; ++k)
		{
			QuadHalfEdge new_he;
			new_he.m_id = m_quad_he_array.size();
			new_he.m_chart_id = chart_id;
			new_he.m_node_id = q_chart.m_node_index_array[k];
			new_he.m_path_id = q_chart.m_path_index_array[k];
			new_he.m_type = GetHalfEdgeType(k);

			m_quad_he_array.push_back(new_he);

			q_chart.m_he_index_array[k] = new_he.m_id;
		}

		for(size_t k=0; k<4; ++k)
		{
			int cur_he_idx = q_chart.m_he_index_array[k];
			int nxt_he_idx = q_chart.m_he_index_array[(k+1)%4];
			int prv_he_idx = q_chart.m_he_index_array[(k+3)%4];
			QuadHalfEdge& q_he = m_quad_he_array[cur_he_idx];
			q_he.m_next_id = nxt_he_idx;
			q_he.m_prev_id = prv_he_idx;
		}

		return 0;
	}

	int QuadParam::FloodFillChart(int face_id, std::set<int>& face_set, std::set<int>& quad_path_set,
		std::vector<bool>& face_visited_flag, const std::map< std::pair<int,int>, vector<int> >& edge_pid_mapping)
	{
		PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		queue<int> q;
		q.push(face_id);
		face_visited_flag[face_id] = true;

		while(!q.empty())
		{
			int fid = q.front(); q.pop();
			face_set.insert(fid);

			IndexArray& faces = face_index_array[fid];
			for(size_t k=0; k<faces.size(); ++k)
			{
				int vtx1 = faces[k];
				int vtx2 = faces[(k+1)%faces.size()];

				pair<int, int> edge = make_pair(vtx1, vtx2);
				map< pair<int, int>, vector<int> >::const_iterator im = edge_pid_mapping.find(edge);
				if(im != edge_pid_mapping.end())
				{
					if(im->second.size() == 1)
					{	
						int pid = (im->second)[0];
						quad_path_set.insert(pid);
					}
					continue;
				}

				const vector<int>& edge_adj_faces = GetEdgeAdjFaces(vtx1, vtx2);
				for(size_t i=0; i<edge_adj_faces.size(); ++i)
				{
					int adj_fid = edge_adj_faces[i];
					if(adj_fid != fid && face_visited_flag[adj_fid] == false)
					{
						q.push(adj_fid);
						face_visited_flag[adj_fid] = true;
					}
				}
			}
		}


		return 0;
	}

	void QuadParam::TransParamCoordBetweenTwoChart(int from_chart_id, int to_chart_id, 
		const ParamCoord& original_coord, ParamCoord& new_coord) const
	{
		double s_coord = original_coord.s_coord;
		double t_coord = original_coord.t_coord;
		TransSTBetweenTwoCharts(from_chart_id, to_chart_id, s_coord, t_coord);
		new_coord.s_coord = s_coord;
		new_coord.t_coord = t_coord;
	}

	zjucad::matrix::matrix<double> QuadParam::GetTransMatrixBetweenTwoCharts(int from_chart_id, int to_chart_id) const
	{
		zjucad::matrix::matrix<double> trans_mat(zjucad::matrix::eye<double>(3));
		if(from_chart_id == to_chart_id ) return trans_mat;
		vector<int> trans_list;
		if(!GetTranslistBetweenTwoCharts(from_chart_id, to_chart_id, trans_list))
		{
			printf("Cannot get the transition list!\n");
		}
		assert(trans_list.size() >=2);
		
		for(size_t k=1; k<trans_list.size(); ++k)
		{	
			zjucad::matrix::matrix<double> temp_trans_mat = trans_mat;
			trans_mat =  GetTransMatrixOfAdjCharts(trans_list[k-1], trans_list[k]) * temp_trans_mat;
		}

		return trans_mat;
	}

	zjucad::matrix::matrix<double> QuadParam::GetTransMatrixOfAdjCharts(int from_chart_id, int to_chart_id) const
	{
		const QuadChart& from_chart = m_chart_array[from_chart_id];
		const QuadChart& to_chart = m_chart_array[to_chart_id];

		/// find the common edge of these two charts
		ptrdiff_t com_edge_idx_1(-1), com_edge_idx_2(-1);
		for(size_t k=0; k<from_chart.m_path_index_array.size(); ++k)
		{
			int edge_idx_1 = from_chart.m_path_index_array[k];
			for(size_t i=0; i< to_chart.m_path_index_array.size(); ++i)
			{
				int edge_idx_2 = to_chart.m_path_index_array[i];
				if(edge_idx_1 == edge_idx_2)
				{
					com_edge_idx_1 = k;
					com_edge_idx_2 = i;
					break;
				}
			}
			if(com_edge_idx_1 != -1 && com_edge_idx_2 != -1) break;
		}

		pair<int, int> old_x_axis_1, new_x_axis_1;
		/// find the from chart and to chart's origin and x-axis
		old_x_axis_1 = make_pair(0, 1); 
		new_x_axis_1 = make_pair(com_edge_idx_1, (com_edge_idx_1+1)%4);

		pair<int, int> old_x_axis_2, new_x_axis_2;
		old_x_axis_2 = make_pair(0, 1);
		new_x_axis_2 = make_pair( (com_edge_idx_2+1)%4, com_edge_idx_2);

		zjucad::matrix::matrix<double> trans_mat_1 = GetTransMatrixInOneChart(from_chart_id,
			old_x_axis_1, new_x_axis_1);
		zjucad::matrix::matrix<double> trans_mat_2 = GetTransMatrixInOneChart(to_chart_id,
			old_x_axis_2, new_x_axis_2);
	   inv(trans_mat_2);

		return trans_mat_2 * trans_mat_1;
	}

	zjucad::matrix::matrix<double> QuadParam::GetTransMatrixInOneChart(int chart_id, pair<int,int> old_x_axis, pair<int,int> new_x_axis) const
	{
		zjucad::matrix::matrix<double> t_mat(zjucad::matrix::eye<double>(3));
		zjucad::matrix::matrix<double> r_mat(zjucad::matrix::eye<double>(3));

	    int old_origin = old_x_axis.first, new_origin = new_x_axis.first;
		if(old_origin != new_origin) 
		{
			double tx = m_chart_nodes_local_coord[chart_id][old_origin][0] - m_chart_nodes_local_coord[chart_id][new_origin][0];
			double ty = m_chart_nodes_local_coord[chart_id][old_origin][1] - m_chart_nodes_local_coord[chart_id][new_origin][1];
			t_mat(0, 2) = tx; t_mat(1, 2) = ty;
		}

		if(old_x_axis != new_x_axis)
		{
			Coord2D old_x_vec = Coord2D(m_chart_nodes_local_coord[chart_id][old_x_axis.second] -
				m_chart_nodes_local_coord[chart_id][old_x_axis.first]);
			Coord2D new_x_vec = Coord2D(m_chart_nodes_local_coord[chart_id][new_x_axis.second] -
				m_chart_nodes_local_coord[chart_id][new_x_axis.first]);

			double cross_v = old_x_vec[0] * new_x_vec[1] - old_x_vec[1] * new_x_vec[0];
			double dot_v = old_x_vec[0] * new_x_vec[0] + old_x_vec[1] * new_x_vec[1];

			double r_angle;
			if(fabs(dot_v) < SMALL_ZERO_EPSILON) 
			{
				if(cross_v < 0) r_angle = PI*3/2;
				else r_angle = PI/2;
			}else if(fabs(cross_v) < SMALL_ZERO_EPSILON)
			{
				if(dot_v < 0) r_angle = PI;
				else r_angle = 0;
			}
			double cos_v = cos(r_angle), sin_v = sin(r_angle);
			r_mat(0, 0) = cos_v; r_mat(0, 1) = sin_v;
			r_mat(1, 0) = -sin_v; r_mat(1, 1) = cos_v;
		}

		return r_mat*t_mat;
	}
	

	bool QuadParam::GetTranslistBetweenTwoCharts(int from_chart_id, int to_chart_id, vector<int>& trans_list) const
	{
		trans_list.clear();

		std::set<int> visited_face_set;
		std::map<int, int> prev_chart_map;

		queue<int> q;
		q.push(to_chart_id);
		visited_face_set.insert(to_chart_id);
		prev_chart_map[to_chart_id] = -1;


		bool flag = false;
		while(!q.empty())
		{
			int cur_chart_id = q.front(); q.pop();
			if(cur_chart_id == from_chart_id) { flag = true; break;}

			vector<int> adj_charts;
			for(size_t k=0; k<4; ++k)
			{
				const ChartNeighFun& func = m_chart_neigh_func[cur_chart_id*4 + k];
				assert(func.m_from_chart_id == cur_chart_id);
				if(func.m_to_chart_id != -1) adj_charts.push_back(func.m_to_chart_id);
			}

			for(size_t k=0; k<adj_charts.size(); ++k)
			{
				int chart_id = adj_charts[k];
				if(visited_face_set.find(chart_id) == visited_face_set.end())
				{
					q.push(chart_id);
					visited_face_set.insert(chart_id);
					prev_chart_map[chart_id] = cur_chart_id;
				}
			}
		}

		if(!flag) return false;

		trans_list.push_back(from_chart_id);
		int prev_chart = prev_chart_map[from_chart_id];

		while(prev_chart != -1)
		{
			trans_list.push_back(prev_chart);
			prev_chart = prev_chart_map[prev_chart];
		}

		return true;
	}

	int QuadParam::GetTransFuncBetweenTwoVertics(int vid1, int vid2, 
		TransMode tran_mode, ChartTransFun& tran_fun)
	{
		int chart_id_1 = m_vertex_group[vid1];
		int chart_id_2 = m_vertex_group[vid2];


		int func_idx(-1);
		for(size_t k=0; k<4; ++k)
		{
			const ChartNeighFun& func = m_chart_neigh_func[chart_id_1*4 + k];
			if(func.m_from_chart_id == chart_id_1 && func.m_to_chart_id == chart_id_2)
			{
				func_idx = chart_id_1*4 + k;
			}
		}


		if(func_idx == -1) 
		{
			printf("Transition between %d and %d\n", chart_id_1, chart_id_2);
			bool flag = false;
			vector<ChartNeighFun> chain_t_funs;
			for (size_t k = 0; k < 4; ++k)
			{
				const ChartNeighFun& func_1 = m_chart_neigh_func[chart_id_1 * 4 + k];

				for (size_t j = 0; j < 4; ++j)
				{
					if(func_1.m_to_chart_id < 0) continue;
				
					const ChartNeighFun& func_2 = m_chart_neigh_func[ func_1.m_to_chart_id*4 + j];
					if(func_2.m_to_chart_id == chart_id_2)
					{
						chain_t_funs.push_back(func_2);
						chain_t_funs.push_back(func_1);
						flag = true;
						break;
					}				
				}
				if(flag) break;
			}
			if(!flag) return -1;

			vector<ChartTransFun> trans_fun_vec;
			TransMode next_trans_mode = tran_mode;

			//
			for (size_t i = 0; i < chain_t_funs.size(); i++)
			{
				ChartNeighFun& n_fun = chain_t_funs[i];
				ChartTransFun t_fun;
				GetTransFunc(next_trans_mode, n_fun.m_type, t_fun);
			
				trans_fun_vec.push_back(t_fun);
				next_trans_mode = t_fun.m_mode == S_MODE ? TRANS_S_MODE : TRANS_T_MODE;
			}
			tran_fun = trans_fun_vec[0];
			for (size_t i = 1; i < trans_fun_vec.size(); i++) 
			{
				tran_fun = tran_fun + trans_fun_vec[i];
			}
		//	printf("%d to %d : %lf %lf\n", chart_id_1, chart_id_2, tran_fun.m_coef, tran_fun.m_div);
		}
		else
		{
			const ChartNeighFun& func = m_chart_neigh_func[func_idx];
			GetTransFunc(tran_mode, func.m_type, tran_fun);
		}

		return 0;
	}

	int QuadParam::TransSTBetweenTwoCharts(int chart_id_1, int chart_id_2, double& s_coord, double& t_coord) const
	{
		zjucad::matrix::matrix<double> origin_pos(3, 1);
		origin_pos(0, 0) = s_coord; origin_pos(1, 0) = t_coord; origin_pos(2, 0) = 1;

		zjucad::matrix::matrix<double> new_pos = GetTransMatrixBetweenTwoCharts(chart_id_1, chart_id_2) * origin_pos;

		s_coord = new_pos(0, 0);
		t_coord = new_pos(1, 0);

		return 0;

		int func_idx(-1);
		for(size_t k=0; k<4; ++k)
		{
			const ChartNeighFun& func = m_chart_neigh_func[chart_id_1*4 + k];
			if(func.m_from_chart_id == chart_id_1 && func.m_to_chart_id == chart_id_2)
			{
				func_idx = chart_id_1*4 + k;
			}
		}

		if(func_idx == -1) 
		{
			//printf("%lf %lf\n", s_coord, t_coord);
			bool flag = false;
			vector<ChartNeighFun> chain_t_funs;
			for (size_t i = 0; i < 4; ++i)
			{
				const ChartNeighFun& func_1 = m_chart_neigh_func[chart_id_1 * 4 + i];

				for (int j = 0; j < 4; j++)
				{
					if(func_1.m_to_chart_id == -1) continue;

					const ChartNeighFun& func_2 = m_chart_neigh_func[ func_1.m_to_chart_id*4 + j];
					if (func_2.m_to_chart_id == chart_id_2)
					{
						chain_t_funs.push_back(func_2);
						chain_t_funs.push_back(func_1);
						flag = true;
						break;
					}
				}
				if(flag) break;
			}
			if(!flag) 
				return -1;

			for (int i = (int)chain_t_funs.size()-1; i >= 0; i--) 
			{
				ChartNeighFun& fun = chain_t_funs[i];
				TransChartCoord(fun.m_type, s_coord, t_coord);
			}
			//printf("%lf %lf\n", s_coord, t_coord);
		}else
		{
			const ChartNeighFun& func = m_chart_neigh_func[func_idx];
			TransChartCoord(func.m_type, s_coord, t_coord);
		}
		return 0;
	}

	int QuadParam::GetTransFunc(TransMode mode, TransType type, ChartTransFun& transFun)
	{
		switch(type)
		{
		case UP_TO_UP:
			if(mode == TRANS_S_MODE)
			{
				transFun.m_coef = -1.0;
			  	transFun.m_mode = S_MODE;
		   		transFun.m_div = 1.0;
			}else if (mode == TRANS_T_MODE)	
			{
				transFun.m_coef = -1.0;
			   	transFun.m_mode = T_MODE;
				transFun.m_div = 2.0;
			}
			break;
		case UP_TO_DOWN:
			if (mode == TRANS_S_MODE)
			{
				transFun.m_coef = 1.0;					
				transFun.m_mode = S_MODE;
				transFun.m_div = 0;
			}else if (mode == TRANS_T_MODE)
			{
				transFun.m_coef = 1.0;
				transFun.m_mode = T_MODE;
				transFun.m_div = -1.0;
			}
			break;
		case UP_TO_LEFT:
			if (mode == TRANS_S_MODE)
			{
				transFun.m_coef = 1.0;
				transFun.m_mode = T_MODE;
				transFun.m_div = -1.0;
			}else if (mode == TRANS_T_MODE)
			{
				transFun.m_coef = -1.0;
				transFun.m_mode = S_MODE;
				transFun.m_div = 1.0;
			}
			break;
		 case UP_TO_RIGHT:
			 if (mode == TRANS_S_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = T_MODE;
				 transFun.m_div = 2.0;
			 }else if (mode == TRANS_T_MODE)
			 {
				 transFun.m_coef = 1.0;
				 transFun.m_mode = S_MODE;
				 transFun.m_div = 0;
			 }
			 break;
		 case DOWN_TO_DOWN:
			 if (mode == TRANS_S_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = S_MODE;
				 transFun.m_div = 1.0;
			 }else if (mode == TRANS_T_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = T_MODE;
				 transFun.m_div = 0;
			 }
			 break;
		 case DOWN_TO_UP:
			 if (mode == TRANS_S_MODE)
			 {
				 transFun.m_coef = 1.0;
				 transFun.m_mode = S_MODE;
				 transFun.m_div = 0;
			 }else if (mode == TRANS_T_MODE)
			 {
				transFun.m_coef = 1.0;
				transFun.m_mode = T_MODE;
				transFun.m_div = 1.0;
			 }
			 break;
		  case DOWN_TO_LEFT:
			 if (mode == TRANS_S_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = T_MODE;
				 transFun.m_div = 0;
			 }else if (mode == TRANS_T_MODE)
			 {
				 transFun.m_coef = 1.0;
				 transFun.m_mode = S_MODE;
				 transFun.m_div = 0;
			 }
			 break;
		
		  case DOWN_TO_RIGHT:
			 if (mode == TRANS_S_MODE)
			 {
				 transFun.m_coef = 1.0;
				 transFun.m_mode = T_MODE;
				 transFun.m_div = 1.0;
			 }else if (mode == TRANS_T_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = S_MODE;
				 transFun.m_div = 1.0;
			 }
			 break;
			
		  case LEFT_TO_LEFT:
			 if (mode == TRANS_S_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = S_MODE;
				 transFun.m_div = 0;
			 }else if (mode == TRANS_T_MODE)
			 {
				 transFun.m_coef = -1.0;
				 transFun.m_mode = T_MODE;
				 transFun.m_div = 1.0;
			 }
			 break;
			
		  case LEFT_TO_RIGHT:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = 1.0;
			  }else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 0;
			  }
			  break;
		  case LEFT_TO_UP:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = -1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 1.0;
			  }
			  else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = 1.0;
			  }
			  break;
		  case LEFT_TO_DOWN:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 0;
			  }
			  else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = -1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = 0;
			  }
			  break;
		  case RIGHT_TO_RIGHT:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = -1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = 2.0;
			  }
			  else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = -1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 1.0;
			  }
			  break;
		  case RIGHT_TO_LEFT:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = -1.0;
			  }else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 0;
			  }
			  break;
		  case RIGHT_TO_UP:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 0;
			  }
			  else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = -1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = 2.0;
			  }
			  break;
		  case RIGHT_TO_DOWN:
			  if (mode == TRANS_S_MODE)
			  {
				  transFun.m_coef = -1.0;
				  transFun.m_mode = T_MODE;
				  transFun.m_div = 1.0;
			  }
			  else if (mode == TRANS_T_MODE)
			  {
				  transFun.m_coef = 1.0;
				  transFun.m_mode = S_MODE;
				  transFun.m_div = -1.0;
			  }
			  break;
		  default:
			  break;
		  }
			
	
		  transFun.m_trans_type_vec.push_back(type);
		  return 0;
	}

	int QuadParam::TransChartCoord(TransType type, double& s_coord, double& t_coord) const
	{
		double origin_s_coord = s_coord;
		double origin_t_coord = t_coord;
		switch(type)
		{
		case UP_TO_UP:
			s_coord = 1 - origin_s_coord; /// 1-s
			t_coord = 2 - origin_t_coord; /// 2-t
			break;
		case UP_TO_DOWN:
			s_coord = origin_s_coord; /// s
			t_coord = -1 + origin_t_coord; /// t-1
			break;
		case UP_TO_LEFT:
			s_coord = -1 + origin_t_coord; /// t-1
			t_coord = 1 - origin_s_coord; /// 1-s
			break;
		case UP_TO_RIGHT:
			s_coord = 2 - origin_t_coord; /// 2-t
			t_coord = origin_s_coord; /// s
			break;

		case DOWN_TO_DOWN:
			s_coord = 1 - origin_s_coord; /// 1-s
			t_coord = 0 - origin_t_coord; /// -t
			break;
		case DOWN_TO_UP:
			s_coord = 0 + origin_s_coord; /// s
			t_coord = 1 + origin_t_coord; /// 1+t
			break;
		case DOWN_TO_LEFT:
			s_coord = 0 - origin_t_coord; /// -t
			t_coord = 0 + origin_s_coord; /// s
			break;
		case DOWN_TO_RIGHT:
			s_coord = 1 + origin_t_coord; /// 1+t
			t_coord = 1 - origin_s_coord; /// 1-s
			break;

		case LEFT_TO_LEFT:
			s_coord = 0 - origin_s_coord; /// -s
			t_coord = 1 - origin_t_coord; /// 1-t
			break;
		case LEFT_TO_RIGHT:
			s_coord = 1 + origin_s_coord; /// 1 + s
			t_coord = 0 + origin_t_coord; /// t
			break;
		case LEFT_TO_UP:
			s_coord = 1 - origin_t_coord; /// 1-t
			t_coord = 1 + origin_s_coord; /// 1+s
			break;
		case LEFT_TO_DOWN:
			s_coord = 0 + origin_t_coord; /// t
			t_coord = 1 - origin_s_coord; /// 1-s
			break;

		case RIGHT_TO_RIGHT:
			s_coord = 2 - origin_s_coord; /// 2-s
			t_coord = 1 - origin_t_coord; /// 1-t
			break;
		case RIGHT_TO_LEFT:
			s_coord = -1 + origin_s_coord; /// s-1
			t_coord = origin_t_coord; /// t
			break;
		case RIGHT_TO_UP:
			s_coord = origin_t_coord; /// t
			t_coord = 2 - origin_s_coord; /// 2-s
			break;
		case RIGHT_TO_DOWN:
			s_coord = 1 - origin_t_coord; /// 1-t
			t_coord = -1 + origin_s_coord; /// s-1
			break;

		default:
			break;
		}
		return 0;
	}

	int QuadParam::SolveParam()
	{
   		FormQuadChart();
		FormQuadHalfEdge();
		FormChartNeighFunc();

		SetCotCoef();
		SetLapMatrix();

		size_t loop_num = 5;
		for(size_t i=0; i<loop_num; ++i)
		{
			SolveQuadParam();

			if(i<loop_num-1)
			{
				AdjustVertexGroup();
			}
		}

		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		m_chart_param_coord.clear(); 
		m_chart_param_coord.resize(vert_num);
		for(int vid = 0; vid<vert_num; ++vid)
		{
			int chart_id = m_vertex_group[vid];
			m_chart_param_coord[vid] = ChartParamCoord(m_param_coord[vid], chart_id);
		}

		SetAllChartTexture();
		SaveSTCoord("st_coord.txt");

		//ComputeDistortion();
		ComputeAllChartDistortion();

		//FaceVaule2VtxColor(m_local_distortion);
		VtxValue2VtxColor(m_vertex_distortion);
		return 0;
	}

	int QuadParam::SolveQuadParam()
	{

		size_t vert_num = p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();
		int nb_variables_ = (int) vert_num*2;

		LinearSolver linear_solver(nb_variables_);

		/// fix conner and boundary vertex
		for(size_t k=0; k<m_chart_array.size(); ++k)
		{
			QuadChart& q_chart = m_chart_array[k];
			
			vector<int>& chart_halfedges = q_chart.m_he_index_array;
			assert(chart_halfedges.size() == 4);

			map< int, pair<double, double> > node_st;
			
			for(size_t k=0; k<4; ++k)
			{
				int he_idx = chart_halfedges[k];
				QuadHalfEdge& he = m_quad_he_array[he_idx];
			
				int node = he.m_node_id;
			
				double s_coord, t_coord;
				GetNodeSTCoord(he.m_type, s_coord, t_coord);

				node_st[ he.m_node_id] = make_pair(s_coord, t_coord);

				if(m_vertex_group[node] != q_chart.m_id) continue;

				int var_index_1 = he.m_node_id*2;
				int var_index_2 = he.m_node_id*2 + 1;
				linear_solver.variable( var_index_1).lock();
				linear_solver.variable( var_index_1).set_value(s_coord);
				linear_solver.variable( var_index_2).lock();
				linear_solver.variable( var_index_2).set_value(t_coord);

			}
				   
			for(size_t k=0; k<4; ++k)
			{
				int he_idx = chart_halfedges[k];
				QuadHalfEdge& he = m_quad_he_array[he_idx];
				/// fix boundary vertex 
				if(he.m_pair_id == -1)
				{
					int q_path_idx = he.m_path_id;
					QuadPath& quad_path = m_chart_path_array[q_path_idx];
				 					
					vector<int>& mesh_path = quad_path.m_path;
					assert(mesh_path.size() >=2);
					int start_node = quad_path.GetStartNode();
					int end_node = quad_path.GetEndNode();

					double start_s_coord = node_st[start_node].first;
					double start_t_coord = node_st[start_node].second;
					double end_s_coord = node_st[end_node].first;
					double end_t_coord = node_st[end_node].second;

					double path_len = GetPathLength(mesh_path, 0, mesh_path.size());
					for(size_t i=1; i<mesh_path.size()-1; ++i)
					{
						int mesh_vtx_idx = mesh_path[i];
						if(m_vertex_group[mesh_vtx_idx] != q_chart.m_id)
						{
							printf("Adjust %d from chart %d to chart %d!\n", 
								mesh_vtx_idx, m_vertex_group[mesh_vtx_idx], q_chart.m_id);
							m_vertex_group[mesh_vtx_idx] = q_chart.m_id;
						}
						double cur_len = GetPathLength(mesh_path, 0, i+1);
						double lambda = cur_len/path_len;
						double s_coord = (1-lambda)*start_s_coord + lambda*end_s_coord;
						double t_coord = (1-lambda)*start_t_coord + lambda*end_t_coord;

						int var_index_1 = mesh_path[i]*2;
						int var_index_2 = mesh_path[i]*2 + 1;
						linear_solver.variable( var_index_1).lock();
						linear_solver.variable( var_index_1).set_value(s_coord);
						linear_solver.variable( var_index_2).lock();
						linear_solver.variable( var_index_2).set_value(t_coord);
					}
				
				}
			}
		
		}

		linear_solver.begin_equation();

		for(size_t i=0; i<vert_num; ++i)
		{
			if(p_mesh->m_BasicOp.IsBoundaryVertex(i) /*|| IsConnerVertex(i)*/) continue;

			int to_vid = (int) i;
			int to_chart_id = m_vertex_group[to_vid];

			vector<int>& row_index = m_lap_matrix.m_RowIndex[i];
			vector<double>& row_data = m_lap_matrix.m_RowData[i];

			for(int st_index = 0; st_index < 2; ++st_index)
			{
				linear_solver.begin_row();
				
				double right_b_ = 0;
				for(size_t j=0; j<row_index.size(); ++j)
				{
					int from_vid = row_index[j];
					int from_chart_id = m_vertex_group[from_vid];
					
					if(m_vertex_group[from_vid] == m_vertex_group[to_vid]) /// transition function is identity
					{
						linear_solver.add_coefficient(from_vid*2 + st_index, row_data[j]);
					}else
					{
						ChartTransFun trans_func;
						TransMode trans_mode = (st_index == 0) ? TRANS_S_MODE : TRANS_T_MODE;
						if(GetTransFuncBetweenTwoVertics(from_vid, to_vid, trans_mode, trans_func) != 0)
						{
							printf("Can't find translation function between %d and %d!\n", from_vid, to_vid);
						}
						int var_index = from_vid*2 + ((trans_func.m_mode == S_MODE) ? 0 : 1);
						linear_solver.add_coefficient(var_index, row_data[j]*trans_func.m_coef);
						right_b_ -= (row_data[j]*trans_func.m_div);
					}
				}

				linear_solver.set_right_hand_side(right_b_);
				linear_solver.end_row();
			}

		}
		linear_solver.end_equation();
		
		linear_solver.solve();

		m_param_coord.clear();
		m_param_coord.resize(vert_num);

		for (size_t i = 0; i < vert_num; ++i)
		{
			int var_index = (int) i*2;
			double s_ = linear_solver.variable(var_index).value();
			double t_ = linear_solver.variable(var_index+1).value();

			m_param_coord[i].s_coord = s_;
			m_param_coord[i].t_coord = t_;
		}

		return 0;
	}

	void QuadParam::SetFaceGroup(std::vector<int>& face_group)
	{
		const PolyIndexArray& faceIndexArray = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		//
		size_t face_num = faceIndexArray.size();
		face_group.resize(face_num);
		fill(face_group.begin(), face_group.end(), -1);
		for (size_t i = 0; i < face_num; ++i)
		{
			const IndexArray& face_ = faceIndexArray[i];
			vector<int> chart_id_vec;
			for (size_t j = 0; j < 3; j++)
			{
				chart_id_vec.push_back(m_vertex_group[face_[j]]);
			}
			int max_times(-1), group_;
			for(size_t k=0; k<chart_id_vec.size(); ++k)
			{
				int g = chart_id_vec[k];
				int cur_times = count(chart_id_vec.begin(), chart_id_vec.end(), g);
				if(cur_times > max_times)
				{
					max_times = cur_times;
					group_ = g;
				}
			}
			face_group[i] = group_;
		}
		
	//	SetFaceGroupColor(face_group);
	}

	void QuadParam::SetAllChartTexture()
	{
		PolyTexCoordArray& faceTexCoord = p_mesh->m_Kernel.GetFaceInfo().GetTexCoord();

		//
		size_t face_num = p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();
		faceTexCoord.resize(face_num);

		m_face_tex_flag.clear();
		m_face_tex_flag.resize(face_num, false);
		
		//
		SetFaceGroup(m_face_tex_group);
		for (size_t i = 0 ; i < m_chart_array.size(); ++i)
		{
			SetChartTexture(i, m_face_tex_group);
		}

		m_face_group = m_face_tex_group;
	}

	void QuadParam::SetChartTexture(int chart_id, std::vector<int>& face_group)
	{
		PolyTexCoordArray& faceTexCoord = p_mesh->m_Kernel.GetFaceInfo().GetTexCoord();
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
 
		size_t i, k;
		double sValue, tValue;

	   size_t face_num = fIndex.size();
	   
		for (i = 0; i < face_num; ++i)
		{
			if (face_group[i] == chart_id && !m_face_tex_flag[i]) 
			{
				const IndexArray& face = fIndex[i];
				TexCoordArray& face_tex = faceTexCoord[i];
				face_tex.clear();

				for (k = 0; k < face.size(); k++)
				{
					const VertexID& vID = face[k];
					sValue = m_param_coord[vID].s_coord;
					tValue = m_param_coord[vID].t_coord;

					if (m_vertex_group[vID] != chart_id)
					{
						if(TransSTBetweenTwoCharts(m_vertex_group[vID], chart_id, sValue, tValue) !=0)
						{
							printf("%d Can't translate charts %d to %d!\n", vID, m_vertex_group[vID], chart_id);
						}
					}
					face_tex.push_back(Coord2D(sValue, tValue));
				}
				m_face_tex_flag[i] = true;
			}
		}
	}

	void QuadParam::GetNodeSTCoord(HalfEdgeType he_type, double& s_coord, double& t_coord)
	{
		switch(he_type)
		{
		case DOWN_HALF_EDGE:
			s_coord = 0.0;
			t_coord = 0.0;
			break;
			
		case RIGHT_HALF_EDGE:
			s_coord = 1.0;
			t_coord = 0.0;
			break;

		case UP_HALF_EDGE:
			s_coord = 1.0;
			t_coord = 1.0;
			break;

		case LEFT_HALF_EDGE:
			s_coord = 0.0;
			t_coord = 1.0;
			break;
		}
	}
	
	double QuadParam::GetPathLength(const std::vector<int>& _path, int from_idx, int end_idx)
	{
	 
		assert(from_idx < end_idx && from_idx >=0 && from_idx < (int)_path.size() -1&&
			end_idx >0 && end_idx <=(int) _path.size());

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		double len = 0.0;
		for(int k=from_idx+1; k!=end_idx; ++k)
		{
			int vtx1 = _path[k-1];
			int vtx2 = _path[k];
			len += (vCoord[vtx2] - vCoord[vtx1]).abs(); 
		}
		return len;
	}

	HalfEdgeType QuadParam::GetHalfEdgeType(int index)
	{
		HalfEdgeType type = UN_DEFINE_HALF_EDGE;
		switch(index)
		{
		case 0:
			type = DOWN_HALF_EDGE;
			break;
		case 1:
			type = RIGHT_HALF_EDGE;
			break;
		case 2:
			type = UP_HALF_EDGE;
			break;
		case 3:
			type = LEFT_HALF_EDGE;
			break;
		default:
			break;
		}
		return type;
	}

	void QuadParam::SetCotCoef()
	{
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		size_t nFace = fIndex.size();
		m_cot_coef_vec.clear();
		m_cot_coef_vec.resize(nFace);

		size_t i, j;
		Coord e[3];
		double a[3];
		for(i = 0; i < nFace; ++ i)
		{
			const IndexArray& f = fIndex[i];
			for(j = 0; j < 3; ++ j)
			{
				e[j] = (vCoord[f[(j+1)%3]] - vCoord[f[j]]).unit();
			}

			a[0] = angle(e[0], -e[2]);
			a[1] = angle(e[1], -e[0]);
			a[2] = PI-a[0]-a[1];

			for(j = 0; j < 3; ++ j)
			{
				// angle < 1 or angle > 179
				if(fabs(a[j]-0.0) < 0.0174 || fabs(a[j]-PI) < 0.0174)
				{
					m_cot_coef_vec[i][j] = 57.289;   // atan(1)
				}
				else
				{
					m_cot_coef_vec[i][j] = 1.0/tan(a[j]);
				}
			}
		}

		int vert_num = (int) p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();
		// Check the alpha+belta < PI is satisfied or not
		const PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		int nAdjust = 0;
		int k, h;
		for(i = 0; i < (size_t)vert_num ; ++i)
		{
			const IndexArray& adjf = vAdjFaces[i];
			size_t n = adjf.size();
			size_t begin_idx = 0, end_idx = n-1;
			if(p_mesh->m_BasicOp.IsBoundaryVertex((int) i))
			{
				begin_idx = 1;
				end_idx = n-1;
			}
			for(j = begin_idx; j <= end_idx; ++ j)
			{
				FaceID fID1 = adjf[(j+n-1)%n];
				FaceID fID2 = adjf[j];
				const IndexArray& f1 = fIndex[fID1];
				const IndexArray& f2 = fIndex[fID2];
				for( k = 0; k < 3; ++ k)
				{
					if(f1[k] == i)
						break;
				}
				for( h = 0; h < 3; ++ h)
				{
					if(f2[h] == i)
						break;
				}
				if(m_cot_coef_vec[fID1][(k+1)%3] + m_cot_coef_vec[fID2][(h+2)%3] < 0.0)
				{
					m_cot_coef_vec[fID1][(k+1)%3] = m_cot_coef_vec[fID2][(h+2)%3] = 0.0;
					++ nAdjust;
				}
			}
		}
		printf("#Cotangent Weight Adjust = %d\n", nAdjust);
	}

	void QuadParam::SetLapMatrix()
	{
		PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		int i, j, k, n;
		int fID, vID;
		int row, col;
		double coef, d;

		size_t vert_num = p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();
		//
		int nb_variables_ = (int) vert_num;
		m_lap_matrix.SetRowCol(nb_variables_, nb_variables_);

		for(i = 0; i < nb_variables_; ++ i)
		{
			IndexArray& adjFaces = vAdjFaces[i];
			row = i;
			n = (int) adjFaces.size();

			for (j = 0; j < n; j++)
			{
				fID = adjFaces[j];
				IndexArray& f = fIndex[fID];

				// Find the position of vertex i in face fID
				for(k = 0; k < 3; ++ k)
				{
					if(f[k] == i)
					{
						break;
					}
				}
				assert(k!=3);

				vID = f[(k+1)%3];
				col = vID;

				coef = m_cot_coef_vec[fID][(k+2)%3];
				m_lap_matrix.GetElement(row, col, d);
				d -= coef;
				m_lap_matrix.SetElement(row, col, d);

				vID = f[(k+2)%3];
				col = vID;

				coef = m_cot_coef_vec[fID][(k+1)%3];
				m_lap_matrix.GetElement(row, col, d);
				d -= coef;
				m_lap_matrix.SetElement(row, col, d);

				m_lap_matrix.GetElement(row, row, d);
				d += m_cot_coef_vec[fID][(k+1)%3] + m_cot_coef_vec[fID][(k+2)%3];
				m_lap_matrix.SetElement(row, row, d);
			}
		}

	}

	vector<int> QuadParam::GetEdgeAdjFaces(int vtx1, int vtx2) const
	{
		vector<int> adj_faces;
		if(p_mesh == NULL)
			return adj_faces;

		const PolyIndexArray& vf_adj_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const IndexArray& adj_faces1 = vf_adj_array[vtx1];
		const IndexArray& adj_faces2 = vf_adj_array[vtx2];

		for(size_t k=0; k<adj_faces1.size(); ++k)
		{
			int fid = adj_faces1[k];
			if(find(adj_faces2.begin(), adj_faces2.end(), fid) != adj_faces2.end())
			{
				adj_faces.push_back(fid);
			}
		}

		return adj_faces;
	}

	int QuadParam::GetChartIdFromQuadPath(const vector<int>& q_path_index_array)
	{
		int ret =-1;
		vector<int> tmp_path(q_path_index_array);
		sort(tmp_path.begin(), tmp_path.end());
		for(size_t k=0; k<m_chart_array.size(); ++k)
		{
			const QuadChart& q_chart = m_chart_array[k];
			vector<int> q_path = q_chart.m_path_index_array;
			sort(q_path.begin(), q_path.end());
			if(q_path == tmp_path) 
			{
				ret = k;
				break;
			}
		}
		if(ret == -1) { printf("Can't find a chart!\n"); }
		return ret;
	}

	

	int QuadParam::AdjustVertexGroup()
	{
		vector<int> out_range_vtx_array;
		GetOutRangeVertex(out_range_vtx_array);
		printf("There are %d out range vertices!\n", out_range_vtx_array.size());

		int adjust_vtx_num = 0;
		for(size_t k=0; k<out_range_vtx_array.size(); ++k)
		{
			int out_range_vtx = out_range_vtx_array[k];
			int cur_chart_id = m_vertex_group[out_range_vtx];
			assert(cur_chart_id != -1);

			vector<int> nb_chart_array;
			FindNeighborCharts2Ring(cur_chart_id, nb_chart_array);

			double s_coord = m_param_coord[out_range_vtx].s_coord;
			double t_coord = m_param_coord[out_range_vtx].t_coord;
			bool flag = false;
			for(size_t i=0; i<nb_chart_array.size(); ++i)
			{
				int to_chart_id = nb_chart_array[i];
				double tmp_s_coord = s_coord;
				double tmp_t_coord = t_coord;

				TransSTBetweenTwoCharts(cur_chart_id, to_chart_id, tmp_s_coord, tmp_t_coord);

				if(InRangeST(tmp_s_coord, tmp_t_coord))
				{
					m_param_coord[out_range_vtx].s_coord = tmp_s_coord;
					m_param_coord[out_range_vtx].t_coord = tmp_t_coord;
					m_vertex_group[out_range_vtx] = to_chart_id;
					flag = true;
					adjust_vtx_num ++;
					break;
				}
			}
		}

		printf("Direct adjust %d vertex!\n", adjust_vtx_num);
		GetOutRangeVertex(out_range_vtx_array);
		printf("After adjust, there still are %d out range vertices!\n", out_range_vtx_array.size());
		return 0;
	}

	void QuadParam::GetOutRangeVertex(std::vector<int>& _out_range_vtx_array)
	{
		_out_range_vtx_array.clear();
		for(size_t k=0; k<m_param_coord.size(); ++k)
		{
			double s_coord = m_param_coord[k].s_coord;
			double t_coord = m_param_coord[k].t_coord;
			if(!InRangeST(s_coord, t_coord))
			{
				_out_range_vtx_array.push_back(k);
			}
		}
	}

	void QuadParam::FindNeighborCharts2Ring(int chart_id, std::vector<int>& chart_array)
	{
		const QuadChart& chart = m_chart_array[chart_id];

		chart_array.clear();
		for(size_t k=0; k<4; ++k)
		{
			const ChartNeighFun& chart_fun = m_chart_neigh_func[chart_id*4 + k];
			assert(chart_id == chart_fun.m_from_chart_id);
			if(chart_fun.m_to_chart_id == -1) continue;
			chart_array.push_back(chart_fun.m_to_chart_id);

			int adj_chart_id = chart_fun.m_to_chart_id;
			for(size_t i=0; i<4; ++i)
			{
				const ChartNeighFun& adj_chart_fun = m_chart_neigh_func[ adj_chart_id*4 + i];
				assert(adj_chart_id == adj_chart_fun.m_from_chart_id);
				int to_chart_id = adj_chart_fun.m_to_chart_id;
				if(to_chart_id == chart_id || to_chart_id == -1) continue;
				if(find(chart_array.begin(), chart_array.end(), to_chart_id) == chart_array.end())
				{
					chart_array.push_back(to_chart_id);
				}
			}
		}
		
	}

	void QuadParam::ComputeDistortion()
	{
		size_t face_num = p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();
		m_local_distortion.clear(); 
		m_local_distortion.resize(face_num);

		for(size_t k=0; k<m_chart_array.size(); ++k)
		{
			ComputeChartDistortion(k, m_face_group);
		}

		FaceVaule2VtxColor(m_local_distortion);
		const PolyIndexArray& vtx_adj_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		size_t vtx_num = vtx_adj_array.size();
		m_vertex_distortion.clear();
		m_vertex_distortion.resize(vtx_num);
		for(size_t k=0; k<vtx_num; ++k)
		{
			m_vertex_distortion[k] = 0.0;
			const IndexArray& adj_faces = vtx_adj_array[k];
			for(size_t i=0; i<adj_faces.size(); ++i)
			{
				int fid = adj_faces[i];
				m_vertex_distortion[k] += m_local_distortion[fid];
			}
			m_vertex_distortion[k] /= (adj_faces.size());
		}
	}

	void QuadParam::ComputeAllChartDistortion()
	{
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const DoubleArray& face_area_array = p_mesh->m_Kernel.GetFaceInfo().GetFaceArea();
		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		size_t face_num = face_index_array.size();
		m_local_distortion.clear(); 
		m_local_distortion.resize(face_num);

		for(size_t k=0; k<face_num; ++k)
		{
			int chart_id = m_face_tex_group[k];
				
			const IndexArray& face_vtx = face_index_array[k];
			double area = face_area_array[k];
			double distort = 0.0;

			vector<Coord> vtx_coord_3d (3);
			for(int i=0; i<3; ++i) vtx_coord_3d[i] = vtx_coord_array[face_vtx[i]];
			vector<Coord2D> vtx_coord_2d = ComputeTriangleLocal2DCoord(vtx_coord_3d);

			for(size_t i=0; i<3; ++i)				
			{
				int v_i = face_vtx[i]; //  index of vertex i 
				int v_i_1 = face_vtx[ (i+1) %3]; // index of vertex i+1
				int v_i_2 = face_vtx[ (i+2) %3]; // index of vertex i+2

				double s_coord_1 (m_param_coord[v_i_1].s_coord), t_coord_1(m_param_coord[v_i_1].t_coord);
				double s_coord_2 (m_param_coord[v_i_2].s_coord), t_coord_2(m_param_coord[v_i_2].t_coord);
				   
				if(m_vertex_group[v_i_1] != chart_id)
			  	{
					TransSTBetweenTwoCharts(m_vertex_group[v_i_1], chart_id, s_coord_1, t_coord_1);
				}

				if(m_vertex_group[v_i_2] != chart_id)
				{
					TransSTBetweenTwoCharts(m_vertex_group[v_i_2], chart_id, s_coord_2, t_coord_2);
				}

				
				
				Coord2D uv ( s_coord_2 - s_coord_1, t_coord_2 - t_coord_1);

				int idx_1 = (i+1)%3, idx_2 = (i+2)%3;

				Coord2D p1 = (vtx_coord_2d[idx_1]-vtx_coord_2d[i]);
				Coord2D p2 = (vtx_coord_2d[idx_2]-vtx_coord_2d[i]);

				double uv_2 = (uv.abs()) *(uv.abs());
				double p = dot(p1, p2);
				
				distort += uv_2*p;
			}
			distort /= (2*area*area);
			m_local_distortion[k] = distort;
		}

		vector<double>::iterator min_it = min_element(m_local_distortion.begin(), m_local_distortion.end());
		vector<double>::iterator max_it = max_element(m_local_distortion.begin(), m_local_distortion.end());

		printf("The min distortion is %lf\n", *min_it);
		printf("The max distortion is %lf\n", *max_it);

//		FaceVaule2VtxColor(m_local_distortion);
		const PolyIndexArray& vtx_adj_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		size_t vtx_num = vtx_adj_array.size();
		m_vertex_distortion.clear();
		m_vertex_distortion.resize(vtx_num);
		for(size_t k=0; k<vtx_num; ++k)
		{
			m_vertex_distortion[k] = 0.0;
			const IndexArray& adj_faces = vtx_adj_array[k];
			for(size_t i=0; i<adj_faces.size(); ++i)
			{
				int fid = adj_faces[i];
				m_vertex_distortion[k] += m_local_distortion[fid];
			}
			m_vertex_distortion[k] /= (adj_faces.size());
		}
	}

	void QuadParam::ComputeChartDistortion(int chart_id, const std::vector<int>& face_group)
	{
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const DoubleArray& face_area_array = p_mesh->m_Kernel.GetFaceInfo().GetFaceArea();
		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		size_t face_num = face_index_array.size();

		for(size_t k=0; k<face_num; ++k)
		{
			if(face_group[k] == chart_id)
			{
				const IndexArray& face_vtx = face_index_array[k];
				double area = face_area_array[k];
				double distort = 0.0;
				for(size_t i=0; i<3; ++i)
				{
					int v_i = face_vtx[i]; //  index of vertex i 
					int v_i_1 = face_vtx[ (i+1) %3]; // index of vertex i+1
					int v_i_2 = face_vtx[ (i+2) %3]; // index of vertex i+2

					double s_coord_1 (m_param_coord[v_i_1].s_coord), t_coord_1(m_param_coord[v_i_1].t_coord);
					double s_coord_2 (m_param_coord[v_i_2].s_coord), t_coord_2(m_param_coord[v_i_2].t_coord);
					if(m_vertex_group[v_i_1] != chart_id)
					{
						TransSTBetweenTwoCharts(m_vertex_group[v_i_1], chart_id, s_coord_1, t_coord_1);
					}

					if(m_vertex_group[v_i_2] != chart_id)
					{
						TransSTBetweenTwoCharts(m_vertex_group[v_i_2], chart_id, s_coord_2, t_coord_2);
					}
				
					Coord2D uv ( s_coord_2 - s_coord_1, t_coord_2 - t_coord_1);
					Coord p1 = (vtx_coord_array[v_i_1] - vtx_coord_array[v_i]);
					Coord p2 = (vtx_coord_array[v_i_2] - vtx_coord_array[v_i]);
					double uv_2 = (uv.abs()) *(uv.abs());
					double p = dot(p1, p2);
					distort += uv_2*p;
				}
				distort /= (2*area*area);
				m_local_distortion[k] = distort;
			}
		}

	}
    int QuadParam::SaveSTCoord(const std::string& file_name)
	{
		ofstream fout(file_name.c_str());
		if(fout.fail())
		{
			printf("save st_coord file fail!\n");
			return -1;
		}

		for(size_t k=0; k<m_param_coord.size(); ++k)
		{
			fout << m_param_coord[k].s_coord << " " << m_param_coord[k].t_coord << endl;
		}

		fout.close();

		return 0;
	}

	void QuadParam::SetFaceGroupColor(const std::vector<int>& face_group)
	{
		
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

		for(size_t i = 0; i < m_face_group.size(); ++ i)
		{
			int index = face_group[i] % 48;
			colorArray[i] = Color(colors[index][0], colors[index][1], colors[index][2]);
		}
	}

	void QuadParam::Draw() const
	{		
		glEnable(GL_LIGHTING);
		glEnable(GL_POLYGON_SMOOTH);
		glPolygonMode(GL_FRONT, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(2.0, 2.0);
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

		DrawPatchConner();

		DrawUnCorrespondingVertices();

		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_POLYGON_SMOOTH);

		glDisable(GL_LIGHTING);
		glDepthFunc(GL_LEQUAL);


		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glLineWidth(3.0f);

		DrawPatchEdge();

		glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glEnable(GL_LIGHTING);
		glDepthFunc(GL_LESS);
		glEnable(GL_LIGHTING);
	}
	
	void QuadParam::DrawPatchConner() const
	{
		if(p_mesh == NULL) return ;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		glColor3ub(0, 64, 128);
		for(size_t k=0; k<m_chart_node_array.size(); ++k)
		{
			int mesh_idx = m_chart_node_array[k].m_mesh_index;
			DrawSphere(vCoord[mesh_idx]);
		}
		//DrawSphere(vCoord[]);
	}

	void QuadParam::DrawPatchEdge() const
	{
		if(p_mesh == NULL) return;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		
		glColor3ub(255, 255, 0);
		for(size_t k=0; k<m_chart_path_array.size(); ++k)
		{
			const vector<int>& path = m_chart_path_array[k].m_path;
			glBegin(GL_LINE_STRIP);
			for(size_t i=0; i<path.size(); ++i)
			{
				const Coord& vtxCoord = vCoord[path[i]];
				glVertex3d(vtxCoord[0], vtxCoord[1], vtxCoord[2]);
			}
			glEnd();
		}
	}

	void QuadParam::DrawUnCorrespondingVertices() const
	{
		if(p_mesh == NULL) return;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		glColor3ub(255, 0, 0);
		for(size_t k=0; k<m_uncorresponding_vtx_array.size(); ++k)
		{
			int mesh_idx = m_uncorresponding_vtx_array[k];
			DrawSphere(vCoord[mesh_idx]);
		}
	}

	void QuadParam::DrawSphere(const Coord& center, double point_size) const
	{
		double bRadius;
		Coord c;
		p_mesh->m_Kernel.GetModelInfo().GetBoundingSphere(c, bRadius);
		GLUquadricObj* qobj = gluNewQuadric();
		gluQuadricDrawStyle(qobj, GLU_FILL);
		gluQuadricNormals(qobj, GLU_SMOOTH);
		glTranslatef((float)center[0], (float)center[1], (float)center[2]);

		gluSphere(qobj, point_size *0.02 * bRadius , 15, 10);
		glTranslatef((float)-center[0], (float)-center[1], (float)-center[2]);
	}

	void QuadParam::DrawDistortion() const
	{
	}

	void QuadParam::VtxValue2VtxColor(const vector<double>& vertexValue)
	{
		ColorArray& colorArray = p_mesh->m_Kernel.GetVertexInfo().GetColor();
		colorArray.clear();

		CHSVColor color;
		color.RGBtoHSV(1, 0, 0);
		color.m_S = 0.9f;
		color.m_V = 0.9f;

		size_t vNum = p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();

		double min = 1e20, max = -1e20;
		for(size_t i = 0; i < vNum; ++i)
		{
			if(vertexValue[i] < min)	min = vertexValue[i];
			if(vertexValue[i] > max)	max = vertexValue[i];
		}

		double range = (max - min) * 1.1;

		bool eq = ALMOST_EQUAL_LARGE(range, 0.0);

		for(size_t i = 0; i < vertexValue.size();++i)
		{
			float R, G, B;
			if (eq)
			{
				color.m_H = (float) 0.5 * 255;
			}
			else
			{
				double prop = (vertexValue[i] - min) / range;

				//color.m_S = prop;
				color.m_H = (float)(1 - prop) * 255;
			}
			color.HSVtoRGB(&R, &G, &B);
			Color c(R, G, B);
			colorArray.push_back(c);
		}
	}


	void QuadParam::FaceVaule2VtxColor(const std::vector<double>& face_value)
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

	void QuadParam::SetMeshTexCoordForGrid16()
	{
		const CoordArray& vtx_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		
		TexCoordArray& tex_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetTexCoord();
		tex_coord_array.clear();
		tex_coord_array.resize(vert_num);

		if(vert_num < 255) return;

		vector<int> conner_vtx;
		conner_vtx.push_back(0); conner_vtx.push_back(15);
		conner_vtx.push_back(240); conner_vtx.push_back(255);
// 		conner_vtx.push_back(0); conner_vtx.push_back(3);
// 		conner_vtx.push_back(12); conner_vtx.push_back(15);

		tex_coord_array[conner_vtx[0]] = Coord2D(1, 1);
		tex_coord_array[conner_vtx[1]] = Coord2D(-1, 1);
		tex_coord_array[conner_vtx[2]] = Coord2D(-1, -1);
		tex_coord_array[conner_vtx[3]] = Coord2D(1, -1);

		for(int k=0; k<vert_num; ++k)
		{
			if(find(conner_vtx.begin(), conner_vtx.end(), k) == conner_vtx.end())
			{
				double x = vtx_coord_array[k][0];//2 + 0.5;
				double y = vtx_coord_array[k][2];//2 + 0.5;
				tex_coord_array[k] = Coord2D(x, y);
			}
		}	   
	}

	void QuadParam::SetMeshTexCoord(const std::vector<TexCoord>& _tex_coord_array)
	{
		TexCoordArray& tex_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetTexCoord();
	    tex_coord_array = _tex_coord_array;
	}

	bool QuadParam::ComputeSurfaceCoord(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const
	{
		std::set<int> chart_set;
		queue<int> q;

		int chart_id = chart_param_coord.chart_id;
		q.push(chart_id);
		chart_set.insert(chart_id);

		while(!q.empty())
		{
			int cur_chart_id = q.front(); q.pop();
		    
			if(ComputeSurfaceCoordOnChart(chart_id, chart_param_coord, surface_coord)) return true;

			const QuadChart& quad_chart = m_chart_array[cur_chart_id];
			const vector<int> neighbor_charts = quad_chart.m_neighbor_chart_array;

			for(size_t k=0; k<neighbor_charts.size(); ++k)
			{
				if(chart_set.find(neighbor_charts[k]) == chart_set.end()) 
				{
					q.push(neighbor_charts[k]);
					chart_set.insert(neighbor_charts[k]);
				}
			}
			                                                                                                                        
		}

		return false;
	}

	bool QuadParam::ComputeSurfaceCoordOnChart(int chart_id, const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const
	{
		std::set<int> visited_face_set;

		const PolyIndexArray& vert_adjface_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const PolyIndexArray& face_index_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		
		const vector<int>& chart_vertices = m_chart_array[chart_id];

		for(size_t k=0; k<chart_vertices.size(); ++k)
		{
			int vid = chart_vertices[k];
			const vector<int>& nb_faces = vert_adjface_array[vid];
			for(size_t i=0; i<nb_faces.size(); ++i)
			{
				int fid = nb_faces[i];
				if(visited_face_set.find(fid) != visited_face_set.end()) continue;
				visited_face_set.insert(fid);
				const IndexArray& face_vertices = face_index_array[fid];
				vector<ParamCoord> node_param_coord(face_vertices.size());
				for(int j=0; j<3; ++j)
				{
					int cur_vid = face_vertices[j];
					node_param_coord[j] = m_param_coord[cur_vid];
					int cur_chart_id = m_vertex_group[cur_vid];
					assert(cur_chart_id == chart_id);
					if(cur_chart_id != chart_param_coord.chart_id)
					{
						TransParamCoordBetweenTwoChart(cur_chart_id, chart_param_coord.chart_id, 
							m_param_coord[cur_vid], node_param_coord[j]);
					}
				}

				Barycentrc baryc_coord = ComputeVertexBarycentric(node_param_coord, chart_param_coord.param_coord);
				if(IsValidBarycentic(baryc_coord))
				{
					surface_coord = SurfaceCoord(fid, baryc_coord);
					return true;
				}

			}

		}

		return false;
	}

} // end namespace 