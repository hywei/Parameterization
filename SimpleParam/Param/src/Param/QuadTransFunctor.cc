#include "QuadTransFunctor.h"
#include "QuadChartCreator.h"
#include "../hj_3rd/include/math/blas_lapack.h"
#include "../hj_3rd/include/zjucad/matrix/lapack.h"

#include <map>
#include <set>
#include <queue>
#include <iostream>


namespace PARAM
{
	QuadTransFunctor::QuadTransFunctor(boost::shared_ptr<QuadChartCreator> quad_chart_creator) 
		: p_quad_chart_creator(quad_chart_creator) {}
	QuadTransFunctor::~QuadTransFunctor(){}

	zjucad::matrix::matrix<double> QuadTransFunctor::GetTransMatrix(int from_chart_id, int to_chart_id) const
	{
		zjucad::matrix::matrix<double> trans_mat(zjucad::matrix::eye<double>(3));
		if(from_chart_id == to_chart_id ) return trans_mat;
		std::vector<int> trans_list;
		if(!GetTranslistBetweenTwoCharts(from_chart_id, to_chart_id, trans_list))
		{
			std::cout<< "Cannot get the transition list!\n";
			return trans_mat;
		}
		assert(trans_list.size() >=2);

		for(size_t k=1; k<trans_list.size(); ++k)
		{	
			zjucad::matrix::matrix<double> temp_trans_mat = trans_mat;
			trans_mat =  GetTransMatrixOfAdjCharts(trans_list[k-1], trans_list[k]) * temp_trans_mat;
		}

		return trans_mat;
	}

	bool QuadTransFunctor::GetTranslistBetweenTwoCharts(int from_chart_id, int to_chart_id, std::vector<int>& trans_list) const
	{
		trans_list.clear();

		const std::vector<QuadPatch>& patch_array = p_quad_chart_creator->GetQuadPatchArray();

		std::set<int> visited_face_set;
		std::map<int, int> prev_chart_map;

		std::queue<int> q;
		q.push(to_chart_id);
		visited_face_set.insert(to_chart_id);
		prev_chart_map[to_chart_id] = -1;


		bool flag = false;
		while(!q.empty())
		{
			int cur_chart_id = q.front(); q.pop();
			if(cur_chart_id == from_chart_id) { flag = true; break;}

			const std::vector<int>& adj_charts = patch_array[cur_chart_id].m_neighbor_patch_array;
			
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

	zjucad::matrix::matrix<double> QuadTransFunctor::GetTransMatrixInOneChart(int chart_id, std::pair<int,int> old_x_axis,
		std::pair<int,int> new_x_axis) const
	{
		const std::vector<QuadChart>& chart_array = p_quad_chart_creator->GetQuadChartArray();
		const QuadChart& quad_chart = chart_array[chart_id];
		const std::vector<ParamCoord> conner_param_coord_vec = quad_chart.m_conner_param_coord_array;

		zjucad::matrix::matrix<double> t_mat(zjucad::matrix::eye<double>(3));
		zjucad::matrix::matrix<double> r_mat(zjucad::matrix::eye<double>(3));

		int old_origin = old_x_axis.first, new_origin = new_x_axis.first;
		if(old_origin != new_origin) 
		{
			double tx = conner_param_coord_vec[old_origin].s_coord - conner_param_coord_vec[new_origin].s_coord;
			double ty = conner_param_coord_vec[old_origin].t_coord - conner_param_coord_vec[new_origin].t_coord;
			t_mat(0, 2) = tx; t_mat(1, 2) = ty;
		}

		if(old_x_axis != new_x_axis)
		{
			double x1, y1, x2, y2;
			x1 = conner_param_coord_vec[old_x_axis.second].s_coord - conner_param_coord_vec[old_x_axis.first].s_coord;
			y1 = conner_param_coord_vec[old_x_axis.second].t_coord - conner_param_coord_vec[old_x_axis.first].t_coord;
			x2 = conner_param_coord_vec[new_x_axis.second].s_coord - conner_param_coord_vec[new_x_axis.first].s_coord;
			y2 = conner_param_coord_vec[new_x_axis.second].t_coord - conner_param_coord_vec[new_x_axis.first].t_coord;

			double cross_v = x1*y2 - y1*x2;
			double dot_v = x1*x2 + y1*y2;

			double r_angle;
			if(fabs(dot_v) < LARGE_ZERO_EPSILON) 
			{
				if(cross_v < 0) r_angle = PI*3/2;
				else r_angle = PI/2;
			}else if(fabs(cross_v) < LARGE_ZERO_EPSILON)
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

	zjucad::matrix::matrix<double> QuadTransFunctor::GetTransMatrixOfAdjCharts(int from_chart_id, int to_chart_id) const
	{
		const std::vector<QuadPatch>& patch_array = p_quad_chart_creator->GetQuadPatchArray();
		const QuadPatch& from_patch = patch_array[from_chart_id];
		const QuadPatch& to_patch = patch_array[to_chart_id];

		const std::vector<PatchEdge>& patch_edge_array = p_quad_chart_creator->GetPatchEdgeArray();

		/// find the common edge of these two charts
		int com_edge_idx_1(-1), com_edge_idx_2(-1);
		for(size_t k=0; k<from_patch.m_edge_index_array.size(); ++k)
		{
			int edge_idx_1 = from_patch.m_edge_index_array[k];
			const PatchEdge& patch_edge = patch_edge_array[edge_idx_1];
			const std::vector<int>& pe_neighbors = patch_edge.m_neighbor_patch_array;
			if(find(pe_neighbors.begin(), pe_neighbors.end(), to_chart_id) != pe_neighbors.end())
			{
				for(size_t i=0; i<to_patch.m_edge_index_array.size(); ++i)
				{
					int edge_idx_2 = to_patch.m_edge_index_array[i];
					if(edge_idx_1 == edge_idx_2)
					{
						com_edge_idx_1 = k;
						com_edge_idx_2 = i;
						break;
					}
				}
			}
			if(com_edge_idx_1 != -1 && com_edge_idx_2 != -1) break;
		}

		std::pair<int, int> old_x_axis_1, new_x_axis_1;
		/// find the from chart and to chart's origin and x-axis
		old_x_axis_1 = std::make_pair(0, 1); 
		new_x_axis_1 = std::make_pair(com_edge_idx_1, (com_edge_idx_1+1)%4);

		std::pair<int, int> old_x_axis_2, new_x_axis_2;
		old_x_axis_2 = std::make_pair(0, 1);
		new_x_axis_2 = std::make_pair( (com_edge_idx_2+1)%4, com_edge_idx_2);

		zjucad::matrix::matrix<double> trans_mat_1 = GetTransMatrixInOneChart(from_chart_id,
			old_x_axis_1, new_x_axis_1);
		zjucad::matrix::matrix<double> trans_mat_2 = GetTransMatrixInOneChart(to_chart_id,
			old_x_axis_2, new_x_axis_2);
		inv(trans_mat_2);

		return trans_mat_2 * trans_mat_1;
	}
} 