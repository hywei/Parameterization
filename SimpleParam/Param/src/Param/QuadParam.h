#ifndef HYWEI_QUADPARAM_H_
#define HYWEI_QUADPARAM_H_

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

#include "QuadChart.h"
#include "../ModelMesh/MeshModel.h"
#include "../Numerical/MeshSparseMatrix.h"
#include <boost/shared_ptr.hpp>
#include "../hj_3rd/include/zjucad/matrix/matrix.h"

namespace PARAM
{
	class QuadParam
	{

	public:
		QuadParam();
		~QuadParam();

	public:
		int AttchMesh(const boost::shared_ptr<MeshModel> mesh);
		int LoadQuadFile(const std::string& quad_file_name);
	
		int SolveParam();

		void Draw() const;

	public:
		//! Translate parameter coordinate between two charts
		/*
		 \param from_chart_id the from chart's id
		 \param to_chart_id the to chart's id
		 \param original_coord the original parameter coordinate
		 \param new_coord the new parameter coordinate
		 */
		void TransParamCoordBetweenTwoChart(int from_chart_id, int to_chart_id, 
			const ParamCoord& original_coord, ParamCoord& new_coord) const;

		bool IsConnerVertex(int vid) const;

		//! Get transition matrix between two charts
		/*
		 \param from_chart_id the from chart's id
		 \param to_chart_id the to chart's id
		 \return the transition matrix between these two charts
		 */
		zjucad::matrix::matrix<double> GetTransMatrixBetweenTwoCharts(int from_chart_id,
			int to_chart_id) const;


		bool ComputeSurfaceCoord(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const;
		bool ComputeSurfaceCoordOnChart(int chart_id, const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const;	

	public:
		/// Get Method
		const std::vector<ParamCoord>& GetParamCoord() const
		{
			return m_param_coord;
		}

		PARAM::ParamCoord GetParamCoord(int vtx_index) const 
		{
			return m_param_coord[vtx_index];
		}

		const boost::shared_ptr<MeshModel>& GetMeshModel() const
		{
			return p_mesh;
		}

		int GetFaceChartID(int fid) const
		{
			return m_face_group[fid];
		}

		int GetVertexChartID(int vid) const
		{
			return m_vertex_group[vid];
		}

		const CMeshSparseMatrix& GetLaplaceMatrix() const
		{
			return m_lap_matrix;
		}


		void SetUnCorrespondingVtxArray(const std::vector<int>& uncorresponding_vtx_array)
		{
			m_uncorresponding_vtx_array = uncorresponding_vtx_array;
		}

		void SetVtxParamCoord(int vtx, const ParamCoord& param_coord)
		{
			m_param_coord[vtx] = param_coord;
		}


	public:
		int FormQuadChart();
		int FormQuadHalfEdge();
		int FormChartNeighFunc();

	    void SetChartNeighbors();

		int FloodFillChart(int face_id, std::set<int>& face_set, std::set<int>& quad_path_array,
			std::vector<bool>& face_visited_flag, 
			const std::map< std::pair<int,int>, std::vector<int> >& edge_pid_mapping);

		void SetVertexGroup();
		
		int FormChartHalfEdge(int chart_id);

		void SetCotCoef();
		void SetLapMatrix();

		int SolveQuadParam();
		

		int AdjustVertexGroup();
		void GetOutRangeVertex(std::vector<int>& _out_range_vtx_array);
		void FindNeighborCharts2Ring(int chart_id, std::vector<int>& chart_array);
		bool InRangeST(double s, double t)
		{
			return GreaterEqual(s, 0.0, SMALL_ZERO_EPSILON) && LessEqual(s, 1.0, SMALL_ZERO_EPSILON) && 
				GreaterEqual(t, 0.0, SMALL_ZERO_EPSILON) && LessEqual(t, 1.0, SMALL_ZERO_EPSILON);
		}


		/// help function
		//! get an edge's adjcent  faces
		std::vector<int> GetEdgeAdjFaces(int vtx1, int vtx2) const;
        int GetChartIdFromQuadPath(const vector<int>& q_path_index_array);

		//! set each face's color for view and debug
		void SetFaceGroupColor(const std::vector<int>& face_group);

		void GetNodeSTCoord(HalfEdgeType he_type, double& s_coord, double& t_coord);
		double GetPathLength(const std::vector<int>& _path, int from_idx, int end_idx);

		HalfEdgeType GetHalfEdgeType(int idx);
	
		int GetTransFuncBetweenTwoVertics(int vid1, int vid2, TransMode tran_mode, ChartTransFun& tran_fun);
		int GetTransFunc(TransMode mode, TransType type, ChartTransFun& transFun);

		int TransSTBetweenTwoCharts(int chart_id1, int chart_id2, double& s_coord, double& t_coord) const;
		int TransChartCoord(TransType type, double& s_coord, double& t_coord) const;
	 

		bool GetTranslistBetweenTwoCharts(int from_chart_id, int to_chart_id, std::vector<int>& trans_list) const;

		zjucad::matrix::matrix<double> GetTransMatrixOfAdjCharts(int from_chart, int to_chart) const;
		//! Get the transition matrix in one chart
		zjucad::matrix::matrix<double> GetTransMatrixInOneChart(int chart_id, 
			pair<int,int> old_x_axis, pair<int,int> new_x_axis) const;

		void SetFaceGroup(std::vector<int>& group);
		void SetAllChartTexture();
		void SetChartTexture(int chart_id, std::vector<int>& face_group);

		int SaveSTCoord(const std::string& file_name);


		void ComputeDistortion();
		void ComputeAllChartDistortion();
		void ComputeChartDistortion(int chart_id, const std::vector<int>& face_group);

		void DrawPatchConner() const;
		void DrawPatchEdge() const;
		void DrawSphere(const Coord& center, double point_size=1.0) const;
		void DrawDistortion() const;
		void DrawUnCorrespondingVertices() const;
		
		

    public:
		void VtxValue2VtxColor(const std::vector<double>& vertexValue);
		void FaceVaule2VtxColor(const std::vector<double>& face_value);

		/// just for model grid_16x16.obj
		void SetMeshTexCoordForGrid16();
		void SetMeshTexCoord(const std::vector<TexCoord>& tex_coord_array);
	  
	private:
		boost::shared_ptr<MeshModel> p_mesh;

		std::vector<QuadNode> m_chart_node_array;
		std::vector<QuadPath> m_chart_path_array;
		std::vector<QuadChart> m_chart_array;

		std::vector<QuadHalfEdge> m_quad_he_array;

		std::vector<int> m_vertex_group;
		std::vector<int> m_face_group;

		CMeshSparseMatrix m_lap_matrix;
		std::vector<Coord> m_cot_coef_vec;

		std::vector<ChartNeighFun> m_chart_neigh_func;

		std::vector<ParamCoord> m_param_coord;
		std::vector<ChartParamCoord> m_chart_param_coord;

		std::vector<bool> m_face_tex_flag;
		std::vector<int> m_face_tex_group;

		std::vector<double> m_local_distortion;
		std::vector<double> m_vertex_distortion;

		std::vector<int> m_uncorresponding_vtx_array;

		std::vector<int> m_neighbor_charts;

		std::vector< vector<Coord2D> > m_chart_nodes_local_coord; /// each chart's four nodes' local coordinate

	};
}

#endif