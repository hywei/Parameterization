#ifndef QUADPARAMETER_H_
#define QUADPARAMETER_H_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "Parameterization.h"

class MeshModel;
class LinearSolver;
class CMeshSparseMatrix;

namespace PARAM
{
	class QuadChartCreator;

	class QuadParameter
	{
	public:
		QuadParameter(boost::shared_ptr<MeshModel> _p_mesh);
		~QuadParameter();
		
		bool ComputeParamCoord();

	public:
		//! IO Function
		bool LoadQuadFile(const std::string& file_name);

		int GetVertexChartID(int vid) const { return m_vert_chart_array[vid]; }
		ParamCoord GetVertexParamCoord(int vid)  const { return m_vert_param_coord_array[vid]; }

		const std::vector<int>& GetVertexChartArray() const { return m_vert_chart_array; }
		const std::vector<ParamCoord>& GetVertexParamCoordArray() const { return m_vert_param_coord_array; }

		boost::shared_ptr<MeshModel> GetMeshModel() const { return p_mesh; }
		boost::shared_ptr<QuadChartCreator> GetChartCreator() const { return p_quad_chart_creator; }

		const std::vector<int>& GetOutRangeVertArray() const { return m_out_range_vert_array; }
		const std::vector<int>& GetUnSetFaceArray() const { return m_unset_layout_face_array; }		

	private:
		void SetVariIndexMapping(std::vector<int>& vari_index_mapping);
		void SetBoundaryVertexParamValue();

		void SetConstraint(LinearSolver& linear_solver);
		void SolveQuadParameter(const CMeshSparseMatrix& lap_mat);
		void VertexRelatextion();

		//! after each iterator, we need reassign vertices's chart  
		void AdjustPatchBoundary();

		void SetInitVertChartLayout();
		void SetInitFaceChartLayout();
		
	  
		void GetOutRangeVertices(std::vector<int>& out_range_vert_array) const;
		bool FindValidChartForOutRangeVertex(int our_range_vert, int max_ringe_num = 5);
		double ComputeOutRangeError(ParamCoord param_coord)
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

		//! get the length of a mesh path
		double ComputeMeshPathLength(const std::vector<int>& mesh_path, int start_idx, int end_idx) const;

		//! get a conner's index in a patch/chart
		int GetConnerIndexInPatch(int conner_id, int patch_id);

		//! reset face chart layout
		void ResetFaceChartLayout();

		//! set mesh texture 
		void SetMeshFaceTextureCoord();

		//! set each chart's vertices 
		void SetChartVerticesArray();

		//! find the corresponding surface position on chart with the chart parameter coordinate 
		bool FindCorrespondingInChart(const ChartParamCoord& chart_param_coord, 
			int chart_id, SurfaceCoord& surface_coord) const;

	public:
		void TransParamCoordBetweenCharts(int from_chart_id, int to_chart_id, 
			const ParamCoord& from_param_coord, ParamCoord& to_param_coord) const;

		//! for one face's three vertices, they may be in different chart, and we need transite them to same chart.
		//! so there are several choose for the common chart, and our standart is the best triangle shape in parameter domain. 
		int FindBestChartIDForTriShape(int fid) const;

		//! compute harmonic maps and isometric maps distortion
		void ComputeDistortion();

		//! find the corresponding surface position with the chart paramerter coordinate
		bool FindCorrespondingOnSurface(const ChartParamCoord& chart_param_coord,
			SurfaceCoord& surface_coord) const;
	
		void FaceValue2VtxColor(const std::vector<double>& face_value);

	private:
		boost::shared_ptr<MeshModel> p_mesh;
		boost::shared_ptr<QuadChartCreator> p_quad_chart_creator;

		std::vector<int> m_vert_chart_array;
		std::vector<int> m_face_chart_array;

		std::vector<ParamCoord> m_vert_param_coord_array;

		std::vector< std::vector<int> > m_chart_vertices_array; //! 

		//! for debug
		std::vector<int> m_out_range_vert_array;
		std::vector<int> m_unset_layout_face_array;
		
	};
}

#endif
