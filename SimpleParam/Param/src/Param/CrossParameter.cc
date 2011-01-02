#include "CrossParameter.h"
#include "QuadParameter.h"
#include "../ModelMesh/MeshModel.h"

namespace PARAM
{

	CrossParameter::CrossParameter(const QuadParameter& quad_parameter_1, const QuadParameter& quad_parameter_2)
		: m_quad_parameter_1(quad_parameter_1), m_quad_parameter_2(quad_parameter_2){}
	CrossParameter::~CrossParameter(){}

	bool CrossParameter::GetSurfaceCoordOnA(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const
	{
		if(!m_quad_parameter_1.FindCorrespondingOnSurface(chart_param_coord, surface_coord))
		{
			printf("@@@Error: Can't find corresponding on surface A!\n");
			return false;
		}
		return true;
	}

	bool CrossParameter::GetSurfaceCoordOnB(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const
	{				
		if(!m_quad_parameter_2.FindCorrespondingOnSurface(chart_param_coord, surface_coord))
		{
			printf("@@@Error: Can't find corresponding on surface B!\n");
			return false;
		}
		return true;
	}


	void CrossParameter::FindCorrespondingAB()
	{
		const boost::shared_ptr<MeshModel> p_mesh_1 = m_quad_parameter_1.GetMeshModel();
		int vert_num = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		m_uncorresponding_vert_array_A.clear();
		for(int vid=0; vid < vert_num; ++vid)
		{
			int chart_id = m_quad_parameter_1.GetVertexChartID(vid);
			PARAM::ParamCoord param_coord = m_quad_parameter_1.GetVertexParamCoord(vid);

			PARAM::SurfaceCoord surface_coord;
			if(!GetSurfaceCoordOnB(PARAM::ChartParamCoord(param_coord, chart_id), surface_coord))
			{
				m_uncorresponding_vert_array_A.push_back(vid);
			}
		}
	}

	void CrossParameter::FindCorrespondingBA()
	{
		const boost::shared_ptr<MeshModel> p_mesh_2 = m_quad_parameter_2.GetMeshModel();
		int vert_num = p_mesh_2->m_Kernel.GetModelInfo().GetVertexNum();

		m_uncorresponding_vert_array_B.clear();
		for(int vid=0; vid < vert_num; ++vid)
		{
			int chart_id = m_quad_parameter_2.GetVertexChartID(vid);
			PARAM::ParamCoord param_coord = m_quad_parameter_2.GetVertexParamCoord(vid);

			PARAM::SurfaceCoord surface_coord;
			if(!GetSurfaceCoordOnB(PARAM::ChartParamCoord(param_coord, chart_id), surface_coord))
			{
				m_uncorresponding_vert_array_B.push_back(vid);
			}
		}
	}
}