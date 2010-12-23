#include "CrossParam.h"
#include "ParamDistortion.h"

#include "../ModelMesh/MeshModel.h"
#include <boost/shared_ptr.hpp>
#include <numeric>

#include "../hj_3rd/include/math/blas_lapack.h"
#include "../hj_3rd/include/zjucad/matrix/lapack.h"

namespace PARAM
{

	CrossParam::CrossParam(QuadParam& param1, QuadParam& param2) 
		: m_param_1(param1), m_param_2(param2)
	{}

	CrossParam::~CrossParam(){}

	void CrossParam::ComputeUintedDistortion()
	{
		printf("Compute United Distortion!\n");

		FindVtxMappingForVtxInSurface1();
		FindVtxMappingForVtxInSurface2();

		ParamDistortion param_distortion_1(m_param_1);
		ParamDistortion param_distortion_2(m_param_2);

		const boost::shared_ptr<MeshModel> p_mesh_1 = m_param_1.GetMeshModel();
		const PolyIndexArray& face_index_array_1 = p_mesh_1->m_Kernel.GetFaceInfo().GetIndex();
		const CoordArray& vtx_coord_array_1 = p_mesh_1->m_Kernel.GetVertexInfo().GetCoord();

		const boost::shared_ptr<MeshModel> p_mesh_2 = m_param_2.GetMeshModel();
	    const PolyIndexArray& face_index_array_2 = p_mesh_2->m_Kernel.GetFaceInfo().GetIndex();
		int face_num_2 = p_mesh_2->m_Kernel.GetModelInfo().GetFaceNum();
		
		m_united_distortion.clear();
		m_united_distortion.resize(face_num_2);

		for(int k=0; k<face_num_2; ++k)
		{
			const IndexArray& face_index = face_index_array_2[k];
			/// the three vertices's position coordinate in mesh 1
			vector<Coord> new_vtx_pos_coord(3);
			/// the three vertices's parameter  coordinate, same in two mesh 
			vector<ParamCoord> vtx_param_coord(3);

			/// for computing the jacobi matrix, we need three vertices's parameter coordinate in the same chart
			//int chart_id = m_param_2.GetFaceChartID(k);
			int chart_id = m_param_2.GetVertexChartID(face_index[0]);
		
			for(int i=0; i<3; ++i)
			{
				int vtx_idx = face_index[i];
				ParamCoord param_coord = m_param_2.GetParamCoord(vtx_idx);
				int cur_chart_id = m_param_2.GetVertexChartID(vtx_idx);
				if(cur_chart_id != chart_id)
				{
					m_param_2.TransParamCoordBetweenTwoChart(cur_chart_id, chart_id, param_coord, vtx_param_coord[i]);
				}else
				{
					vtx_param_coord[i] = param_coord;
				}

				int face_id_in_mesh_1 = m_vtx_face_in_mesh_1[vtx_idx];
				Coord barycentric = m_vtx_barycentric_in_mesh_1[vtx_idx];

				const IndexArray& face_in_mesh_1 = face_index_array_1[face_id_in_mesh_1];
				vector<Coord> vtx_coord_in_mesh_1(3);
				for(size_t j=0; j<3; ++j)
				{
					int vtx = face_in_mesh_1[j];
					vtx_coord_in_mesh_1[j] = vtx_coord_array_1[vtx];
				}
				new_vtx_pos_coord[i] = ParamDistortion::ComputePosWithBarycentric(vtx_coord_in_mesh_1, barycentric);				
			}
		
			/// J_m = (J_m_2)^-1 * J_m_1
			const zjucad::matrix::matrix<double>& jacobi_1 = param_distortion_1.ComputeTriJacobiMatrix(new_vtx_pos_coord, vtx_param_coord);
			zjucad::matrix::matrix<double> inv_jacobi_1(jacobi_1);
			if(!inv(inv_jacobi_1))
			{
				printf("compute inv jacobi matrix 1 fail.%d\n", k);
			}
			const zjucad::matrix::matrix<double>& jacobi_2 = param_distortion_2.ComputeTriJacobiMatrix(k);

			//printf("sum: %d: %lf, %lf\n", k, sum(jacobi_2), sum(inv_jacobi_1));
			zjucad::matrix::matrix<double> united_jacobi = jacobi_2*inv_jacobi_1;
			zjucad::matrix::matrix<double> J_tJ = trans(united_jacobi)*united_jacobi;	
			//inv(J_tJ);
			
			zjucad::matrix::matrix<double> i_mat (zjucad::matrix::eye<double>(2));
			
			/// isometric maps
			m_united_distortion[k] = norm(J_tJ - i_mat); 

			/// harmonic maps
			m_united_distortion[k] = 0.5*(J_tJ(0, 0) + J_tJ(1, 1));

		}
		   	
		vector<double>::iterator min_it = min_element(m_united_distortion.begin(), m_united_distortion.end());
		vector<double>::iterator max_it = max_element(m_united_distortion.begin(), m_united_distortion.end());

		printf("The min distortion is %d %lf\n", (min_it - m_united_distortion.begin()), *min_it);
		printf("The max distortion is %d %lf\n", (max_it - m_united_distortion.begin()), *max_it);

		printf("The total distortion is %lf\n", accumulate(m_united_distortion.begin(), m_united_distortion.end(), 0.0));

		//copy(m_united_distortion.begin(), m_united_distortion.end(), ostream_iterator<double>(cout, " "));
	}

	void CrossParam::FindVtxMappingForVtxInSurface1()
	{
		ParamDistortion param_distortion(m_param_2);

		const boost::shared_ptr<MeshModel> p_mesh_1 = m_param_1.GetMeshModel();
		int vert_num_1 = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		m_vtx_face_in_mesh_2.clear();
		m_vtx_face_in_mesh_2.resize(vert_num_1, -1);

		m_vtx_barycentric_in_mesh_2.clear();
		m_vtx_barycentric_in_mesh_2.resize(vert_num_1);

		m_uncorresponding_vtx_array.clear();

		for(int vid=0; vid < vert_num_1; ++vid)
		{
			ChartParamCoord chart_param_coord (m_param_1.GetParamCoord(vid), m_param_1.GetVertexChartID(vid));
			SurfaceCoord surface_coord;
			bool ret = param_distortion.ComputeSurfaceCoord(chart_param_coord, surface_coord);

			if(ret == false)
			{
				printf("@@@@@@@@Error, can't find the correspoinding vertex!@@@@@@@@\n");
				m_uncorresponding_vtx_array.push_back(vid);
				printf("barycentric : %lf %lf %lf\n", surface_coord.barycentric[0], 
					surface_coord.barycentric[1], surface_coord.barycentric[2]);
			}

			m_vtx_face_in_mesh_2[vid] = surface_coord.face_index;
			m_vtx_barycentric_in_mesh_2[vid] = surface_coord.barycentric;
		}
	}

	void CrossParam::FindVtxMappingForVtxInSurface2()
	{
		ParamDistortion param_distortion(m_param_1);

		const boost::shared_ptr<MeshModel> p_mesh_2 = m_param_2.GetMeshModel();
		int vert_num_2 = p_mesh_2->m_Kernel.GetModelInfo().GetVertexNum();

		m_vtx_face_in_mesh_1.clear(); 
		m_vtx_face_in_mesh_1.resize(vert_num_2, -1);
		m_vtx_barycentric_in_mesh_1.clear(); 
		m_vtx_barycentric_in_mesh_1.resize(vert_num_2);

		m_uncorresponding_vtx_array.clear();
		
		for(int vid = 0; vid < vert_num_2; ++vid)
		{		   
			ChartParamCoord chart_param_coord (m_param_2.GetParamCoord(vid), m_param_2.GetVertexChartID(vid));
			SurfaceCoord surface_coord;
			bool ret = param_distortion.ComputeSurfaceCoord(chart_param_coord, surface_coord);
			
			if(ret == false)
			{
				printf("@@@@@@@@Error, can't find the correspoinding vertex!@@@@@@@@\n");
				m_uncorresponding_vtx_array.push_back(vid);
				printf("barycentric : %lf %lf %lf\n", surface_coord.barycentric[0], 
					surface_coord.barycentric[1], surface_coord.barycentric[2]);
			}

			m_vtx_face_in_mesh_1[vid] = surface_coord.face_index;
			m_vtx_barycentric_in_mesh_1[vid] = surface_coord.barycentric;
		}
		
	}
	
	void CrossParam::ComputeTexCoordOnSurface1()
	{
		const boost::shared_ptr<MeshModel> p_mesh_2 = m_param_2.GetMeshModel();
		const PolyIndexArray& face_index_array_2 = p_mesh_2->m_Kernel.GetFaceInfo().GetIndex();
		const TexCoordArray& vtx_tex_coord_array_2 = p_mesh_2->m_Kernel.GetVertexInfo().GetTexCoord();

		const boost::shared_ptr<MeshModel> p_mesh_1 = m_param_1.GetMeshModel();
		int vert_num_1 = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();

		m_tex_coord_1.clear();
		m_tex_coord_1.resize(vert_num_1);

		for(int vid=0; vid<vert_num_1; ++vid)
		{
			ParamCoord param_coord = m_param_1.GetParamCoord(vid);
			int cur_chart_id = m_param_1.GetVertexChartID(vid);

			int fid_in_mesh_2 = m_vtx_face_in_mesh_2[vid];
			Coord barycentric = m_vtx_barycentric_in_mesh_2[vid];

			const IndexArray& face_index = face_index_array_2[fid_in_mesh_2];
			vector<TexCoord> tex_coord(3);

			for(size_t k=0; k<3; ++k)
			{
				int v = face_index[k];
				tex_coord[k] = vtx_tex_coord_array_2[v];
			}
			m_tex_coord_1[vid] = tex_coord[0]*barycentric[0] + tex_coord[1]*barycentric[1] + tex_coord[2]*barycentric[2];
		}

	}

	void CrossParam::ComputeTexCoordOnSurface2()
	{
		const boost::shared_ptr<MeshModel> p_mesh_1 = m_param_1.GetMeshModel();
		const PolyIndexArray& face_index_array_1 = p_mesh_1->m_Kernel.GetFaceInfo().GetIndex();
		const TexCoordArray& vtx_tex_coord_array_1 = p_mesh_1->m_Kernel.GetVertexInfo().GetTexCoord();

		const boost::shared_ptr<MeshModel> p_mesh_2 = m_param_2.GetMeshModel();
		
		int vert_num_2 = p_mesh_2->m_Kernel.GetModelInfo().GetVertexNum();
		m_tex_coord_2.clear();
		m_tex_coord_2.resize(vert_num_2);
		for(int vid=0; vid<vert_num_2; ++vid)
		{
			ParamCoord param_coord = m_param_2.GetParamCoord(vid);
			int cur_chart_id = m_param_2.GetVertexChartID(vid);

			int fid_in_mesh_1 = m_vtx_face_in_mesh_1[vid];
			Coord barycentric = m_vtx_barycentric_in_mesh_1[vid];
						
			const IndexArray& face_index = face_index_array_1[fid_in_mesh_1];
			vector<TexCoord> tex_coord(3);
			for(size_t k=0; k<3; ++k)
			{
				int v = face_index[k];
				tex_coord[k] = vtx_tex_coord_array_1[v];
			}

			m_tex_coord_2[vid] = tex_coord[0]*barycentric[0] + tex_coord[1]*barycentric[1] + tex_coord[2]*barycentric[2];
		}
	}
}