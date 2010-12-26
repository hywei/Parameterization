#include "QuadDistortion.h"
#include "QuadParameter.h"
#include "Barycentric.h"
#include "../ModelMesh/MeshModel.h"
#include <boost/shared_ptr.hpp>

#include <hj_3rd/hjlib/math/blas_lapack.h>
#include <hj_3rd/zjucad/matrix/lapack.h>
#include <hj_3rd/zjucad/matrix/io.h>

namespace PARAM
{
	QuadDistortion::QuadDistortion(const QuadParameter& quad_parameter) 
		: m_quad_parameter(quad_parameter){}
	QuadDistortion::~QuadDistortion(){}

	void QuadDistortion::ComputeDistortion()
	{
		boost::shared_ptr<MeshModel> p_mesh = m_quad_parameter.GetMeshModel();
		assert(p_mesh != NULL);
		
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		size_t face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();

		m_face_harmonic_distortion.clear(); m_face_harmonic_distortion.resize(face_num);
		m_face_isometric_distortion.clear(); m_face_isometric_distortion.resize(face_num);

		for(size_t fid = 0; fid < face_num; ++fid)
		{
			zjucad::matrix::matrix<double> jacobi_mat = ComputeParamJacobiMatrix(fid);
			m_face_harmonic_distortion[fid] = ComputeHarmonicDistortion(jacobi_mat);
			m_face_isometric_distortion[fid] = ComputeIsometricDistortion(jacobi_mat);
		}
		
	}

	zjucad::matrix::matrix<double> QuadDistortion::ComputeParamJacobiMatrix(int fid) const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_quad_parameter.GetMeshModel();
		const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const IndexArray& faces = face_list_array[fid];

		// Algorithm : Sig2007 parameterization course, p40, equation(4.8)  
		std::vector<Coord> tri_vert_coord_3d(3);
		for(size_t k=0; k<3; ++k) tri_vert_coord_3d[k] = vert_coord_array[faces[k]];

		std::vector<Coord2D> local_coord = ComputeTriangleLocal2DCoord(tri_vert_coord_3d);

		double area_2 = 2*(p_mesh->m_Kernel.GetFaceInfo().GetFaceArea())[fid];

		/// get these three vertices's parameter coordinate
		int chart_id = m_quad_parameter.FindBestChartIDForTriShape(fid);
		std::vector<ParamCoord> vert_param_corod(3);
		for(size_t k=0; k<3; ++k)
		{
			int vid = faces[k];
			int cur_chart_id = m_quad_parameter.GetVertexChartID(vid);
			ParamCoord cur_param_coord = m_quad_parameter.GetVertexParamCoord(vid);
			vert_param_corod[k] = cur_param_coord;
			if(cur_chart_id != chart_id)
			{
				m_quad_parameter.TransParamCoordBetweenCharts(cur_chart_id, chart_id, 
					cur_param_coord, vert_param_corod[k]);
			}
		}

		zjucad::matrix::matrix<double> tm_0(2, 2); 
		tm_0(0, 0) = 0; tm_0(0, 1) = -1; tm_0(1, 0) = 1; tm_0(1, 1) = 0;

		zjucad::matrix::matrix<double> tm_1(2, 3);
		tm_1(0, 0) = local_coord[2][0] - local_coord[1][0]; tm_1(1, 0) = local_coord[2][1] - local_coord[1][1];
		tm_1(0, 1) = local_coord[0][0] - local_coord[2][0]; tm_1(1, 1) = local_coord[0][1] - local_coord[2][1]; 
		tm_1(0, 2) = local_coord[1][0] - local_coord[0][0]; tm_1(1, 2) = local_coord[1][1] - local_coord[0][1];


		zjucad::matrix::matrix<double> tm_2(3, 2);
		for(int k=0; k<3; ++k) { tm_2(k, 0) = vert_param_corod[k].s_coord; tm_2(k, 1) = vert_param_corod[k].t_coord;}
		
		zjucad::matrix::matrix<double> tmp = tm_0 * tm_1;
		tmp = tmp* (1.0 /area_2);
		
		return tmp * tm_2;
	}

	double QuadDistortion::ComputeHarmonicDistortion(const zjucad::matrix::matrix<double>& jacobi_mat) const
	{
		zjucad::matrix::matrix<double> J_tJ = trans(jacobi_mat) * jacobi_mat;

		return 0.5*(J_tJ(0, 0) + J_tJ(1, 1));
	}

	double QuadDistortion::ComputeIsometricDistortion(const zjucad::matrix::matrix<double>& jacobi_mat) const
	{
		return norm( trans(jacobi_mat) * jacobi_mat - zjucad::matrix::eye<double>(2));
	}
}
