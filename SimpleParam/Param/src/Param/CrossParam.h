#ifndef CROSSPARAM_H_
#define CROSSPARAM_H_

#include "QuadParam.h"
#include <hj_3rd/zjucad/matrix/matrix.h>

namespace PARAM
{
	class CrossParam
	{
	public:
		CrossParam(QuadParam& param1, QuadParam& param2);
		~CrossParam();


		void ComputeUintedDistortion();
		

		void FindVtxMappingForVtxInSurface2();		 
		void FindVtxMappingForVtxInSurface1();

		void ComputeTexCoordOnSurface2();
		void ComputeTexCoordOnSurface1();
	

	public:
		QuadParam& GetQuadParam1()  { return m_param_1; }
		QuadParam& GetQuadParam2()  { return m_param_2; }
	//	const QuadParam& GetQuadParam1() const { return m_param_1; }
	//	const QuadParam& GetQuadParam2() const { return m_param_2; }

		const std::vector<double>& GetUnitedDistortion() const 
		{
			return m_united_distortion;		
		}

		const std::vector<int>& GetUnCorrespondingVtxArray() const
		{
			return m_uncorresponding_vtx_array;
		}

		const std::vector<TexCoord>& GetTexCoordOnSurface1()
		{
			return m_tex_coord_1;
		}

		const std::vector<TexCoord>& GetTexCoordOnSurface2()
		{
			return m_tex_coord_2;
		}
		
	private:
		QuadParam& m_param_1;
		QuadParam& m_param_2;

		std::vector<double> m_united_distortion;

		std::vector<int> m_uncorresponding_vtx_array;

		std::vector<TexCoord> m_tex_coord_2;
		std::vector<int> m_vtx_face_in_mesh_1;
		std::vector<Coord> m_vtx_barycentric_in_mesh_1;

		std::vector<TexCoord> m_tex_coord_1;
		std::vector<int> m_vtx_face_in_mesh_2;
		std::vector<Coord> m_vtx_barycentric_in_mesh_2;
	};

}

#endif
