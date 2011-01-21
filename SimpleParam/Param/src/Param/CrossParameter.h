#ifndef CROSSPARAMETER_H_
#define CROSSPARAMETER_H_

#include "Parameterization.h"
#include <string>

namespace PARAM
{
	class Parameter;

	class CrossParameter
	{
	public:
		CrossParameter(const Parameter& parameter_1, const Parameter& parameter_2);
		~CrossParameter();

		bool LoadCorrespondingFile(const std::string& corresponding_file);


		void ComputeUnitedDistortionAB();
		void ComputeUnitedDistortionBA();
	    
		void FindCorrespondingAB();
		void FindCorrespondingBA();

		void VertTextureTransferAB();
		void VertTextureTransferBA();
		void FaceTextureTransferAB();
		void FaceTextureTransferBA();

		//! 
		bool GetSurfaceCoordOnA(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const;
		bool GetSurfaceCoordOnB(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const;

		//! get parameter mapping between Surface A and the common domain
		const Parameter& GetParameterA() const { return m_parameter_1; }
		const Parameter& GetParameterB() const { return m_parameter_2; }

		const std::vector<TexCoord>& GetTransferedVertTexArrayA() const { return m_transfer_vert_tex_array_A; }
		const std::vector<TexCoord>& GetTransferedVertTexArrayB() const { return m_transfer_vert_tex_array_B; }
		const std::vector<TexCoordArray>& GetTransferedFaceTexArrayA() const { return m_transfer_face_tex_array_A; }
		const std::vector<TexCoordArray>& GetTransferedFaceTexArrayB() const { return m_transfer_face_tex_array_B; }

		const std::vector<int>& GetUnCorrespondingVertArrayOnA() const { return m_uncorresponding_vert_array_A; }
		const std::vector<int>& GetUnCorrespondingVertArrayOnB() const { return m_uncorresponding_vert_array_B; }

		std::vector<SurfaceCoord> m_corresponding_AB;

	private:
		ChartParamCoord GetChartParamCoord4CorrespondingChartOnA(const ChartParamCoord& chart_param_coord_onB) const;
		ChartParamCoord GetChartParamCoord4CorrespondingChartOnB(const ChartParamCoord& chart_param_coord_onA) const;

	private:
		const Parameter& m_parameter_1;
		const Parameter& m_parameter_2;

		std::vector<SurfaceCoord> m_corresponding_BA;

		std::vector<TexCoord> m_transfer_vert_tex_array_A;
		std::vector<TexCoord> m_transfer_vert_tex_array_B;
		std::vector<TexCoordArray> m_transfer_face_tex_array_A;
		std::vector<TexCoordArray> m_transfer_face_tex_array_B;

		//! for debug
		std::vector<int> m_uncorresponding_vert_array_A;
		std::vector<int> m_uncorresponding_vert_array_B;

		std::map<int, int> m_conner_correspondingAB;
		std::map<int, int> m_conner_correspondingBA;
		std::map<int, int> m_edge_correspondingAB;
		std::map<int, int> m_edge_correspondingBA;
		std::map<int, int> m_patch_correspondingAB;
		std::map<int, int> m_patch_correspondingBA;
	};
}

#endif //CROSSPARAMETER_H_