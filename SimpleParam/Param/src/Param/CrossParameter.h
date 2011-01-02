#ifndef CROSSPARAMETER_H_
#define CROSSPARAMETER_H_

#include "Parameterization.h"

namespace PARAM
{
	class QuadParameter;

	class CrossParameter
	{
	public:
		CrossParameter(const QuadParameter& quad_parameter_1, const QuadParameter& quad_parameter_2);
		~CrossParameter();

		void ComputeUnitedDistortionAB();
		void ComputeUnitedDistortionBA();
	    
		void FindCorrespondingAB();
		void FindCorrespondingBA();

		//! 
		bool GetSurfaceCoordOnA(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const;
		bool GetSurfaceCoordOnB(const ChartParamCoord& chart_param_coord, SurfaceCoord& surface_coord) const;

		//! get parameter mapping between Surface A and the common domain
		const QuadParameter& GetQuadParameterA() const { return m_quad_parameter_1; }
		const QuadParameter& GetQuadParameterB() const { return m_quad_parameter_2; }

		const std::vector<int>& GetUnCorrespondingVertArrayOnA() const { return m_uncorresponding_vert_array_A; }
		const std::vector<int>& GetUnCorrespondingVertArrayOnB() const { return m_uncorresponding_vert_array_B; }

	private:
		const QuadParameter& m_quad_parameter_1;
		const QuadParameter& m_quad_parameter_2;

		//! for debug
		std::vector<int> m_uncorresponding_vert_array_A;
		std::vector<int> m_uncorresponding_vert_array_B;

	};
}

#endif //CROSSPARAMETER_H_