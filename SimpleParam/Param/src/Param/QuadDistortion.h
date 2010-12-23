#ifndef QUADDISTORTION_H_
#define QUADDISTORTION_H_

#include "..\hj_3rd\include\zjucad\matrix\matrix.h"
#include <vector>

namespace PARAM
{
	class QuadParameter;

	class QuadDistortion
	{
	public:
		QuadDistortion(const QuadParameter& quad_parameter);
		~QuadDistortion();

		void ComputeDistortion();

	public:
		//! IO 
		const std::vector<double>& GetFaceHarmonicDistortion() const { return m_face_harmonic_distortion; }
		const std::vector<double>& GetFaceIsometricDistortion() const{ return m_face_isometric_distortion; }

	private:
		double ComputeHarmonicDistortion(const zjucad::matrix::matrix<double>& jacobi_mat) const;
		double ComputeIsometricDistortion(const zjucad::matrix::matrix<double>& jacobi_mat) const;

		//! compute the jacobi matrix from surface to parameter domain
		zjucad::matrix::matrix<double> ComputeParamJacobiMatrix(int fid) const;

	private:
		const QuadParameter& m_quad_parameter;
		std::vector<double> m_face_harmonic_distortion; //! each face's harmonic map distortion
		std::vector<double> m_face_isometric_distortion; //! each face's isometric map distortion

	};
}

#endif //QUADDISTORTION