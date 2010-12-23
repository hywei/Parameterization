#ifndef QUADTRANSFUNCTOR_H_
#define QUADTRANSFUNCTOR_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include "../hj_3rd/include/zjucad/matrix/matrix.h"

namespace PARAM
{
	class QuadChartCreator;

	class QuadTransFunctor
	{
	public:
		QuadTransFunctor(boost::shared_ptr<QuadChartCreator> quad_chart_creator);
		~QuadTransFunctor();

		//! get the transition matrix between two charts
		zjucad::matrix::matrix<double> GetTransMatrix(int from_chart_id, int to_chart_id) const;
	private:
		//! get the transition list between two charts
		bool GetTranslistBetweenTwoCharts(int from_chart_id, int to_chart_id, std::vector<int>& trans_list) const;

		//! get the matrix between two adjacent charts
		zjucad::matrix::matrix<double> GetTransMatrixOfAdjCharts(int from_chart, int to_chart) const;
		//! get the transition matrix in one chart from a frame to another frame
		zjucad::matrix::matrix<double> GetTransMatrixInOneChart(int chart_id, 
			std::pair<int,int> old_x_axis, std::pair<int,int> new_x_axis) const;

	private:
		boost::shared_ptr<QuadChartCreator> p_quad_chart_creator;

	};
}

#endif //QUADTRANSFUNCTOR_H_