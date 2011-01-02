#ifndef TRANSFUNCTOR_H_
#define TRANSFUNCTOR_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <hj_3rd/zjucad/matrix/matrix.h>

namespace PARAM
{
    class ChartCreator;
	class ParamCoord;

    class TransFunctor
    {
    public:
        TransFunctor(boost::shared_ptr<ChartCreator> _p_chart_creator);
        ~TransFunctor();

        //! get the transition matrix between two charts
		zjucad::matrix::matrix<double> GetTransMatrix(int from_chart_id, int to_chart_id) const;

		//! get the transition matrix between two ambiguity charts, the ambiguity charts is
		//! that two charts share two or more edges
		zjucad::matrix::matrix<double> GetTransMatrixBetweenAmbiguityCharts(int from_vert, 
			int from_chart_id, int to_vert, int to_chart_id) const;

		void TransParamCoordBetweenAmbiguityCharts(int from_vert, int from_chart_id, int to_chart_id, 
			const ParamCoord& from_coord, ParamCoord& to_coord) const;
        
	private:
		//! get the transition list between two charts
		bool GetTranslistBetweenTwoCharts(int from_chart_id, int to_chart_id, std::vector<int>& trans_list) const;

		//! get the edge transition lists between two vertices, if the two vertices in
		//! same chart, then the lists is empty, else find the shortest path between them 
		bool GetPatchEdgeTransLists(int from_vid, int to_vid, std::vector<int>& edge_trans_list) const;

		//! get the matrix between two adjacent charts
		zjucad::matrix::matrix<double> GetTransMatrixOfAdjCharts(int from_chart, int to_chart) const;
		//! get the transition matrix in one chart from a frame to another frame
		zjucad::matrix::matrix<double> GetTransMatrixInOneChart(int chart_id, 
			std::pair<int,int> old_x_axis, std::pair<int,int> new_x_axis) const;

		//! choose which edge to transition for ambiguity charts 
		int ChooseEdgeForAmbiguityChart( const std::vector<int>& common_edge_index_array,
			 int from_vert, int to_vert=-1) const;

	private:
		boost::shared_ptr<ChartCreator> p_chart_creator;

    };
}

#endif
