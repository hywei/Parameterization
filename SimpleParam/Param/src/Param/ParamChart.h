#ifndef HYWEI_PARAMCHART_H_
#define HYWEI_PARAMCHART_H_

#include "Parameterization.h"
#include <iostream>
#include <cmath>

namespace PARAM
{
	class ParamChart
	{
	public:
        ParamChart(int patch_id=-1)
            :m_patch_id(patch_id), m_valid_range( std::make_pair(0.0, 1.0) ) {}
        
        ~ParamChart(){}

	public:
		bool InValidRangle(const ParamCoord& param_coord) const
		{
            if(m_conner_param_coord_array.size() == 4){
                return GreaterEqual(param_coord.s_coord, m_valid_range.first, LARGE_ZERO_EPSILON)
                    && LessEqual(param_coord.s_coord, m_valid_range.second, LARGE_ZERO_EPSILON)
                    && GreaterEqual(param_coord.t_coord, m_valid_range.first, LARGE_ZERO_EPSILON)
                    && LessEqual(param_coord.t_coord, m_valid_range.second, LARGE_ZERO_EPSILON);
            }else if(m_conner_param_coord_array.size() == 3){
				Coord2D p(param_coord.s_coord, param_coord.t_coord);
				Coord2D a(m_conner_param_coord_array[0].s_coord, m_conner_param_coord_array[0].t_coord);
				Coord2D b(m_conner_param_coord_array[1].s_coord, m_conner_param_coord_array[1].t_coord);
				Coord2D c(m_conner_param_coord_array[2].s_coord, m_conner_param_coord_array[2].t_coord);

				double dist = DistanceToTriangle(p, a, b, c);
				return (dist <= LARGE_ZERO_EPSILON) ;
            }
			return false;			
			
		}

	public:
		int m_patch_id;
		std::vector<ParamCoord> m_conner_param_coord_array;		
		std::pair<double, double> m_valid_range;
	};

    typedef ParamChart QuadChart;
}

#endif
