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
                double x = param_coord.s_coord, y = param_coord.t_coord;
                //std::cout<< x <<" " << y << std::endl;
                bool flag1 = (GreaterEqual(x, 0 , LARGE_ZERO_EPSILON) && LessEqual( x, 1, LARGE_ZERO_EPSILON));
                bool flag2 = LessEqual( y - sqrt(3.0)*x, 0, LARGE_ZERO_EPSILON);
                bool flag3 = LessEqual( y + sqrt(3.0)*(x - 1), 0, LARGE_ZERO_EPSILON);
                bool flag4 = GreaterEqual(y, 0, LARGE_ZERO_EPSILON);
                return flag1 && flag2 && flag3 && flag4;
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
