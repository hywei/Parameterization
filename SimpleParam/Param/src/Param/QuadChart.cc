#include "QuadChart.h"

namespace PARAM
{
	QuadChart::QuadChart(int patch_id/* =-1 */) 
		: m_patch_id(patch_id), m_valid_range(std::make_pair(0.0, 1.0)){}

	QuadChart::~QuadChart(){}
}
