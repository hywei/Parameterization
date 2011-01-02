#ifndef CROSSPARAMETER_H_
#define CROSSPARAMETER_H_

//#include "../Param/CrossParam.h"
#include "../Param/QuadParameter.h"

class CrossParamWrapper
{
public:
	CrossParamWrapper();
	~CrossParamWrapper();
	
	void ComputeUnitedDistortion();

public:
	//PARAM::QuadParameter m_param_1;
	//PARAM::QuadParameter m_param_2;

private:
	//PARAM::CrossParam m_cross_param;

};

#endif