#ifndef PARAMETERIZATION_H_
#define PARAMETERIZATION_H_

#include "../Common/BasicDataType.h"
#include <boost/shared_ptr.hpp>

class MeshModel;
class CMeshSparseMatrix;

namespace PARAM
{
	class ParamCoord
	{
	public:
		ParamCoord(){}
		ParamCoord(double s, double t): s_coord(s), t_coord(t){}
		~ParamCoord(){}

		double s_coord;
		double t_coord;
	};

	class ChartParamCoord 
	{
	public:
		ChartParamCoord(double u=0, double v=0, int _chart_id = -1) : param_coord(u, v),
			chart_id(_chart_id){}
		ChartParamCoord(const ParamCoord& _param_coord, int _chart_id = -1) : 
		chart_id(_chart_id), param_coord(_param_coord) {}
		~ChartParamCoord(){}

	public:
		int chart_id;		
		ParamCoord param_coord;
	};

	class SurfaceCoord
	{
	public:
		SurfaceCoord() : face_index(-1), barycentric(0, 0, 0){}
		SurfaceCoord(int fid, const Coord& _coord) : face_index(fid), barycentric(_coord){}
		SurfaceCoord(int fid, double b0, double b1, double b2) : face_index(fid), barycentric(b0, b1, b2) {}
		~SurfaceCoord(){}

		int face_index;
		Coord barycentric;
	};

	void SetCotCoef(const boost::shared_ptr<MeshModel> p_mesh, std::vector<Coord>& cot_coef_vec);
	void SetLapMatrixCoef(const boost::shared_ptr<MeshModel> p_mesh, CMeshSparseMatrix& lap_mat);
	
}

#endif // PARAMTERIZATION_H_