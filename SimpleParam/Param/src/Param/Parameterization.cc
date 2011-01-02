#include "Parameterization.h"
#include "../ModelMesh/MeshModel.h"
#include "../Numerical/MeshSparseMatrix.h"
#include "../Common/HSVColor.h"
#include <limits>
#include <set>
#include <queue>

namespace PARAM
{
	void SetCotCoef(const boost::shared_ptr<MeshModel> p_mesh, std::vector<Coord>& cot_coef_vec)
	{
		if(p_mesh == NULL) return;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		size_t nFace = fIndex.size();
		cot_coef_vec.clear();
		cot_coef_vec.resize(nFace);

		size_t i, j;
		Coord e[3];
		double a[3];
		for(i = 0; i < nFace; ++ i)
		{
			const IndexArray& f = fIndex[i];
			for(j = 0; j < 3; ++ j)
			{
				e[j] = (vCoord[f[(j+1)%3]] - vCoord[f[j]]).unit();
			}

			a[0] = angle(e[0], -e[2]);
			a[1] = angle(e[1], -e[0]);
			a[2] = PI-a[0]-a[1];

			for(j = 0; j < 3; ++ j)
			{
				// angle < 1 or angle > 179
				if(fabs(a[j]-0.0) < 0.0174 || fabs(a[j]-PI) < 0.0174)
				{
					cot_coef_vec[i][j] = 57.289;   // atan(1)
				}
				else
				{
					cot_coef_vec[i][j] = 1.0/tan(a[j]);
				}
			}
		}

		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();
		// Check the alpha+belta < PI is satisfied or not
		const PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		int nAdjust = 0;
		int k, h;
		for(i = 0; i < (size_t)vert_num ; ++i)
		{
			const IndexArray& adjf = vAdjFaces[i];
			size_t n = adjf.size();
			size_t begin_idx = 0, end_idx = n-1;
			if(p_mesh->m_BasicOp.IsBoundaryVertex((int) i))
			{
				begin_idx = 1;
				end_idx = n-1;
			}
			for(j = begin_idx; j <= end_idx; ++ j)
			{
				FaceID fID1 = adjf[(j+n-1)%n];
				FaceID fID2 = adjf[j];
				const IndexArray& f1 = fIndex[fID1];
				const IndexArray& f2 = fIndex[fID2];
				for( k = 0; k < 3; ++ k)
				{
					if(f1[k] == i)
						break;
				}
				for( h = 0; h < 3; ++ h)
				{
					if(f2[h] == i)
						break;
				}
				if(cot_coef_vec[fID1][(k+1)%3] + cot_coef_vec[fID2][(h+2)%3] < 0.0)
				{
					cot_coef_vec[fID1][(k+1)%3] = cot_coef_vec[fID2][(h+2)%3] = 0.0;
					++ nAdjust;
				}
			}
		}
		printf("#Cotangent Weight Adjust = %d\n", nAdjust);
	}

	void SetLapMatrixCoef(const boost::shared_ptr<MeshModel> p_mesh, CMeshSparseMatrix& lap_matrix)
	{
		std::vector<Coord> cot_coef_vec;
		SetCotCoef(p_mesh, cot_coef_vec);

		const PolyIndexArray& vAdjFaces = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		int i, j, k, n;
		int fID, vID;
		int row, col;
		double coef, d;

		size_t vert_num = p_mesh->m_Kernel.GetVertexInfo().GetCoord().size();
		//
		int nb_variables_ = (int) vert_num;
		lap_matrix.SetRowCol(nb_variables_, nb_variables_);

		for(i = 0; i < nb_variables_; ++ i)
		{
			const IndexArray& adjFaces = vAdjFaces[i];
			row = i;
			n = (int) adjFaces.size();

			for (j = 0; j < n; j++)
			{
				fID = adjFaces[j];
				const IndexArray& f = fIndex[fID];

				// Find the position of vertex i in face fID
				for(k = 0; k < 3; ++ k)
				{
					if(f[k] == i) break;
				}
				assert(k!=3);

				vID = f[(k+1)%3];
				col = vID;

				coef = cot_coef_vec[fID][(k+2)%3];
				lap_matrix.GetElement(row, col, d);
				d -= coef;
				lap_matrix.SetElement(row, col, d);

				vID = f[(k+2)%3];
				col = vID;

				coef = cot_coef_vec[fID][(k+1)%3];
				lap_matrix.GetElement(row, col, d);
				d -= coef;
				lap_matrix.SetElement(row, col, d);

				lap_matrix.GetElement(row, row, d);
				d += cot_coef_vec[fID][(k+1)%3] + cot_coef_vec[fID][(k+2)%3];
				lap_matrix.SetElement(row, row, d);
			}
		}
	}
	

	void SetLapMatrixCoefWithMeanValueCoord(const boost::shared_ptr<MeshModel> p_mesh, 
		CMeshSparseMatrix& lap_matrix)
	{
		if(p_mesh == NULL) return;
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		std::vector< std::vector<double> > vert_tan_coef(vert_num);
		const PolyIndexArray& vert_adj_vertices = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		for(int vid=0; vid < vert_num; ++vid){
			if(p_mesh->m_BasicOp.IsBoundaryVertex(vid)) continue;
			const IndexArray& adj_vert_array = vert_adj_vertices[vid];
			for(size_t i=0; i<adj_vert_array.size(); ++i)
			{
				int cur_vtx = adj_vert_array[i];
				int nxt_vtx = adj_vert_array[(i+1)%adj_vert_array.size()];
				Coord e1 = vCoord[cur_vtx] - vCoord[vid];
				Coord e2 = vCoord[nxt_vtx] - vCoord[vid];

				double a = angle(e1, e2)/2;
				vert_tan_coef[vid].push_back(tan(a));
			}
		}	   
	 
		//
		int nb_variables_ = vert_num;
		lap_matrix.SetRowCol(nb_variables_, nb_variables_);

		for(int i = 0; i < nb_variables_; ++ i)
		{
			if(p_mesh->m_BasicOp.IsBoundaryVertex(i)) continue;
			const IndexArray& adjVertices = vert_adj_vertices[i];
			int row = i;

			int adj_num = adjVertices.size();
			for(int j=0; j<adj_num; ++j)
			{
				int col = adjVertices[j];
				double edge_len = (vCoord[row]-vCoord[col]).abs();
				double coef = (vert_tan_coef[i][j] + vert_tan_coef[i][(j+adj_num-1)%adj_num])/edge_len;
				double cur_value;
				lap_matrix.GetElement(row, col, cur_value);
				cur_value -= coef;
				lap_matrix.SetElement(row, col, cur_value);

				lap_matrix.GetElement(row, row, cur_value);
				cur_value += coef;
				lap_matrix.SetElement(row, row, cur_value);

			}
		
		}

	}
		
	double xmult(double x1,double y1,double x2,double y2,double x0,double y0){
			return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
	}
	
	double area_triangle(double x1,double y1,double x2,double y2,double x3,double y3){
		return fabs(xmult(x1,y1,x2,y2,x3,y3))/2;
	}
	double dis_ptoline(double x1,double y1,double x2,double y2,double ex,double ey)
	{ 
		double dis,tem1,tem2,t1,t2,yd=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		t2=sqrt((x2-ex)*(x2-ex)+(y2-ey)*(y2-ey));
		t1=sqrt((x1-ex)*(x1-ex)+(y1-ey)*(y1-ey));
		dis=area_triangle(x1,y1,x2,y2,ex,ey)*2/yd;
		tem1=sqrt(t1*t1-dis*dis);
		tem2=sqrt(t2*t2-dis*dis);

		if (tem1>yd||tem2>yd) {
			if (t1>t2) dis = t2;
			else dis = t1;
		}
		return dis;
	}

	void FaceValue2VtxColor(boost::shared_ptr<MeshModel> p_mesh, const std::vector<double>& face_value)
	{
		ColorArray& face_color_array = p_mesh->m_Kernel.GetFaceInfo().GetColor();
		face_color_array.clear();

		CHSVColor color;
		color.RGBtoHSV(1, 0, 0);
		color.m_S = 0.9f;
		color.m_V = 0.9f;

		size_t face_num = p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();

		double min = 1e20, max = -1e20;
		for(size_t i = 0; i < face_num; ++i)
		{
			if(face_value[i] < min) min = face_value[i];
			if(face_value[i] > max) max = face_value[i];
		}

		double range = (max - min) * 1.1;

		bool eq = ALMOST_EQUAL_LARGE(range, 0.0);

		for(size_t i = 0; i < face_value.size();++i)
		{
			float R, G, B;
			if (eq)
			{
				color.m_H = (float) 0.5 * 255;
			}
			else
			{
				double prop = (face_value[i] - min) / range;

				//color.m_S = prop;
				color.m_H = (float)(1 - prop) * 255;
			}
			color.HSVtoRGB(&R, &G, &B);
			Color c(R, G, B);
			face_color_array.push_back(c);
		}
	}

	double GetNearestVertexOnPath(boost::shared_ptr<MeshModel> p_mesh, int from_vert, const std::vector<int>& path, int& nearest_vid)
	{
		if(p_mesh == NULL) return -1;

// 		const PolyIndexArray& adj_vertices = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
// 
// 		std::set<int> path_vert_set(path.begin(), path.end());
// 		std::set<int> visited_vert_set;
// 		std::queue<int> q; 
// 		q.push(from_vert); 
// 		visited_vert_set.insert(from_vert);
// 
// 		while(!q.empty()){
// 			int cur_vtx = q.front(); 
// 		}




		double shortest_dist = std::numeric_limits<double>::infinity();
		nearest_vid = -1;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		for(size_t k=0; k<path.size(); ++k){
			int vid = path[k];
			std::vector<int> s_path;
			p_mesh->m_BasicOp.GetShortestPath(from_vert, vid, s_path);
			double cur_dist=0;
			if(s_path.size() > 1){
				for(size_t i=1; i<s_path.size(); ++i){
					cur_dist += (vCoord[s_path[i]] - vCoord[s_path[i-1]]).abs();
				}
			}
			if(shortest_dist > cur_dist){
				shortest_dist = cur_dist;
				nearest_vid = vid;
			}
		}
		assert(nearest_vid != -1);
		return shortest_dist;
	}
}
