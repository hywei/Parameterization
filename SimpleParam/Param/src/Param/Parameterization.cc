#include "Parameterization.h"
#include "../ModelMesh/MeshModel.h"
#include "../Numerical/MeshSparseMatrix.h"

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
}
