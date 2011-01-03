#include "Parameterization.h"
#include "../ModelMesh/MeshModel.h"
#include "../Numerical/MeshSparseMatrix.h"
#include "../Common/HSVColor.h"
#include <limits>
#include <set>
#include <map>
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

		const PolyIndexArray& adj_vertices = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();

		std::set<int> path_vert_set;
		if(path.size() <= 2){
			for(size_t k=0; k<path.size(); ++k) path_vert_set.insert(path[k]);
		}else{
			for(size_t k=1; k<path.size()-1; ++k) path_vert_set.insert(path[k]);
		}
		
		std::set<int> visited_vert_set;
		std::queue<int> q; 
		int cur_step =0;
		q.push(from_vert); 
		q.push(cur_step);
		visited_vert_set.insert(from_vert);

		int cur_vert;
		while(!q.empty()){
			cur_vert = q.front(); q.pop();
			cur_step = q.front(); q.pop();

			if(path_vert_set.find(cur_vert) != path_vert_set.end()) { break; }

			const IndexArray& adj_vert_array = adj_vertices[cur_vert];

			for(size_t i=0; i<adj_vert_array.size(); ++i){
				int adj_vert = adj_vert_array[i];
				if(visited_vert_set.find(adj_vert) == visited_vert_set.end()){
					q.push(adj_vert);
					q.push(cur_step+1);
					visited_vert_set.insert(adj_vert);
				}
			}

		}

		nearest_vid = cur_vert;
		return cur_step*1.0;
	}

	/************************************************************************/
	/* \fn int FindShortestPathInRegion(int start_vid, int end_vid, 
	const std::set<size_t>& regoin_face_set, PATH& path)
	*  \brief FindShortestPathInRegion
	*  Find the shortest mesh surface path between two vertexes and this path's
	*  edges should be in one fixed region
	*  @param start_vid The start vertex id
	*  @param end_vid The end vertex id
	*  @param regoin_edge_set The edge set  make up this region
	*  @param path The path we find
	*  @return return 0 if success else return Error_Code
	/************************************************************************/
	int FindShortestPathInRegion(boost::shared_ptr<MeshModel> p_mesh, int start_vid, int end_vid, 
		const std::set< std::pair<int, int> >& region_edge_set, std::vector<int>& path)
	{
		assert(p_mesh);
		path.clear();
		if(start_vid==end_vid) {
			path.push_back(start_vid); return 0;
		}

		PolyIndexArray& adjVtxArray = p_mesh->m_Kernel.GetVertexInfo().GetAdjVertices();
		CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		const std::set<std::pair<int, int> >& edge_set = region_edge_set;

		std::map<int,int> preVtx; //
		preVtx[end_vid] = -1;
		std::map<int,double> dist; //
		dist[end_vid] = 0.0;
		std::set<int> vtxInQueue; // vertexes in queue
		queue<int> q;
		q.push(end_vid);
		vtxInQueue.insert(end_vid);

		double curLen;
		map<int,double>::iterator im;
		while(!q.empty())
		{
			int u = q.front(); q.pop();
			vtxInQueue.erase( vtxInQueue.find(u));
			if(u==(int)start_vid) break;

			double len = dist[u];
			IndexArray& adjVtxes = adjVtxArray[u];
			int nv = (int)adjVtxes.size();
			for(int k=0; k<nv; ++k)
			{
				int vid = adjVtxes[k];
				std::pair<int, int> e = MakeEdge(u, vid);
				if(edge_set.find(e) == edge_set.end()) continue;

				im = dist.find(vid);
				if(im == dist.end()) curLen=INFINITE_DISTANCE;
				else curLen = im->second;
				double d = (vCoord[vid] -vCoord[u]).abs();
				if( len + d < curLen){
					dist[vid] = len+d;
					preVtx[vid] = u;
					if(vtxInQueue.find(vid) == vtxInQueue.end()){
						vtxInQueue.insert(vid);
						q.push(vid);
					}
				}
			}
		}

		int curVID = start_vid, preVID = preVtx[curVID];
		while(preVID!= -1){
			path.push_back(curVID);
			curVID = preVID;
			preVID = preVtx[curVID];
		}
		path.push_back(end_vid);

		return 0;
	}

	HalfEdge::HalfEdge(){}

	HalfEdge::~HalfEdge()
	{
		for(std::vector<_HE_edge*>::iterator it = hf_edge.begin(); it != hf_edge.end(); ++it)
		{
			if( (*it) != NULL) delete (*it);
		}
	}

	int HalfEdge::CreateHalfEdge(boost::shared_ptr<MeshModel> mesh)
	{
		assert(mesh);
		if(mesh == NULL) return -1;

		PolyIndexArray& fIndex = mesh->m_Kernel.GetFaceInfo().GetIndex();

		size_t faceNum = fIndex.size();
		for(size_t k=0; k<faceNum; ++k)
		{
			IndexArray& f = fIndex[k];
			size_t vtxNum = f.size();

			size_t vid1 , vid2;

			_HE_edge* prev_hf = NULL;
			_HE_edge* first_hf = NULL;

			for(size_t i=0; i<vtxNum; ++i)
			{
				vid1 = f[i];
				vid2 = f[(i+1)%vtxNum];

				_HE_edge* edge = new _HE_edge(k, vid1);
				if(i==0) first_hf = edge;

				hf_edge.push_back(edge);
				size_t edgeID = hf_edge.size() - 1;

				std::pair<size_t, size_t> e = make_pair(vid2, vid1);
				std::map< std::pair<size_t, size_t>, size_t>::iterator im = edge_map.find(e); 	

				if(im != edge_map.end())
				{
					// the pair half edge exsits
					edge->pair = hf_edge[im->second];
					hf_edge[im->second]->pair = edge;
				}
				edge_map[std::make_pair(vid1, vid2)] = edgeID;
				if(prev_hf != NULL) prev_hf->next = edge;
				if(i==vtxNum-1) edge->next = first_hf;

				prev_hf = edge;
			}

		}

		return 0;
	}

	_HE_edge* HalfEdge::GetHalfEdge(std::pair<size_t, size_t> e)
	{
		std::map<std::pair<size_t, size_t>, size_t>::iterator im = edge_map.find(e);
		if(im == edge_map.end()) return NULL;
		size_t edgeId = im->second;
		assert(hf_edge[edgeId] != NULL);
		return hf_edge[edgeId];
	}

	const _HE_edge* HalfEdge::GetHalfEdge(std::pair<size_t, size_t> e) const
	{
		std::map<std::pair<size_t, size_t>, size_t>::const_iterator im = edge_map.find(e);
		if(im == edge_map.end()) return NULL;
		size_t edgeId = im->second;
		assert(hf_edge[edgeId] != NULL);
		return hf_edge[edgeId];
	}

	std::vector<int> GetMeshEdgeAdjFaces(boost::shared_ptr<MeshModel> p_mesh, int vtx1, int vtx2)
	{
		std::vector<int> adj_faces;
		if(p_mesh == NULL)
			return adj_faces;

		const PolyIndexArray& vf_adj_array = p_mesh->m_Kernel.GetVertexInfo().GetAdjFaces();
		const IndexArray& adj_faces1 = vf_adj_array[vtx1];
		const IndexArray& adj_faces2 = vf_adj_array[vtx2];

		for(size_t k=0; k<adj_faces1.size(); ++k)
		{
			int fid = adj_faces1[k];
			if(find(adj_faces2.begin(), adj_faces2.end(), fid) != adj_faces2.end())
			{
				adj_faces.push_back(fid);
			}
		}

		return adj_faces;
	}

	/************************************************************************/
	/** \brief FindInnerFace
	* Find the inner face of one region which is surrounded 
	* by the boundary paths
	* @param boundary_path The boundaries that surround this region, the paths
	* should have been sorted by cw and be a circle
	* @param face_set The faces that this region includes
	* @return Return 0 if there is no error happen, else return Error Code
	/************************************************************************/
	int FindInnerFace(boost::shared_ptr<MeshModel> p_mesh, const std::vector<int>& boundary_path,  
		std::vector<int>& face_set, const HalfEdge& half_edge)
	{
	    assert(p_mesh);
		size_t faceNum = p_mesh->m_Kernel.GetFaceInfo().GetIndex().size();
		face_set.clear();
		vector<bool> faceVisitFlag(faceNum, false);
		// find the boundary edges
		std::set<std::pair<size_t, size_t> > edges;
		size_t bdVtxNum = boundary_path.size();
		size_t vid1, vid2;
		for(size_t k=1; k<bdVtxNum; ++k)
		{
			vid1 = boundary_path[k-1];
			vid2 = boundary_path[k];
			edges.insert(std::make_pair(vid1, vid2));
		}

		PolyIndexArray& fIndex = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		for(size_t k=1; k<bdVtxNum; ++k)
		{
			vid1 = boundary_path[k-1];
			vid2 = boundary_path[k];

			if(p_mesh->m_BasicOp.IsBoundaryVertex(vid1) || p_mesh->m_BasicOp.IsBoundaryFace(vid2))
			{
				continue;
			}

			std::pair<size_t, size_t> e = std::make_pair(vid1, vid2);
			std::pair<size_t, size_t> e_ = std::make_pair(vid2, vid1);
			if(edges.find(e_) == edges.end())
			{
				// the msc edge is cw, but the half edge is ccw
				const HE_edge* hf_e = half_edge.GetHalfEdge(e_);
				assert(hf_e != NULL);
				if(hf_e == NULL) {printf("he error!\n"); return -1; }

				size_t fid = hf_e->fid;

				if(faceVisitFlag[fid] == true) continue;

				queue<size_t> q;
				q.push(fid);
				faceVisitFlag[fid] = true;

				while(!q.empty())
				{
					fid = q.front(); q.pop();
					face_set.push_back(fid);

					for(size_t k=0; k<3; ++k)
					{
						size_t vid1 = fIndex[fid][k];
						size_t vid2 = fIndex[fid][(k+1)%3];

						if(edges.find(std::make_pair(vid1, vid2)) == edges.end() &&
							edges.find(std::make_pair(vid2, vid1)) == edges.end())
						{
							std::pair<size_t, size_t> e = MakeEdge(vid1, vid2);
							std::vector<int> adjFaces= GetMeshEdgeAdjFaces(p_mesh, e.first, e.second);
							if(adjFaces[0] == adjFaces[1]) continue;
							size_t nxtF = (fid==adjFaces[0] && adjFaces.size() == 2)?adjFaces[1]:adjFaces[0];
							if(faceVisitFlag[nxtF] == false)
							{
								faceVisitFlag[nxtF] = true;
								q.push(nxtF);
							}
						}

					}// end for
				}// end while
			}
		}

		return 0;
	}
}
