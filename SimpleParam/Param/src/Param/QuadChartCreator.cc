#include "QuadChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>
#include <iostream>
#include <queue>

namespace PARAM
{
	QuadChartCreator::QuadChartCreator(boost::shared_ptr<MeshModel> _p_mesh): p_mesh(_p_mesh){}
	QuadChartCreator::~QuadChartCreator(){}

	bool QuadChartCreator::LoadQuadFile(const std::string& quad_file)
	{
		std::ifstream fin(quad_file.c_str());
		if(fin.fail())
		{
			std::cout<<"Error : Load Quad File Fail!" << std::endl;
			return false;
		}

		m_patch_conner_array.clear(); 
		m_patch_edge_array.clear();
		m_quad_patch_array.clear();

		size_t patch_conner_num, patch_edge_num, patch_num;

		//! read patch conner infomation
		fin >> patch_conner_num; 
		m_patch_conner_array.resize(patch_conner_num);

		for(size_t k=0; k<patch_conner_num; ++k)
		{
			fin >> m_patch_conner_array[k];
		}

		fin >> patch_edge_num;
		m_patch_edge_array.resize(patch_edge_num);

		for(size_t k=0; k<patch_edge_num; ++k)
		{
			std::vector<int>& path = m_patch_edge_array[k].m_mesh_path;
			size_t path_vert_num;
			fin >> path_vert_num;
			path.resize(path_vert_num);
			for(size_t i=0; i<path_vert_num; ++i) 
			{
				fin >> path[i];
			}
		}

		fin >> patch_num;
		m_quad_patch_array.resize(patch_num);

		for(size_t k=0; k<patch_num; ++k)
		{
			QuadPatch& patch = m_quad_patch_array[k];
			patch.m_conner_index_array.resize(4);
			for(int i=0; i<4; ++i) fin >> patch.m_conner_index_array[i];
			patch.m_edge_index_array.resize(4);
			for(int i=0; i<4; ++i) fin >> patch.m_edge_index_array[i];
		}

		fin.close();

		return true;
	}

	void QuadChartCreator::SetPatchEdgeNeighborPatch()
	{
		for(size_t pid = 0; pid < m_quad_patch_array.size(); ++pid)
		{
			const QuadPatch& quad_patch = m_quad_patch_array[pid];
			for(int k=0; k<4; ++k)
			{
				int pe_idx = quad_patch.m_edge_index_array[k];
				PatchEdge& patch_edge = m_patch_edge_array[pe_idx];
				patch_edge.m_neighbor_patch_array.push_back(pid);
			}
		}
	}

	void QuadChartCreator::SetQuadPatchNeighborPatch()
	{
		for(size_t pid = 0; pid < m_quad_patch_array.size(); ++pid)
		{
			QuadPatch& quad_patch = m_quad_patch_array[pid];
			quad_patch.m_neighbor_patch_array.clear();
			for(int k=0; k<4; ++k)
			{
				int pe_idx = quad_patch.m_edge_index_array[k];
				const PatchEdge& patch_edge = m_patch_edge_array[pe_idx];
				for(size_t i=0; i<patch_edge.m_neighbor_patch_array.size(); ++i)
				{
					int patch_id = patch_edge.m_neighbor_patch_array[i];
					if(patch_id != pid) quad_patch.m_neighbor_patch_array.push_back(patch_id);
				}
			}
		}
	}

	bool QuadChartCreator::FloodFillFaceForAllPatchs()
	{
		std::map< std::pair<int, int>, std::vector<int> > me_pe_mapping;
		SetMeshEdgePatchEdgeMapping(me_pe_mapping);

		int face_num = p_mesh->m_Kernel.GetModelInfo().GetFaceNum();

		std::vector<bool> face_visited_flag(face_num, false);

		for(int fid = 0; fid < face_num; ++fid)
		{
			if(!face_visited_flag[fid])
			{
				if(!FloodFillFaceAPatch(fid, face_visited_flag, me_pe_mapping)) return false;
			}
		}
		
		return true;
	}


	bool QuadChartCreator::FloodFillFaceAPatch(int init_fid, std::vector<bool>& face_visited_flag,
		const std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping)	
	{
		const PolyIndexArray& face_list = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		std::set<int> visited_face_set;
		std::queue<int> q;
		q.push(init_fid);
		visited_face_set.insert(init_fid);

		std::set<int> patch_edge_set;

		while(!q.empty())
		{
			int cur_fid = q.front(); q.pop();

			const IndexArray& face_vertices = face_list[cur_fid];
			for(int k=0; k<3; ++k)
			{
				int vid1 = face_vertices[k], vid2 = face_vertices[(k+1)%3];
				std::pair<int, int> mesh_edge = (vid1 < vid2) ? make_pair(vid1, vid2) : make_pair(vid2, vid1);
				const std::map< std::pair<int, int>, std::vector<int> >::const_iterator im = me_pe_mapping.find(mesh_edge);
				if(im == me_pe_mapping.end()) 
				{
					const vector<int>& adj_faces = GetMeshEdgeAdjFaces(vid1, vid2);
					for(size_t i=0; i<adj_faces.size(); ++i)
					{
						int adj_fid = adj_faces[i];
						if(adj_fid != cur_fid && visited_face_set.find(adj_fid) == visited_face_set.end())
						{
							q.push(adj_fid);
							visited_face_set.insert(adj_fid);
						}
					}
				}else
				{
					const std::vector<int>& patch_edge_index_array = im->second;
					for(size_t k=0; k<patch_edge_index_array.size(); ++k)
					{
						int pe_idx = patch_edge_index_array[k];
						patch_edge_set.insert(pe_idx);
					
					}
				}
			}
		}

		/// decide which patch
		//: TODO: possible performance problem, we should really use a hash table but now we'll use a linear search
		int patch_id(-1);
		for(size_t pid = 0; pid < m_quad_patch_array.size(); ++pid)
		{
			const QuadPatch& quad_patch = m_quad_patch_array[pid];
			bool flag = true;
			for(int k=0; k<4; ++k)
			{
				int pe_idx = quad_patch.m_edge_index_array[k];
				if(patch_edge_set.find(pe_idx) == patch_edge_set.end()) { flag = false; break;}
				
			}
			if(flag == true) { patch_id = (int)pid; break;}
		}
		if(patch_id == -1) 
		{
			cout<<"Error: Can't find valid patch!\n";
			return false;
		}

		QuadPatch& quad_patch = m_quad_patch_array[patch_id];
		quad_patch.m_face_index_array.clear();
		for(std::set<int>::const_iterator is = visited_face_set.begin(); is!=visited_face_set.end(); ++is)
		{
			face_visited_flag[*is] = true;
			quad_patch.m_face_index_array.push_back(*is);
		}
		return true;
	}

	void QuadChartCreator::SetMeshEdgePatchEdgeMapping(std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping) const
	{
		me_pe_mapping.clear();
		for(size_t pe_idx = 0; pe_idx < m_patch_edge_array.size(); ++pe_idx)
		{
			const std::vector<int>& mesh_path = m_patch_edge_array[pe_idx].m_mesh_path;
			assert(mesh_path.size() >= 2);
			for(size_t k = 1; k<mesh_path.size(); ++k)
			{
				int vid1 = mesh_path[k-1];
				int vid2 = mesh_path[k];
				std::pair<int, int> mesh_edge = (vid1 < vid2) ? make_pair(vid1 , vid2) : make_pair(vid2, vid1);

				std::vector<int>& pe_array = me_pe_mapping[mesh_edge];
				if(find(pe_array.begin(), pe_array.end(), pe_idx) == pe_array.end()) pe_array.push_back(pe_idx);
			}
		}
	}

	bool QuadChartCreator::FormParamQuadCharts()
	{
		std::cout<<"Form Quad Chart Layout...\n";
		SetPatchEdgeNeighborPatch();
		SetQuadPatchNeighborPatch();

		if(!FloodFillFaceForAllPatchs()) return false;

		m_quad_chart_array.clear();
		m_quad_chart_array.resize(m_quad_patch_array.size());
		for(size_t k=0; k<m_quad_patch_array.size(); ++k)
		{
			const QuadPatch& quad_patch = m_quad_patch_array[k];
			QuadChart& quad_chart = m_quad_chart_array[k];
			quad_chart.m_patch_id = k;
			quad_chart.m_conner_param_coord_array.clear();
			quad_chart.m_conner_param_coord_array.resize(4);
		    quad_chart.m_conner_param_coord_array[0] = ParamCoord(0, 0);
			quad_chart.m_conner_param_coord_array[1] = ParamCoord(1, 0);
			quad_chart.m_conner_param_coord_array[2] = ParamCoord(1, 1);
			quad_chart.m_conner_param_coord_array[3] = ParamCoord(0, 1);
		}

		std::cout<<"Form " << m_quad_chart_array.size() << " Charts!\n";
		return true;
	}

	std::vector<int> QuadChartCreator::GetMeshEdgeAdjFaces(int vtx1, int vtx2) const
	{
		vector<int> adj_faces;
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

}
