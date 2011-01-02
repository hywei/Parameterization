#include "ChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>
#include <iostream>
#include <queue>
#include <set>

namespace PARAM
{
    ChartCreator::ChartCreator(boost::shared_ptr<MeshModel> _p_mesh):
        p_mesh(_p_mesh) {}

    ChartCreator::~ChartCreator(){}

    bool ChartCreator::LoadPatchFile(const std::string& patch_file)
    {
        std::ifstream fin(patch_file.c_str());
        if(fin.fail()){
            std::cerr << "@@@Eroor : Load Patch File Fail!" << std::endl;
        }

        //! clear data
        m_patch_conner_array.clear();
        m_patch_edge_array.clear();
        m_patch_array.clear();
        
        int patch_conner_num(0), patch_edge_num(0), patch_num(0);

        //! read patch conner info
        fin >> patch_conner_num;
        m_patch_conner_array.resize(patch_conner_num);

        for(int k=0; k<patch_conner_num; ++k){
            PatchConner& patch_conner = m_patch_conner_array[k];
            fin >> patch_conner.m_conner_type >> patch_conner.m_mesh_index;
            //! read neighbor info
            int neighbor_num;
            fin >> neighbor_num;
            patch_conner.m_nb_conner_index_array.resize(neighbor_num);
            patch_conner.m_nb_edge_index_array.resize(neighbor_num);
            for(int i=0; i<neighbor_num; ++i){
                fin >> patch_conner.m_nb_conner_index_array[i];
                fin >> patch_conner.m_nb_edge_index_array[i];
            }
            
        }    

        //! read patch edge info
        fin >> patch_edge_num;
        m_patch_edge_array.resize(patch_edge_num);

        for(int k=0; k<patch_edge_num; ++k){
            PatchEdge& patch_edge = m_patch_edge_array[k];
            fin >> patch_edge.m_conner_pair_index.first >> patch_edge.m_conner_pair_index.second;
            
            int neighbor_num;
            fin >> neighbor_num;
            patch_edge.m_nb_patch_index_array.resize(neighbor_num);
            for(int i=0; i<neighbor_num; ++i){
                fin >> patch_edge.m_nb_patch_index_array[i];
            }

            int path_vert_num;
            fin >> path_vert_num;
            patch_edge.m_mesh_path.resize(path_vert_num);
            for(int i=0; i<path_vert_num; ++i){
                fin >> patch_edge.m_mesh_path[i];
            }
        }

        //! read patch info
        fin >> patch_num;
        m_patch_array.resize(patch_num);

        for(int k=0; k<patch_num; ++k){
            ParamPatch& patch = m_patch_array[k];
            //! read face info
            int face_num;
            fin >> face_num;
            patch.m_face_index_array.resize(face_num);
            for(int i=0; i<face_num; ++i){
                fin >> patch.m_face_index_array[i];
            }

            //! read edge info
            int edge_num;
            fin >> edge_num;
            patch.m_edge_index_array.resize(edge_num);
            for(int i=0; i<edge_num; ++i){
                fin >> patch.m_edge_index_array[i];
            }                
        }
        fin.close();
        return true;
    }

    bool ChartCreator::FormParamCharts()
    {
        //! this function should be call after call LoadPatchFile(const std::string&)
        SetPatchConners();
        SetPatchNeighbors();

        // std::vector< std::pair<int, int> > pe_pair;
        // pe_pair.push_back( std::make_pair(0, 2));
        // pe_pair.push_back( std::make_pair(1, 6));
        // pe_pair.push_back( std::make_pair(5, 7));
        // pe_pair.push_back( std::make_pair(3, 4));
        // pe_pair.push_back( std::make_pair(2, 7));
        // pe_pair.push_back( std::make_pair(1, 4));


        // ofstream output("add_edge.txt");
        // for(size_t k=0; k<pe_pair.size(); ++k){
        //     std::vector<int> new_path;
        //     p_mesh->m_BasicOp.GetShortestPath(pe_pair[k].first,pe_pair[k].second, new_path);
        //     output << new_path.size() <<" : " << std::endl;
        //     for(size_t i=0; i<new_path.size(); ++i){
        //         output << new_path[i] <<" ";
        //     }
        //     output << std::endl;
        // }
        // output.close();
        // return true;

        
        //FloodFillFaceForAllPatchs();

		ofstream face_out("faces.txt");
		for(size_t k=0; k<m_patch_array.size(); ++k){
			const ParamPatch& param_patch = m_patch_array[k];
			const std::vector<int>& face_vec = param_patch.m_face_index_array;
			face_out << "Patch " << k <<" : " << face_vec.size() << std::endl;
			for(size_t i=0; i<face_vec.size(); ++i){
				face_out << face_vec[i] <<" ";
			}
			face_out << std::endl;
		}



        ofstream fout("temp_file.txt");
        for(size_t k=0; k<m_patch_array.size(); ++k){
            const ParamPatch& param_patch = m_patch_array[k];
            const std::vector<int>& conners = param_patch.m_conner_index_array;
            const std::vector<int>& edges = param_patch.m_edge_index_array;
            const std::vector<int>& neighbors = param_patch.m_nb_patch_index_array;
            fout << "Patch " <<  k << " : " << std::endl;
                        
            for(size_t i=0; i<conners.size(); ++i){
                fout << conners[i] << " ";
            }
            fout << std::endl;

            for(size_t i=0; i<edges.size(); ++i){
                fout << edges[i] << " ";
            }
            fout << std::endl;
            
            for(size_t i=0; i<neighbors.size(); ++i){
                fout << neighbors[i] << " ";
            }
            fout << std::endl;
        }
        fout.close();
            
        
        m_chart_array.clear();
        m_chart_array.resize(m_patch_array.size());

        for(size_t k=0; k<m_chart_array.size(); ++k){
            const ParamPatch& patch = m_patch_array[k];
            ParamChart& chart = m_chart_array[k];
            chart.m_patch_id = k;
            chart.m_conner_param_coord_array.resize(patch.m_conner_index_array.size());
            chart.m_conner_param_coord_array[0] = ParamCoord(0, 0);
            chart.m_conner_param_coord_array[1] = ParamCoord(1, 0);
            if(chart.m_conner_param_coord_array.size() == 3){
                chart.m_conner_param_coord_array[2] = ParamCoord(0.5, sqrt(3.0)/2);
            }else if(chart.m_conner_param_coord_array.size() == 4){
                chart.m_conner_param_coord_array[2] = ParamCoord(1, 1);
                chart.m_conner_param_coord_array[3] = ParamCoord(0, 1);
            }                
        }

        std::cout << "Form " << m_chart_array.size() <<" Charts!" << std::endl;
        return true;
    }

    void ChartCreator::SetPatchConners()
    {
        //! we assume thate a patch's edges have been sorted by ccw
        //! then we can get conners with the same order

        for(size_t k=0; k<m_patch_array.size(); ++k){
            ParamPatch& patch = m_patch_array[k];
            size_t edge_num = patch.m_edge_index_array.size();
            patch.m_conner_index_array.resize(edge_num);

            std::vector<int> un_decide_conner;
            for(size_t i=0; i<patch.m_edge_index_array.size(); ++i){
                int cur_idx = patch.m_edge_index_array[i];
                int nxt_idx = patch.m_edge_index_array[(i+1)%edge_num];
                const PatchEdge& cur_edge = m_patch_edge_array[cur_idx];
                const PatchEdge& nxt_edge = m_patch_edge_array[nxt_idx];

                //! find the current edge's corresponding conner
                std::vector<int> conners_1(2), conners_2(2);
                conners_1[0] = cur_edge.m_conner_pair_index.first;
                conners_1[1] = cur_edge.m_conner_pair_index.second;
                conners_2[0] = nxt_edge.m_conner_pair_index.first;
                conners_2[1] = nxt_edge.m_conner_pair_index.second;

                if(conners_1[0]>conners_1[1]) swap(conners_1[0], conners_1[1]);
                if(conners_2[0]>conners_2[1]) swap(conners_2[0], conners_2[1]);

                if(conners_1 == conners_2){
                    std::cout << "Warning : these two adjcent edges' endponts are all the same!" << std::endl;
                    un_decide_conner.push_back(i);
                }else{
                    if(conners_1[0] == conners_2[0] || conners_1[0] == conners_2[1]) patch.m_conner_index_array[i] = conners_1[1];
                    else patch.m_conner_index_array[i] = conners_1[0];
                }
            }

            if(un_decide_conner.size() == patch.m_conner_index_array.size()){
                std::cerr << "Error: can't set conner for this patch " << k << std::endl;
            }else{
				if(un_decide_conner.size() !=0 )
					std::cout << "Warning: there are undecide conner for patch " << k << std::endl;
                //! TODO:
                for(size_t k=0; k<un_decide_conner.size(); ++k){
                    
                }
                    
            }
                
        }            
    }

    void ChartCreator::SetPatchNeighbors()
    {
        for(size_t k=0; k<m_patch_array.size(); ++k){
            ParamPatch& patch = m_patch_array[k];
            patch.m_nb_patch_index_array.clear();
            for(size_t i=0; i<patch.m_edge_index_array.size(); ++i){
                int edge_idx = patch.m_edge_index_array[i];
                const PatchEdge& patch_edge = m_patch_edge_array[edge_idx];
                if(patch_edge.m_nb_patch_index_array.size() == 2){
                    int nb_patch_idx = (patch_edge.m_nb_patch_index_array[0] == k) ? patch_edge.m_nb_patch_index_array[1] : patch_edge.m_nb_patch_index_array[0];
                    patch.m_nb_patch_index_array.push_back(nb_patch_idx);
                }
                    
            }
        }
        
    }
    
    void ChartCreator::SetMeshEdgePatchEdgeMapping(std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping) const
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

  	bool ChartCreator::FloodFillFaceForAllPatchs()
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


	bool ChartCreator::FloodFillFaceAPatch(int init_fid, std::vector<bool>& face_visited_flag,
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
		for(size_t pid = 0; pid < m_patch_array.size(); ++pid)
		{
			const ParamPatch& param_patch = m_patch_array[pid];
			bool flag = true;
			for(int k=0; k<param_patch.m_edge_index_array.size(); ++k)
			{
				int pe_idx = param_patch.m_edge_index_array[k];
				if(patch_edge_set.find(pe_idx) == patch_edge_set.end()) { flag = false; break;}
				
			}
			if(flag == true) { patch_id = (int)pid; break;}
		}
		if(patch_id == -1) 
		{
			cout<<"Error: Can't find valid patch!\n";
			return false;
		}else{
            std::cout << "Form patch " << patch_id << " faces" << std::endl;
            const ParamPatch& param_patch = m_patch_array[patch_id];
            const std::vector<int>& conners = param_patch.m_conner_index_array;
            for(size_t k=0; k<conners.size(); ++k){
                std::cout << conners[k] <<" ";
            }
            std::cout << std::endl;
                
        }

       
		ParamPatch& param_patch = m_patch_array[patch_id];
		param_patch.m_face_index_array.clear();
		for(std::set<int>::const_iterator is = visited_face_set.begin(); is!=visited_face_set.end(); ++is)
		{
			face_visited_flag[*is] = true;
			param_patch.m_face_index_array.push_back(*is);
		}
       
		return true;
	}

    	std::vector<int> ChartCreator::GetMeshEdgeAdjFaces(int vtx1, int vtx2) const
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
