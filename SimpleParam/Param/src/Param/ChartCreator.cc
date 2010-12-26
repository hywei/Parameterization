#include "ChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>
#include <iostream>
#include <queue>

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
        return true;
    }

    bool ChartCreator::FormParamCharts()
    {
        //! this function should be call after call LoadPatchFile(const std::string&)
        SetPatchConners();
        SetPatchNeighbors();

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
}
