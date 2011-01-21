#ifndef PARAMPATCH_H_
#define PARAMPATCH_H_

#include <vector>

namespace PARAM
{
    class PatchConner
    {
    public:
        PatchConner(){}
        ~PatchConner() {}

    public:
        int m_conner_type; //! 0: minimum, 1: saddle, 2:maximum
        int m_mesh_index;
        std::vector<int> m_nb_conner_index_array; //! neighbor conners index
        std::vector<int> m_nb_edge_index_array; //! neighbor edges index
    };

    class PatchEdge
    {
    public:
        PatchEdge(){}
        ~PatchEdge(){}

    public:
        std::pair< int, int> m_conner_pair_index; //! two endpoint conners' index
        std::vector<int> m_mesh_path; //! the mesh vertices index
        std::vector<int> m_nb_patch_index_array; //! the neighbor patch index
        
    };

    class ParamPatch
    {
    public:
        ParamPatch(){}
        ~ParamPatch(){}

    public:
        std::vector<int> m_conner_index_array; //! index in patch conner array
        std::vector<int> m_edge_index_array; //! index in patch edge array
        std::vector<int> m_face_index_array; //! index in mesh face array
        std::vector<int> m_nb_patch_index_array; //! neighbor patch index
        
    };
}

#endif //PARAMPATCH_H_
