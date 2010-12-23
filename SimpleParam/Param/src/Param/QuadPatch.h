#ifndef QUADPATCH_H_
#define QUADPATCH_H_

namespace PARAM
{

	typedef int PatchConner;
	
	//! Patch Edge 
	class PatchEdge
	{
	public:
		PatchEdge(){}		
		~PatchEdge(){}

		std::pair<int, int> GetConnerPair() const
		{
			return std::make_pair(m_mesh_path[0], m_mesh_path[m_mesh_path.size()-1]);
		}
	public:
		std::vector<int> m_mesh_path; //! the mesh vertices index  
		std::vector<int> m_neighbor_patch_array; //! the neighbor patches, this number should be 2 or 1(bounary edge)

	};

	class QuadPatch
	{
	public:
		QuadPatch(){}
		~QuadPatch(){}
	    
	public:
		std::vector<int> m_conner_index_array;
		std::vector<int> m_edge_index_array;
		std::vector<int> m_face_index_array;
		std::vector<int> m_neighbor_patch_array;
	};
}

#endif // QUADPATCH_H_