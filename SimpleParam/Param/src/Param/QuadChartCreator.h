#ifndef CHART_DIVIDER_H_
#define CHART_DIVIDER_H_

#include "Parameterization.h"
#include "QuadPatch.h"
#include "QuadChart.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <boost/shared_ptr.hpp>

class MeshModel;

namespace PARAM
{
	class QuadChartCreator
	{
	public:
		QuadChartCreator(boost::shared_ptr<MeshModel> _p_mesh);
		~QuadChartCreator();

		bool LoadQuadFile(const std::string& quad_file);

		//! form the parameter charts
		bool FormParamQuadCharts();


	public:
		const std::vector<QuadPatch>& GetQuadPatchArray() const { return m_quad_patch_array; }
		const std::vector<QuadChart>& GetQuadChartArray() const { return m_quad_chart_array; }
		const std::vector<PatchConner>& GetPatchConnerArray() const { return m_patch_conner_array; }
		const std::vector<PatchEdge>& GetPatchEdgeArray() const { return m_patch_edge_array; }
		
		int GetPatchNumber() const { return m_quad_patch_array.size(); }
		int GetChartNumber() const { return m_quad_chart_array.size(); }

	private:
		//! find each patch's face set by flood fill
		bool FloodFillFaceForAllPatchs(); 

		//! set the mapping between the mesh edge and patch edge
		void SetMeshEdgePatchEdgeMapping(std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping) const;

		//! flood fill to find a patch's inner faces by a initial fill face
		bool FloodFillFaceAPatch(int init_fid, std::vector<bool>& face_visited_flag, 
			const std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping);
			

		//! set each patch edge's neighbor patch
		void SetPatchEdgeNeighborPatch();

		//! set each quad patch's neighbor patch
		void SetQuadPatchNeighborPatch();

		//! get a mesh edge's adjacent faces
		std::vector<int> GetMeshEdgeAdjFaces(int vid1, int vid2) const;
	private:
		boost::shared_ptr<MeshModel> p_mesh;

		std::vector<PatchConner> m_patch_conner_array;
		std::vector<PatchEdge> m_patch_edge_array;
		std::vector<QuadPatch> m_quad_patch_array;

		std::vector<QuadChart> m_quad_chart_array;
	};
}

#endif