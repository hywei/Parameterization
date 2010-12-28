#ifndef CHARTCREATOR_H_
#define CHARTCREATOR_H_

#include "Parameterization.h"
#include "ParamPatch.h"
#include "ParamChart.h"

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include <map>

class MeshModel;

namespace PARAM
{    
    class ChartCreator
    {
    public:
        ChartCreator(boost::shared_ptr<MeshModel> _p_mesh);
        ~ChartCreator();

        bool LoadPatchFile(const std::string& patch_file);

        bool FormParamCharts();

    public:
        //! Get Methods
        
        const std::vector<ParamPatch>& GetPatchArray() const { return m_patch_array; }
        const std::vector<ParamChart>& GetChartArray() const { return m_chart_array; }
        const std::vector<PatchConner>& GetPatchConnerArray() const { return m_patch_conner_array; }
        const std::vector<PatchEdge>& GetPatchEdgeArray() const { return m_patch_edge_array; }
        int GetPatchNumber() const { return (int) m_patch_array.size(); }
        int GetChartNumber() const { return (int) m_chart_array.size(); }

    private:

        //! set each patch's conners
        void SetPatchConners();
        
        //! set each patch's neighbor info
        void SetPatchNeighbors();
        

        void SetMeshEdgePatchEdgeMapping(std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping) const;
        bool FloodFillFaceForAllPatchs();
        bool FloodFillFaceAPatch(int init_fid, std::vector<bool>& face_visited_flag, const std::map< std::pair<int, int>, std::vector<int> >& me_pe_mapping);

        std::vector<int> GetMeshEdgeAdjFaces(int, int) const;
    private:
        boost::shared_ptr<MeshModel> p_mesh;
        
        std::vector<ParamPatch> m_patch_array;

        std::vector<PatchConner> m_patch_conner_array;
        std::vector<PatchEdge> m_patch_edge_array;
        std::vector<ParamChart> m_chart_array;
        
    };
}

#endif //CHARTCREATOR_H_