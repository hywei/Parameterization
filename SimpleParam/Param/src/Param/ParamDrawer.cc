#include "ParamDrawer.h"
#include "ParamPatch.h"
#include "Parameter.h"
#include "ChartCreator.h"
#include "../ModelMesh/MeshModel.h"
#include <vector>
#include <limits>

#include <GL/gl.h>
#include <GL/glu.h>

namespace PARAM
{
	ParamDrawer::ParamDrawer(const Parameter& _param) 
		: m_parameter(_param) { m_selected_vert_coord = -1; }
	ParamDrawer::~ParamDrawer(){}

	void ParamDrawer::Draw() const
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_POLYGON_SMOOTH);
		glPolygonMode(GL_FRONT, GL_FILL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(2.0, 2.0);
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

		if(m_draw_mode & DRAWPATCHCONNER) DrawPatchConner();
		if(m_draw_mode & DRAWSELECTION) DrawSelectedVert();

		DrawPatchConner();
        //		DrawSelectedVert();
        // std::cout << "Select Vertex : " << m_selected_vert_id << " : "<<
        //     m_selected_vert_coord[0] << " " << m_selected_vert_coord[1] << ' '
        //           << m_selected_vert_coord[2] << std::endl;
		
		DrawUnCorrespondingVertex();	   


        //
        //        DrawPatchFace();
        
		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_POLYGON_SMOOTH);

		glDisable(GL_LIGHTING);
		glDepthFunc(GL_LEQUAL);


		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glLineWidth(3.0f);

		if(m_draw_mode & DRAWPATCHEDGE) DrawPatchEdge();
		DrawPatchEdge();

		glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glEnable(GL_LIGHTING);
		glDepthFunc(GL_LESS);
		glEnable(GL_LIGHTING);
	
	
	}

	void ParamDrawer::DrawPatchConner() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();

		if(p_mesh == NULL) return ;
		if(p_chart_creator == NULL) return;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		int colors[56][3] = 
		{
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};


		const std::vector<PatchConner>& conner_array = p_chart_creator->GetPatchConnerArray();
		for(size_t k=0; k<conner_array.size(); ++k)
		{            
			int c = k%56;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			int mesh_idx = conner_array[k].m_mesh_index;
			DrawSphere(vCoord[mesh_idx]);
		}

	}

	void ParamDrawer::DrawPatchEdge() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();

		if(p_mesh == NULL) return;
		if(p_chart_creator == NULL) return;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const std::vector<PatchEdge>& patch_edge_array = p_chart_creator->GetPatchEdgeArray();

		int colors[56][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};

		for(size_t k=0; k<patch_edge_array.size(); ++k)
		{
			int c = k%56;
			glColor3ub(colors[c][0], colors[c][1], colors[c][2]);
			const std::vector<int>& path = patch_edge_array[k].m_mesh_path;
			glBegin(GL_LINE_STRIP);
			for(size_t i=0; i<path.size(); ++i)
			{
				const Coord& vtxCoord = vCoord[path[i]];
				glVertex3d(vtxCoord[0], vtxCoord[1], vtxCoord[2]);
			}
			glEnd();
		}
	}

    void ParamDrawer::DrawPatchFace() const
    {
        boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
        boost::shared_ptr<ChartCreator> p_chart_creator = m_parameter.GetChartCreator();
        
        if(p_mesh == NULL) return;
        if(p_chart_creator == NULL) return;

        int colors[56][3] = 
		{			  
			{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 0}, 
			{ 255, 0 , 255}, {0, 255, 255}, {255, 255, 255}, {0, 0, 0},

			{255, 128, 128}, {0, 64, 128},  {255, 128, 192}, {128, 255, 128}, 
			{0, 255, 128}, {128, 255, 255}, {0, 128, 255}, {255, 128, 255}, 

			{255, 0, 0}, {255, 255, 0}, {128, 255, 0}, {0, 255, 64}, 
			{0, 255, 255}, {0, 128, 192}, {128, 128, 192}, {255, 0, 255}, 

			{128, 64, 64}, {255, 128, 64}, {0, 255, 0}, {0, 128, 128}, 
			{0, 64, 128}, {128, 128, 255}, {128, 0, 64}, {255, 0, 128}, 

			{128, 0, 0}, {0, 128, 0}, {0, 128, 64}, {255, 255, 128},
			{0, 0, 255}, {0, 0, 160}, {128, 0, 128}, {128, 0, 255}, 

			{64, 0, 0}, {128, 64, 0}, {0, 64, 0}, {0, 64, 64}, 
			{0, 0, 128}, {0, 0, 64}, {64, 0, 64}, {64, 0, 128}, 
		};

        
        const std::vector<ParamPatch>& patch_array = p_chart_creator->GetPatchArray();
        const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
        const NormalArray& vert_norm_array = p_mesh->m_Kernel.GetVertexInfo().GetNormal();
        const CoordArray& vert_coord_array = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

        glBegin(GL_TRIANGLES);
        for(size_t k=0; k<patch_array.size(); ++k){           
            const ParamPatch& patch = patch_array[k];
            const std::vector<int>& faces = patch.m_face_index_array;
            glColor3ub(colors[k%56][0], colors[k%56][1], colors[k%56][2]);
            for(size_t i=0; i<faces.size(); ++i){
                int fid = faces[i];
                const IndexArray& vertices = face_list_array[fid];
                for(int j=0; j<3; ++j){
                    const Coord& v = vert_coord_array[vertices[j]];
                    const Normal& n = vert_norm_array[vertices[j]];
                    glNormal3d(n[0], n[1], n[2]);
                    glVertex3d(v[0], v[1], v[2]);
                }
            }
        }
        glEnd();
                    
    }

	void ParamDrawer::DrawOutRangeVertex() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();		
		if(p_mesh == NULL) return ;
		const std::vector<int>& out_range_vert_array = m_parameter.GetOutRangeVertArray();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		glColor3ub(255, 0, 0);
		for(size_t k=0; k<out_range_vert_array.size(); ++k)
		{
			int vid = out_range_vert_array[k];
			const Coord& vtxCoord = vCoord[vid];
			DrawSphere(vtxCoord, 0.5);
		}
	}

	void ParamDrawer::DrawUnCorrespondingVertex() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		glColor3ub(128, 128, 0);
		for(size_t k=0; k<m_uncorrespondnig_vert_array.size(); ++k)
		{
			int vid = m_uncorrespondnig_vert_array[k];
			const Coord& vtxCoord = vCoord[vid];
			DrawSphere(vtxCoord, 1.0);
		}
	}

	void ParamDrawer::DrawUnSetFace() const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();		
		if(p_mesh == NULL) return ;
		const std::vector<int>& unset_face_array = m_parameter.GetUnSetFaceArray();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();

		glColor3ub(0, 0, 255);
		for(size_t k=0; k<unset_face_array.size(); ++k)
		{
			int fid = unset_face_array[k];
			
			const IndexArray& fIndex = face_list_array[fid];
			glBegin(GL_TRIANGLES);

			for(int j = 0; j < 3; ++ j)
			{
					int vID = fIndex[j];
					const Coord& v = vCoord[vID];
					glVertex3d(v[0], v[1], v[2]);
			}
			
			glEnd();
		}

	}

	void ParamDrawer::DrawSelectedVert() const
	{
		glColor3ub(64, 0, 128);
		DrawSphere(m_selected_vert_coord, 2.0);
	}

	int ParamDrawer::FindSelectedVertId(const Coord& select_coord)
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return -1;
		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		int vert_num = p_mesh->m_Kernel.GetModelInfo().GetVertexNum();

		double min_dist = numeric_limits<double>::infinity();
		
		int nearest_vert_id(-1);

		for(int k=0; k<vert_num; ++k)
		{
			const Coord& vtx_coord = vCoord[k];
			double cur_dist = (vtx_coord - select_coord).abs();
			if(min_dist > cur_dist)
			{
				min_dist = cur_dist;
				nearest_vert_id = k;
			}
		}

		return nearest_vert_id;
	}

	void ParamDrawer::SetSelectedVertCoord(const Coord& select_coord)
	{
		m_selected_vert_coord = select_coord;
		m_selected_vert_id = FindSelectedVertId(m_selected_vert_coord);
	}

	void ParamDrawer::SetSelectedVertCoord(const SurfaceCoord& select_coord)
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();
		if(p_mesh == NULL) return;

		const CoordArray& vCoord = p_mesh->m_Kernel.GetVertexInfo().GetCoord();
		const PolyIndexArray& face_list_array = p_mesh->m_Kernel.GetFaceInfo().GetIndex();

		const IndexArray& face_vert_array = face_list_array[select_coord.face_index];
		Coord pos_coord(0, 0, 0);
		for(int k=0; k<3; ++k)
		{
			int vid = face_vert_array[k];
			pos_coord += vCoord[vid] * (select_coord.barycentric[k]);
		}
		m_selected_vert_coord = pos_coord;

		m_selected_vert_id = FindSelectedVertId(m_selected_vert_coord);
	}

	void ParamDrawer::DrawSphere(const Coord& center, double point_size) const
	{
		boost::shared_ptr<MeshModel> p_mesh = m_parameter.GetMeshModel();

		double bRadius;
		Coord c;
		p_mesh->m_Kernel.GetModelInfo().GetBoundingSphere(c, bRadius);
		GLUquadricObj* qobj = gluNewQuadric();
		gluQuadricDrawStyle(qobj, GLU_FILL);
		gluQuadricNormals(qobj, GLU_SMOOTH);
		glTranslatef((float)center[0], (float)center[1], (float)center[2]);

		gluSphere(qobj, point_size *0.02 * bRadius , 15, 10);
		glTranslatef((float)-center[0], (float)-center[1], (float)-center[2]);
	}
}
