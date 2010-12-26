#ifndef PARAMDRAWER_H_
#define PARAMDRAWER_H_

#include "../Common/BasicDataType.h"

namespace PARAM
{
	class Parameter;
	class SurfaceCoord;

	class ParamDrawer
	{
	public:
		enum DrawMode
		{
			DRAWNOTHING = 0x00000000,

			//! draw patch info
			DRAWLAYOUT = 0x00000001,
			DRAWPATCHCONNER = 0x00000011, 
			DRAWPATCHEDGE = 0x000000021,
			DRAWPATCHFACE = 0x000000041,

			//! draw distortion info
			DRAWDISTORTION = 0x00000002,
			DRAWFACEHARMONICDISTORTION = 0x00000012,
			DRAWVERTEXHARMONICDISTORTION = 0x00000022,
			
			//! Draw texture
			DRAWFACETEXTURE = 0x00000004,

			//! Draw Corresponding
			DRAWCORRESPONGING = 0x00000008,
			DRAWSELECTION = 0x000000018
		};
	public:
		void SetUnCorrespondingVertArray(const std::vector<int>& uncorresponding_vert_array)
		{
			m_uncorrespondnig_vert_array = uncorresponding_vert_array;
		}

		void SetSelectedVertCoord(const Coord& select_coord);
		void SetSelectedVertCoord(const SurfaceCoord& select_coord);
			

		int FindSelectedVertId(const Coord& select_coord);

		void SetDrawMode(DrawMode mode)
		{
			m_draw_mode = mode;
		}

		int GetSelectedVertID() const { return m_selected_vert_id; }

	public:
		ParamDrawer(const Parameter& quad_param);
		~ParamDrawer();

		void Draw() const;

	private:
		void DrawPatchConner() const;
		void DrawPatchEdge() const;
		void DrawPatchFace() const;

		void DrawOutRangeVertex() const;
		void DrawUnSetFace() const;

		void DrawFaceDistortion() const;
		void DrawFaceTexture() const;

		void DrawCorresponding() const;
		void DrawUnCorrespondingVertex() const;

		void DrawSelectedVert() const;

		void DrawSphere(const Coord& center, double point_size = 1.0) const;

	private:
		const Parameter& m_parameter;

		std::vector<int> m_uncorrespondnig_vert_array;

		int m_selected_vert_id;
		Coord m_selected_vert_coord;

		DrawMode m_draw_mode;
	};
}

#endif // PARAMDRAWER_H_
