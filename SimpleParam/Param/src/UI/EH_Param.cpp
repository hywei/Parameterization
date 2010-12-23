#include "EH_Param.h"
#include "../MainWindow/QGLViewer.h"
#include "../OpenGL/GLElement.h"
#include "../Param/ParamDrawer.h"

EH_Param::EH_Param()
{
	this->m_MouseMode = MOUSE_MODE_UNDEFINED;
	this->m_KeyMode = KEY_MODE_UNDEFINED;

	m_ValidMouseMode.insert(MOUSE_MODE_SELECT_VERTEX);
}

EH_Param::~EH_Param()
{
	glViewer = NULL;
}

void EH_Param::Init(QGLViewer* _glViewer)
{
	glViewer = _glViewer;
}

bool EH_Param::OnLButtonDown(unsigned int nFlags, int point_x, int point_y)
{
	Coord hit_coord;
	switch(m_MouseMode)
	{
	case MOUSE_MODE_SELECT_VERTEX:
		glViewer->p_opengl->GetWorldCoord(point_x, point_y, hit_coord);
		glViewer->p_param_drawer->SetSelectedVertCoord(hit_coord);
		glViewer->EmitSelectVertexSignal();
		break;

	default:
		break;
	}
	return true;
}

bool EH_Param::OnMouseMove(unsigned int nFlags, int point_x, int point_y)
{
	Coord hit_coord;
	switch(m_MouseMode)
	{
	case MOUSE_MODE_SELECT_VERTEX:
		glViewer->p_opengl->GetWorldCoord(point_x, point_y, hit_coord);
		glViewer->p_param_drawer->SetSelectedVertCoord(hit_coord);
		glViewer->EmitSelectVertexSignal();
		break;

	default:
		break;
	}
	return true;
}

bool EH_Param::IsMouseNeedHandle(MouseMode mouse_mode)
{
	this->m_MouseMode = mouse_mode;
	return m_ValidMouseMode.find(mouse_mode) != m_ValidMouseMode.end();
}