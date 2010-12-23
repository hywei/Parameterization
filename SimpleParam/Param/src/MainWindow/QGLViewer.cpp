#include <iostream>
#include <assert.h>

// --------------------

#include <QMouseEvent>
#include "QGLViewer.h"
#include <QGLFramebufferObject>

#include <Windows.h>
#include <string>
#include <fstream>
// --------------------

#include "../ModelMesh/MeshModel.h"
#include "../OpenGL/GLContext.h"
#include "../UI/UIHandler.h"
#include "../Common/Utility.h"
#include "../Param/QuadParameter.h"
#include "../Param/ParamDrawer.h"
using namespace Qt;


//== IMPLEMENTATION ==========================================================


//----------------------------------------------------------------------------

QGLViewer::QGLViewer( QWidget* _parent )
: QGLWidget( _parent )
{
	init();	
}

QGLViewer::~QGLViewer()
{
}

//----------------------------------------------------------------------------

QGLViewer::QGLViewer( QGLFormat& _fmt, QWidget* _parent )
: QGLWidget( _fmt, _parent )
{
	init();

}


//----------------------------------------------------------------------------

void QGLViewer::init(void)
{
	// qt stuff
	setAttribute(Qt::WA_NoSystemBackground, true);
	setFocusPolicy(Qt::StrongFocus);
	setAcceptDrops( true );
	setCursor(PointingHandCursor);

	createPopMenu();
	
	p_mesh = boost::shared_ptr<MeshModel> (new MeshModel);

	p_opengl = boost::shared_ptr<COpenGL> (new COpenGL);
	p_opengl->OnCreate(this);
	p_opengl->OnInit();
	
	p_UIHander = boost::shared_ptr<CUIHandler> (new CUIHandler);
	p_UIHander->Init(this);
	p_UIHander->GetMouseMode() = MOUSE_MODE_SPIN;	

	p_util = boost::shared_ptr<Utility> (new Utility);

	p_param = boost::shared_ptr<PARAM::QuadParameter> (new PARAM::QuadParameter(p_mesh));

	p_param_drawer = boost::shared_ptr<PARAM::ParamDrawer> (new 
		PARAM::ParamDrawer(*p_param.get()));
}


void QGLViewer::createPopMenu()
{

	// popup menu
	popup_menu_ = new QMenu(this);

	QActionGroup* renderAction = new QActionGroup(this);
	QAction *a ;

	a = popup_menu_->addAction(tr("Solid Smooth"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingSolidSmooth()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("Solid Flat"));
	a->setCheckable(true);
	a->setChecked(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingSolidFlat()));
	renderAction->addAction(a);


	a = popup_menu_->addAction(tr("Solid Wireframe"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingSolidWireframe()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("Transparent"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingTransparent()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("Wireframe"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingWireframe()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("Vertices"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingVertices()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("Texture"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingTexture()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("Parameter Texture"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingParamTexture()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("VertexColor"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingVertexColor()));
	renderAction->addAction(a);

	a = popup_menu_->addAction(tr("FaceColor"));
	a->setCheckable(true);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingFaceColor()));
	renderAction->addAction(a);

	popup_menu_->addSeparator();
	a = new QAction(tr("Bounding Box"), this);
	connect(a, SIGNAL(triggered()), this, SLOT(RenderingBoundingBox()));
	popup_menu_->addAction(a);
	a->setCheckable(true);
	a->setChecked(true);
}

// ----------------- public slots ------------------
void QGLViewer::loadMeshModel()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open mesh file"),
		tr(""),
		tr("OBJ Files (*.obj);;"
		"OFF Files (*.off);;"
		"STL Files (*.stl);;"
		"All Files (*)"));
	std::string f = std::string((const char *)fileName.toLocal8Bit());
	if(f.size()!=0)
	{
		// add your codes here
		p_mesh->ClearData();
		p_mesh->AttachModel(f);
	 	 		
		Coord center;
		double radius;
		p_mesh->m_Kernel.GetModelInfo().GetBoundingSphere(center, radius);
		p_opengl->OnInit();
		p_opengl->SetObjectInfo(center[0], center[1], center[2], radius);
		p_opengl->OnResize(this->width(), this->height());



// 		p_param->SetMeshTexCoordForGrid16();

		p_opengl->Create2DTexture(1);

		updateGL();
	}
}

int QGLViewer::loadTextureImage()
{
	QString fileName = QFileDialog::getOpenFileName(
		this, 
		tr("Open Texture Image File"),
		tr(""),
		tr("Image file(*.bmp);;"));

	std::string f = std::string((const char *)fileName.toLocal8Bit());
	if(f.size() != 0)
	{
		if(p_mesh )
		{
			p_opengl->OnBeginPaint();
			p_mesh->CreateTexture(f);
			p_opengl->OnEndPaint();

			updateGL();
		}else
			return -1;
	}else{
		return -1;
	}
	return 0;
}

int QGLViewer::loadQuadFile()
{
	QString fileName = QFileDialog::getOpenFileName(
		this, 
		tr("Open Quad File"), 
		tr(""),
		tr("Quad file(*.quad);;")
		);

	std::string f = std::string((const char *) fileName.toLocal8Bit());
	if(f.size() != 0)
	{
		if(p_mesh)
		{			
			p_param->LoadQuadFile(f);
			p_param->ComputeParamCoord();

		}else
			return -1;
	}else{
		return -1;
	}

	return 0;
}


void QGLViewer::CreateSquareTexture()
{
	if(p_opengl == NULL) return;
	p_opengl->Create2DTexture(1);
	updateGL();
}

void QGLViewer::CreateLineTexture()
{
	if(p_opengl == NULL) return;
	p_opengl->Create2DTexture(2);
	updateGL();
}

void QGLViewer::CreateBoundaryTexture()
{
	if(p_opengl == NULL) return;
	p_opengl->Create2DTexture(3);
	updateGL();
}

void QGLViewer::SetParamDrawerSelectVertMode()
{
	if(p_param_drawer == NULL) return;
	p_UIHander->GetMouseMode() = MOUSE_MODE_SELECT_VERTEX;
	p_param_drawer->SetDrawMode(PARAM::ParamDrawer::DrawMode::DRAWSELECTION);

	updateGL();
}

void QGLViewer::SetParamDrawerCorrespondMode()
{
	if(p_param_drawer == NULL) return;
	p_param_drawer->SetDrawMode(PARAM::ParamDrawer::DrawMode::DRAWSELECTION);
	updateGL();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void QGLViewer::RenderingSolidSmooth()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_SOLID_SMOOTH;
}

void QGLViewer::RenderingSolidFlat()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_SOLID_FLAT;
}

void QGLViewer::RenderingSolidWireframe()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_SOLID_AND_WIREFRAME;
}

void QGLViewer::RenderingTransparent()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_TRANSPARENT;
}

void QGLViewer::RenderingWireframe()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_WIREFRAME;
}

void QGLViewer::RenderingVertices()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_VERTICES;
}

void QGLViewer::RenderingVertexColor()
{
	// add your codes here
	p_mesh->m_Render.Mode() = RENDER_MODEL_VERTEX_COLOR;	
}

void QGLViewer::RenderingTexture()
{
	p_mesh->m_Render.Mode() = RENDER_MODEL_TEXTURE_MAPPING;
}

void QGLViewer::RenderingParamTexture()
{
	p_mesh->m_Render.Mode() = RENDER_MODEL_PARAM_TEXTURE;
}

void QGLViewer::RenderingFaceColor()
{
	p_mesh->m_Render.Mode() = RENDER_MODEL_FACE_COLOR;
}

void QGLViewer::RenderingBoundingBox()
{
	// add your codes here
	p_util->ToggleFlag(p_mesh->m_Render.State(), RENDER_MODEL_BOUNDING_BOX);
}

void QGLViewer::mouseSpin()
{
	// add your codes here
	p_UIHander->GetMouseMode() = MOUSE_MODE_SPIN;
	updateGL();
}

void QGLViewer::mouseMove()
{
	// add your codes here
	p_UIHander->GetMouseMode() = MOUSE_MODE_MOVE;
	updateGL();
}

void QGLViewer::mouseZoom()
{
	// add your codes here
	p_UIHander->GetMouseMode() = MOUSE_MODE_ZOOM;
	updateGL();;
}

void QGLViewer::getBoundingSphere(Coord& center, double& radius)
{
	p_mesh->m_Kernel.GetModelInfo().GetBoundingSphere(center, radius);
}


//----------------------------------------------------------------------------

void QGLViewer::setDefaultMaterial(void)
{
	GLfloat mat_a[] = {0.1, 0.1, 0.1, 1.0};
	GLfloat mat_d[] = {0.7, 0.7, 0.5, 1.0};
	GLfloat mat_s[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat shine[] = {120.0};

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
}


//----------------------------------------------------------------------------

void QGLViewer::setDefaultLight(void)
{
	GLfloat pos1[] = { 1.0,  1.0, -0.2, 0.0};
	GLfloat pos2[] = {-1.0,  1.0, -0.2, 0.0};
	GLfloat pos3[] = { 0.0,  0.0,  1.0,  0.0};


	GLfloat col1[] = { 0.7,  0.7,  0.8,  1.0};
	GLfloat col2[] = { 0.8,  0.7,  0.7,  1.0};
	GLfloat col3[] = { 1.0,  1.0,  1.0,  1.0};

	glEnable(GL_LIGHT0);    
	glLightfv(GL_LIGHT0,GL_POSITION, pos1);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
	glLightfv(GL_LIGHT0,GL_SPECULAR, col1);

	glEnable(GL_LIGHT1);  
	glLightfv(GL_LIGHT1,GL_POSITION, pos2);
	glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
	glLightfv(GL_LIGHT1,GL_SPECULAR, col2);

	glEnable(GL_LIGHT2);  
	glLightfv(GL_LIGHT2,GL_POSITION, pos3);
	glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
	glLightfv(GL_LIGHT2,GL_SPECULAR, col3);

	// local viewer
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	// two side light
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}


//----------------------------------------------------------------------------

void QGLViewer::initializeGL()
{  
	// OpenGL state
	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glDisable( GL_DITHER );
	glEnable( GL_DEPTH_TEST );

	// Material
	setDefaultMaterial();

	// Lighting
	glLoadIdentity();	
	setDefaultLight(); 

	// Fog
//	glEnable(GL_FOG);
//	GLfloat fogColor[4] = { 0.3, 0.3, 0.4, 1.0 };
	GLfloat fogColor[4] = { 0.0, 0.0, 0.25, 1.0};
	glFogi(GL_FOG_MODE,    GL_LINEAR);
	glFogfv(GL_FOG_COLOR,  fogColor);
	glFogf(GL_FOG_DENSITY, 0.35);
	glHint(GL_FOG_HINT,    GL_DONT_CARE);
	glFogf(GL_FOG_START,    5.0f);
	glFogf(GL_FOG_END,     25.0f);
   
	// scene pos and size
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}




void QGLViewer::resizeGL( int _w, int _h )
{
	assert(p_opengl!=NULL);
	p_opengl->OnResize(_w, _h);
}


//----------------------------------------------------------------------------


void QGLViewer::paintGL()
{
	p_opengl->OnBeginPaint();
	p_opengl->DetectOpenGLError();
	
	p_opengl->SetModelView();

	if(p_mesh) p_mesh->DrawModel();

//	if(p_param) p_param->Draw();
	if(p_param_drawer) p_param_drawer->Draw();

	p_opengl->DetectOpenGLError();

	p_opengl->OnEndPaint();
}


//----------------------------------------------------------------------------


void
QGLViewer::mousePressEvent( QMouseEvent* _event )
{
	int x = _event->x();
	int y = _event->y();
	switch(_event->button())
	{
	case RightButton:
		popup_menu_->exec(QCursor::pos());
		break;

	case LeftButton:
		assert(p_UIHander);
		p_UIHander->OnLButtonDown(_event->buttons(), x, y);
		break;

	default:
		break;
	}
	
	updateGL();
}


//----------------------------------------------------------------------------

void QGLViewer::mouseMoveEvent( QMouseEvent* _event)
{
	int x = _event->x();
	int y = _event->y();
	
	p_UIHander->OnMouseMove(_event->buttons(), x, y);
	updateGL();
	
}

//----------------------------------------------------------------------------


void QGLViewer::mouseReleaseEvent( QMouseEvent*  _event  )
{
	int x = _event->x();
	int y = _event->y();
	switch(_event->button())
	{
	case RightButton:
		popup_menu_->exec(QCursor::pos());
		break;

	case LeftButton:
		assert(p_UIHander);
		p_UIHander->OnLButtonUp(_event->buttons(), x, y);
		break;

	default:
		break;
	}
	updateGL();
}


//-----------------------------------------------------------------------------


void QGLViewer::wheelEvent(QWheelEvent* _event)
{
	int x = _event->x();
	int y = _event->y();
	p_UIHander->OnMouseWheel(_event->modifiers(), _event->delta(), x, y);

	updateGL();
}


//----------------------------------------------------------------------------
void QGLViewer::keyPressEvent( QKeyEvent* _event)
{
	switch( _event->key() )
	{
	case Key_Escape:
		qApp->quit();
		break;
	default: break;
	}
	p_UIHander->OnKeyDown(_event->key());
	updateGL();
}
