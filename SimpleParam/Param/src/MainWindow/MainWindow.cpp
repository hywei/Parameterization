#include "MainWindow.h"

#include <QtGui>
#include <sstream>
#include <iostream>

#include "CrossParamControl.h"
#include "QGLViewer.h"

#include "../Param/Parameter.h"
#include "../Param/ParamDrawer.h"

MainWindow::MainWindow(QWidget *parent):QMainWindow(parent)
{

	setMinimumSize(20, 20);
	resize(1500, 900); 

	glViewer_1 = new QGLViewer;
	glViewer_2 = new QGLViewer;

	m_cross_param_control = new CrossParamControl(glViewer_1, glViewer_2, this);	

	createActions();
	createMenus();
	createToolBars();

	QList <int> sizes1;
	sizes1<<600<<600;

	QSplitter* splitter = new QSplitter(Qt::Horizontal);
	splitter->addWidget(glViewer_1);
	splitter->addWidget(glViewer_2);
	splitter->setSizes(sizes1);

	QList <int> sizes2;
	sizes2 << 1300 << 200;
	QSplitter* mainSplitter = new QSplitter(Qt::Horizontal);
	mainSplitter->addWidget(splitter);
	mainSplitter->addWidget(m_cross_param_control);
	mainSplitter->setSizes(sizes2);
	mainSplitter->setStretchFactor(0, 1);

	setCentralWidget(mainSplitter);

	setWindowTitle("Cross Parameter");

	connect(glViewer_1, SIGNAL(select_vertex()), this, SLOT(findCorrespondingOnB()));
    //	connect(glViewer_2, SIGNAL(select_vertex()), this, SLOT(findCorrespondingOnA()));
}

void MainWindow::setStyle(QStyle* style)
{
	QWidget::setStyle(style);
	if(m_cross_param_control != 0)
	{
		m_cross_param_control->setStyle(style);

		QList<QWidget*> widgets = qFindChildren<QWidget*> (m_cross_param_control);
		foreach(QWidget* w, widgets)
			w->setStyle(style);
	}
}

void MainWindow::mouseSpin()
{
	glViewer_1->mouseSpin();
	glViewer_1->setCursor(Qt::PointingHandCursor);

	glViewer_2->mouseSpin();
	glViewer_2->setCursor(Qt::PointingHandCursor);
}

void MainWindow::mouseMove()
{
	glViewer_1->mouseMove();
	glViewer_1->setCursor(Qt::ArrowCursor);

	glViewer_2->mouseMove();
	glViewer_2->setCursor(Qt::ArrowCursor);
}

void MainWindow::mouseZoom()
{
	glViewer_1->mouseZoom();
	glViewer_1->setCursor(Qt::SizeAllCursor);

	glViewer_2->mouseZoom();
	glViewer_2->setCursor(Qt::SizeAllCursor);
}

void MainWindow::findCorrespondingOnA()
{
// 	int vid = glViewer_2->p_param_drawer->GetSelectedVertID();
// 	if(vid == -1) return;
// 	int chart_id = glViewer_2->p_param->GetVertexChartID(vid);
// 	PARAM::ParamCoord param_coord = glViewer_2->p_param->GetVertexParamCoord(vid);
// 	PARAM::ChartParamCoord chart_param_coord(param_coord, chart_id);
// 	PARAM::SurfaceCoord surface_coord;
// 	glViewer_1->p_param->FindCorrespondingOnSurface(chart_param_coord, surface_coord);  
// 	glViewer_1->p_param_drawer->SetSelectedVertCoord(surface_coord);
// 	glViewer_1->updateGL();
	m_cross_param_control->FindCorrespondingOnA();
}

void MainWindow::findCorrespondingOnB()
{
// 	int vid = glViewer_1->p_param_drawer->GetSelectedVertID();
// 	if(vid == -1) return;
// 	int chart_id = glViewer_1->p_param->GetVertexChartID(vid);
// 	PARAM::ParamCoord param_coord = glViewer_1->p_param->GetVertexParamCoord(vid);
//     std::cout << "Surface A --- chart id : " << chart_id << ", param_coord " << param_coord.s_coord <<" " << param_coord.t_coord << std::endl; 
// 	PARAM::ChartParamCoord chart_param_coord(param_coord, chart_id);
// 	PARAM::SurfaceCoord surface_coord;
// 	glViewer_2->p_param->FindCorrespondingOnSurface(chart_param_coord, surface_coord);  
// 	glViewer_2->p_param_drawer->SetSelectedVertCoord(surface_coord);
// 
//     int vid_2 = glViewer_2->p_param_drawer->GetSelectedVertID();
//     std::cout << "Surface B --- chart id : " << chart_id << ", param_coord " << param_coord.s_coord <<" "<<param_coord.t_coord << std::endl;
// 	glViewer_2->updateGL();
	m_cross_param_control->FindCorrespondingOnB();
}

void MainWindow::createActions()
{
	// file menu actions
	openModelAct = new QAction(QIcon("../images/open.png"), tr("&Open Model"), this);
	openModelAct->setShortcut(QKeySequence::New);
	openModelAct->setStatusTip(tr("open a mesh model"));
	connect(openModelAct, SIGNAL(triggered()), this, SLOT(openModel()));

	openTextureFileAct = new QAction(QIcon("../images/open.png"), tr("&Open Texture Image"), this);
	connect(openTextureFileAct, SIGNAL(triggered()), this, SLOT(openTextureImage()));

	openQuadFileAct = new QAction(QIcon("../images/open.png"), tr("&Open Quad File"), this);
	connect(openQuadFileAct, SIGNAL(triggered()), this, SLOT(openQuadFile()));

	saveModelAct = new QAction(QIcon("../images/save.png"), tr("&Save Model"), this);
	saveModelAct->setShortcut(QKeySequence::Save);
	saveModelAct->setStatusTip(tr("save this mesh model"));
	connect(saveModelAct, SIGNAL(triggered()), this, SLOT(saveModel()));

	saveAsBmpAct = new QAction(tr("Save as BMP"), this);
	connect(saveAsBmpAct, SIGNAL(triggered()), this, SLOT(saveAsBmp()));

	recentFileAct = new QAction(tr("Recent Files"), this);
	connect(recentFileAct, SIGNAL(triggered()), this, SLOT(recentFiles()));

	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcuts(QKeySequence::Quit);
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));


	// toolbar menu actions
	toolBarAct = new QAction(tr("ToolBar"), this);
	connect(toolBarAct, SIGNAL(triggered()), this, SLOT(toobBar()));

	stateBarAct = new QAction(tr("StateBar"), this);
	connect(stateBarAct, SIGNAL(triggered()), this, SLOT(stateBar()));

	// quad menu actions
	
	// help menu actions
	aboutAct = new QAction(tr("About"), this);
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));



	// mouse actions
	mouseSpinAct = new QAction(QIcon("../images/rotate-left.png"), tr("Spin"), this);
	connect(mouseSpinAct, SIGNAL(triggered()), this, SLOT(mouseSpin()));

	mouseMoveAct = new QAction(QIcon("../images/move.png"), tr("Move"), this);
	connect(mouseMoveAct, SIGNAL(triggered()), this, SLOT(mouseMove()));

	mouseZoomAct = new QAction(QIcon("../images/zoom.png"), tr("Zoom"), this);
	connect(mouseZoomAct, SIGNAL(triggered()), this, SLOT(mouseZoom()));

}


void MainWindow::createMenus()
{
	// create file menu
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openModelAct);
	fileMenu->addAction(openTextureFileAct);
	fileMenu->addAction(openQuadFileAct);
	fileMenu->addAction(saveModelAct);
	fileMenu->addSeparator();
	fileMenu->addAction(saveAsBmpAct);
	fileMenu->addAction(recentFileAct);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);

	// create view menu
	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(toolBarAct);
	viewMenu->addAction(stateBarAct);
    
	// create help menu
	helpMenu = menuBar()->addMenu(tr("&Help"));
	helpMenu->addAction(aboutAct);
}

void MainWindow::createToolBars()
{
	fileToolBar = addToolBar(tr("File"));
	fileToolBar->addAction(openModelAct);
	fileToolBar->addAction(saveModelAct);

	mouseActToolBar = addToolBar(tr("Mouse Action"));
	mouseActToolBar->addAction(mouseSpinAct);
	mouseActToolBar->addAction(mouseMoveAct);
	mouseActToolBar->addAction(mouseZoomAct);
}



void MainWindow::openModel()
{
	glViewer_1->loadMeshModel();
}

void MainWindow::openTextureImage()
{
	glViewer_1->loadTextureImage();
	glViewer_2->loadTextureImage();
}
void MainWindow::openQuadFile()
{
	glViewer_1->loadQuadFile();
}
void MainWindow::saveModel()
{
	glViewer_2->saveMeshModel();
}
void MainWindow::saveAsBmp(){}
void MainWindow::recentFiles(){}
void MainWindow::toobBar(){}
void MainWindow::stateBar(){}
void MainWindow::about(){}
MainWindow::~MainWindow(){}
