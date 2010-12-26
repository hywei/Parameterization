#include "RightPanel.h"
#include "QGLViewer.h"
#include <QLabel>
#include <QLineEdit>
#include <QLayout>
#include <QGroupBox>
#include <cstdio>
#include <QTabWidget>
#include <QPushButton>
#include <QValidator>
#include <QCheckBox>
#include <sstream>
#include "QGLViewer.h"
#include "MainWindow.h"

RightPanel::RightPanel(MainWindow* _main_window, QWidget* parent) : QWidget(parent)
{
	tabWidget = new QTabWidget;

	tabWidget->addTab(new SimpleCrossParamTab(_main_window, this), tr("Cross Parameterization"));

	QVBoxLayout* layout = new QVBoxLayout;
	layout->addWidget(tabWidget);

	setLayout(layout);
}

SimpleCrossParamTab::SimpleCrossParamTab(MainWindow* _main_window, QWidget* parent/* =0 */)
{
	QPushButton* load_mesh_1 = new QPushButton( tr("Load the first mesh"));
	connect(load_mesh_1, SIGNAL(clicked()), _main_window->glViewer_1, SLOT(loadMeshModel()));

	QPushButton* load_mesh_2 = new QPushButton( tr("Load the second mesh"));
	connect(load_mesh_2, SIGNAL(clicked()),_main_window->glViewer_2, SLOT(loadMeshModel()));

	QPushButton* load_quad_patch_1 = new QPushButton( tr("Load the first quad patch"));
	connect(load_quad_patch_1, SIGNAL(clicked()), _main_window->glViewer_1, SLOT(loadQuadFile()));

	QPushButton* load_quad_patch_2 = new QPushButton( tr("Load the second quad patch"));
	connect(load_quad_patch_2, SIGNAL(clicked()), _main_window->glViewer_2, SLOT(loadQuadFile()));

	QPushButton* compute_cross_param = new QPushButton( tr("Compute Cross Parameterization"));
	connect(compute_cross_param, SIGNAL(clicked()), this , SLOT(ComputeCrossParam()));



	QGridLayout* settingLayout = new QGridLayout;
	settingLayout->addWidget(load_mesh_1, 0, 0, 1, 1);
	settingLayout->addWidget(load_mesh_2, 1, 0, 1, 1);
	settingLayout->addWidget(load_quad_patch_1, 2, 0, 1, 1);
	settingLayout->addWidget(load_quad_patch_2, 3, 0, 1, 1);
	settingLayout->addWidget(compute_cross_param, 4, 0, 1, 1);

	QVBoxLayout* mainLayout = new QVBoxLayout;
	mainLayout->setSpacing(30);
	mainLayout->addLayout(settingLayout);
	mainLayout->addStretch();

	this->setLayout(mainLayout);

	main_window = _main_window;
}

int SimpleCrossParamTab::LoadMesh1()
{
	return 0;
}

int SimpleCrossParamTab::LoadMesh2()
{
	return 0;
}

int SimpleCrossParamTab::LoadQuadPatch1()
{
	return 0;
}

int SimpleCrossParamTab::LoadQuadPatch2()
{
	return 0;
}

void SimpleCrossParamTab::SelectTexture()
{

}

void SimpleCrossParamTab::ComputeCrossParam()
{
	if(main_window->glViewer_1->p_param == NULL ||
		main_window->glViewer_2->p_param == NULL) return;

 	// PARAM::CrossParam cross_param(*(main_window->glViewer_1->p_param),
 	// 	*(main_window->glViewer_2->p_param));
 
 	// cross_param.ComputeUintedDistortion();
	// cross_param.ComputeTexCoordOnSurface2();

	// const vector<int>& un_corresponding_vtx_array = cross_param.GetUnCorrespondingVtxArray();
	// main_window->glViewer_2->p_param->SetUnCorrespondingVtxArray(un_corresponding_vtx_array);

	// main_window->glViewer_2->p_param->FaceVaule2VtxColor(cross_param.GetUnitedDistortion());

	// main_window->glViewer_2->p_param->SetMeshTexCoord(cross_param.GetTexCoordOnSurface2());
};
