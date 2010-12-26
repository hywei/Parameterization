#include "CrossParamControl.h"
#include "../Param/CrossParameter.h"
#include "../Param/QuadParameter.h"
#include "../Param/ParamDrawer.h"
#include "../ModelMesh/MeshModel.h"

CrossParamControl::CrossParamControl(QGLViewer* _gl_viewer_1, QGLViewer* _gl_viewer_2, QWidget* parent)
: QWidget(parent)
{
	m_gl_viewer_1 = _gl_viewer_1;
	m_gl_viewer_2 = _gl_viewer_2;
    
    //	p_cross_parameter = boost::shared_ptr<PARAM::CrossParameter> (new
    //		PARAM::CrossParameter(*(m_gl_viewer_1->p_param), *(m_gl_viewer_2->p_param)));
    
	m_surface_1_group = CreateSurface1Group(this);
	m_surface_2_group = CreateSurface2Group(this);
	m_cross_param_group = CreateCrossParamGroup(this);
	m_texture_setting_group = CreateTextureGroup(this);
	m_visualization_group = CreateVisualizationGroup(this);
    m_corresponding_group = CreateCorrespondingGroup(this);
    
	CreateMainLayout();
}

QGroupBox* CrossParamControl::CreateSurface1Group(QWidget* parent /* = 0 */)
{
	QGroupBox* surface_1_group = new QGroupBox(parent);
	surface_1_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	surface_1_group->setTitle(tr("Surface 1 Setting"));
    
	/// child widgets
	QPushButton* load_surface_1_mesh = new QPushButton(surface_1_group);
	load_surface_1_mesh->setText(tr("Load Surface 1 Mesh"));
	QPushButton* load_surface_1_patch = new QPushButton(surface_1_group);
	load_surface_1_patch->setText(tr("Load Surface 1 Patch"));

	/// layout
	QVBoxLayout* surface_1_layout = new QVBoxLayout(surface_1_group);
	surface_1_layout->addWidget(load_surface_1_mesh);
	surface_1_layout->addWidget(load_surface_1_patch);

	/// connections
	connect(load_surface_1_mesh, SIGNAL(clicked()), m_gl_viewer_1, SLOT(loadMeshModel()));
	connect(load_surface_1_patch, SIGNAL(clicked()), m_gl_viewer_1, SLOT(loadQuadFile()));

	return surface_1_group;
}

QGroupBox* CrossParamControl::CreateSurface2Group(QWidget* parent /* = 0 */)
{
	QGroupBox* surface_2_group = new QGroupBox(parent);
	surface_2_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	surface_2_group->setTitle(tr("Surface 2 Setting"));

	/// child widgets
	QPushButton* load_surface_2_mesh = new QPushButton(surface_2_group);
	load_surface_2_mesh->setText(tr("Load Surface 2 Mesh"));
	QPushButton* load_surface_2_patch = new QPushButton(surface_2_group);
	load_surface_2_patch->setText(tr("Load Surface 2 Patch"));

	/// layout
	QVBoxLayout* surface_2_layout = new QVBoxLayout(surface_2_group);
	surface_2_layout->setSpacing(10);
	surface_2_layout->addWidget(load_surface_2_mesh);
	surface_2_layout->addWidget(load_surface_2_patch);

	/// connections
	connect(load_surface_2_mesh, SIGNAL(clicked()), m_gl_viewer_2, SLOT(loadMeshModel()));
	connect(load_surface_2_patch, SIGNAL(clicked()), m_gl_viewer_2, SLOT(loadQuadFile()));
	
	return surface_2_group;
}

QGroupBox* CrossParamControl::CreateCrossParamGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* cross_param_group = new QGroupBox(parent);
	cross_param_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	cross_param_group->setTitle(tr("Cross Parameterization"));

	/// child widgets
	QPushButton* compute_cross_parameter = new QPushButton(cross_param_group);
	compute_cross_parameter->setText(tr("Compute Cross Parameter"));
	QPushButton* optimizer = new QPushButton(cross_param_group);
	optimizer->setText(tr("Cross Parameter Optimize"));

	/// layout
	QVBoxLayout* cross_parameter_layout = new QVBoxLayout(cross_param_group);
	cross_parameter_layout->addWidget(compute_cross_parameter);
	cross_parameter_layout->addWidget(optimizer);

	/// connections
	connect(compute_cross_parameter, SIGNAL(clicked()), this, SLOT(ComputeCrossParam()));
	connect(optimizer, SIGNAL(clicked()), this, SLOT(OptimizeCrossParam()));

	return cross_param_group;
}

QGroupBox* CrossParamControl::CreateTextureGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* texture_group = new QGroupBox(parent);
	texture_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	texture_group->setTitle(tr("Texture"));

	/// child widgets
	QRadioButton* square_texture = new QRadioButton(texture_group);
	QRadioButton* line_texture = new QRadioButton(texture_group);
	QRadioButton* boundary_texture = new QRadioButton(texture_group);

	square_texture->setText(tr("Square"));
	line_texture->setText(tr("Line"));
	boundary_texture->setText(tr("Boundary"));
	square_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	line_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	boundary_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	QVBoxLayout* texture_layout = new QVBoxLayout(texture_group);
	texture_layout->addWidget(square_texture);
	texture_layout->addWidget(line_texture);
	texture_layout->addWidget(boundary_texture);

	/// connections
	connect(square_texture, SIGNAL(clicked()), m_gl_viewer_1, SLOT(CreateSquareTexture()));
	connect(square_texture, SIGNAL(clicked()), m_gl_viewer_2, SLOT(CreateSquareTexture()));
	connect(line_texture, SIGNAL(clicked()), m_gl_viewer_1, SLOT(CreateLineTexture()));
	connect(line_texture, SIGNAL(clicked()), m_gl_viewer_2, SLOT(CreateLineTexture()));
	connect(boundary_texture, SIGNAL(clicked()), m_gl_viewer_1, SLOT(CreateBoundaryTexture()));
	connect(boundary_texture, SIGNAL(clicked()), m_gl_viewer_2, SLOT(CreateBoundaryTexture()));

	square_texture->setChecked(true);

	return texture_group;
}

QGroupBox* CrossParamControl::CreateVisualizationGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* visual_group = new QGroupBox(parent);
	visual_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	visual_group->setTitle(tr("Visulzation"));

	/// child widgets
	// QRadioButton* patch_layout = new QRadioButton(visual_group);
	QCheckBox* patch_conner = new QCheckBox(visual_group);
	QCheckBox* patch_edge = new QCheckBox(visual_group);
	QCheckBox* patch_face = new QCheckBox(visual_group);
	
// 	QRadioButton* param_texture = new QRadioButton(visual_group);
// 	QCheckBox* single_param_texture = new QCheckBox(param_texture);	
// 	QCheckBox* united_param_texture = new QCheckBox(param_texture);
// 
// 	QRadioButton* param_distortion = new QRadioButton(visual_group);
// 	QCheckBox* single_param_distortion = new QCheckBox(param_distortion);
// 	QCheckBox* united_param_distortion = new QCheckBox(param_distortion);

//	patch_layout->setText(tr("Patch Layout"));
	patch_conner->setText(tr("Patch Conner"));
	patch_edge ->setText(tr("Patch Edge"));
	patch_face ->setText(tr("Patch Face"));

// 	param_texture ->setText(tr("Parameter Texture"));
// 	single_param_texture->setText(tr("Single Parameter Texture"));
//    	united_param_texture->setText(tr("United Param Texture"));
// 
// 	param_distortion->setText(tr("Parameter Distortion"));
// 	single_param_distortion->setText(tr("Single Param Distortion"));
// 	united_param_distortion->setText(tr("United Param Distortion"));

// 	patch_layout->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);	
	patch_conner->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	patch_edge->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	patch_face->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
// 	single_param_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
// 	single_param_distortion->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
// 	united_param_texture->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
// 	united_param_distortion->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	QVBoxLayout* visual_layout = new QVBoxLayout(visual_group);
//	visual_layout->addWidget(patch_layout);
	visual_layout->addWidget(patch_conner);
	visual_layout->addWidget(patch_edge);
	visual_layout->addWidget(patch_face);
//	visual_layout->addWidget(param_texture);
//	visual_layout->addWidget(single_param_texture);
//	visual_layout->addWidget(united_param_texture);
//	visual_layout->addWidget(param_distortion);
//	visual_layout->addWidget(single_param_distortion);
//	visual_layout->addWidget(united_param_distortion);

	/// connects
	//connect(patch_conner, SIGNAL(clicked()), this, SLOT(DisplayPatchConner()));
      //connect(patch_edge, SIGNAL(clicked()), this, SLOT(DisplayPatchEdge()));
      //connect(patch_face, SIGNAL(clicked()), thsi, SLOT(DisplayPatchFace()));

	return visual_group;
}

QGroupBox* CrossParamControl::CreateOptimizerGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* optimizer_group = new QGroupBox(parent);
	optimizer_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	optimizer_group->setTitle("Cross Parameter Optimizer");

	QPushButton* optimizer = new QPushButton(optimizer_group);
	
	return optimizer_group;
}

QGroupBox* CrossParamControl::CreateCorrespondingGroup(QWidget* parent /* = 0 */)
{
	QGroupBox* corresponding_group = new QGroupBox(parent);
	corresponding_group->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
	corresponding_group->setTitle("Corresponding");

	QPushButton* corresponding = new QPushButton(tr("Corresponding"));
	connect(corresponding, SIGNAL(clicked()), m_gl_viewer_1, SLOT(SetParamDrawerSelectVertMode()));
	connect(corresponding, SIGNAL(clicked()), m_gl_viewer_2, SLOT(SetParamDrawerCorrespondMode()));
    
	/// layout
	QVBoxLayout* corresponding_layout = new QVBoxLayout(corresponding_group);
	corresponding_layout->addWidget(corresponding);
	return corresponding_group;
}

void CrossParamControl::ComputeCrossParam()
{
	if(m_gl_viewer_1 == 0 || m_gl_viewer_2 == 0) return;
	if(m_gl_viewer_1->p_param == 0 || m_gl_viewer_2->p_param == 0) return;

	// p_cross_parameter->FindCorrespondingAB();
	
	// m_gl_viewer_1->p_param_drawer->SetUnCorrespondingVertArray(p_cross_parameter->GetUnCorrespondingVertArrayOnA());
}

void CrossParamControl::OptimizeCrossParam()
{
//	if(p_cross_param == NULL) return;
// 	boost::shared_ptr<PARAM::CrossParamHarmonicOptimizer> p_cp_optimizer(new
// 		PARAM::CrossParamHarmonicOptimizer(*p_cross_param.get()));
//
// 	p_cp_optimizer->Optimize();
}

void CrossParamControl::CreateMainLayout()
{
   QGroupBox* main_group = new QGroupBox(this);
   main_group->setFixedWidth(200);
   main_group->setTitle(tr("Cross Parameterization"));

   QVBoxLayout* main_group_layout = new QVBoxLayout(main_group);
   main_group_layout->setMargin(3);
   main_group_layout->setSpacing(30);
   main_group_layout->addWidget(m_surface_1_group);
   main_group_layout->addWidget(m_surface_2_group);
   main_group_layout->addWidget(m_cross_param_group);
   main_group_layout->addWidget(m_texture_setting_group);
   main_group_layout->addWidget(m_visualization_group);
   main_group_layout->addWidget(m_corresponding_group);
   main_group_layout->addStretch(1);

}
