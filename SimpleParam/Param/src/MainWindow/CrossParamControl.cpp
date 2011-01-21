#include "CrossParamControl.h"
#include "../Param/CrossParameter.h"
#include "../Param/Parameter.h"
#include "../Param/ParamDrawer.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>

CrossParamControl::CrossParamControl(QGLViewer* _gl_viewer_1, QGLViewer* _gl_viewer_2, QWidget* parent)
: QWidget(parent)
{
	m_gl_viewer_1 = _gl_viewer_1;
	m_gl_viewer_2 = _gl_viewer_2;    	
    
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
// 	QPushButton* optimize_ambiguity_patch = new QPushButton(surface_1_group);
// 	optimize_ambiguity_patch->setText(tr("Optimize Ambiguity Patch"));
	QPushButton* parameter_1 = new QPushButton(surface_1_group);
	parameter_1->setText(tr("Parameter"));

	/// layout
	QVBoxLayout* surface_1_layout = new QVBoxLayout(surface_1_group);
	surface_1_layout->addWidget(load_surface_1_mesh);
	surface_1_layout->addWidget(load_surface_1_patch);
//	surface_1_layout->addWidget(optimize_ambiguity_patch);
	surface_1_layout->addWidget(parameter_1);

	/// connections
	connect(load_surface_1_mesh, SIGNAL(clicked()), m_gl_viewer_1, SLOT(loadMeshModel()));
	connect(load_surface_1_patch, SIGNAL(clicked()), m_gl_viewer_1, SLOT(loadQuadFile()));
//	connect(optimize_ambiguity_patch, SIGNAL(clicked()), m_gl_viewer_1, SLOT(OptimizeAmbiguityPatch()));
	connect(parameter_1, SIGNAL(clicked()), m_gl_viewer_1, SLOT(SolveParameter()));

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
	QPushButton* parameter_2 = new QPushButton(surface_2_group);
	parameter_2->setText(tr("Parameter"));

	QPushButton* save_face_tex = new QPushButton(surface_2_group);
	save_face_tex->setText(tr("Save face tex coord"));
	QPushButton* load_face_tex = new QPushButton(surface_2_group);
	load_face_tex->setText(tr("Load face tex coord"));

	/// layout
	QVBoxLayout* surface_2_layout = new QVBoxLayout(surface_2_group);
	surface_2_layout->setSpacing(10);
	surface_2_layout->addWidget(load_surface_2_mesh);
	surface_2_layout->addWidget(load_surface_2_patch);
	surface_2_layout->addWidget(parameter_2);
	surface_2_layout->addWidget(save_face_tex);
	surface_2_layout->addWidget(load_face_tex);

	/// connections
	connect(load_surface_2_mesh, SIGNAL(clicked()), m_gl_viewer_2, SLOT(loadMeshModel()));
	connect(load_surface_2_patch, SIGNAL(clicked()), m_gl_viewer_2, SLOT(loadQuadFile()));
	connect(parameter_2, SIGNAL(clicked()), m_gl_viewer_2, SLOT(SolveParameter()));
	connect(save_face_tex, SIGNAL(clicked()), m_gl_viewer_2, SLOT(SaveFaceTexCoord()));
	connect(load_face_tex, SIGNAL(clicked()), m_gl_viewer_2, SLOT(LoadFaceTexCoord()));
	
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
	QCheckBox* outrange_vert = new QCheckBox(visual_group);
	QCheckBox* select_patch = new QCheckBox(visual_group);
	QCheckBox* flipped_triangle = new QCheckBox(visual_group);
	QCheckBox* corresponding = new QCheckBox(visual_group);
	QCheckBox* uncorresponding_vert = new QCheckBox(visual_group);


	patch_conner->setText(tr("Patch Conner"));
	patch_edge ->setText(tr("Patch Edge"));
	patch_face ->setText(tr("Patch Face"));
	outrange_vert->setText(tr("Out Range Vertex"));
	select_patch->setText(tr("Selected Patch"));
	flipped_triangle ->setText("Flipped Triangle");
	corresponding->setText("Corresponding");
	uncorresponding_vert->setText("UnCorresponding Vert");

	patch_conner->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	patch_edge->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	patch_face->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	outrange_vert->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	select_patch->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	flipped_triangle->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	corresponding->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	uncorresponding_vert->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	QVBoxLayout* visual_layout = new QVBoxLayout(visual_group);
	visual_layout->addWidget(patch_conner);
	visual_layout->addWidget(patch_edge);
	visual_layout->addWidget(patch_face);
	visual_layout->addWidget(outrange_vert);
	visual_layout->addWidget(select_patch);
	visual_layout->addWidget(flipped_triangle);
	visual_layout->addWidget(corresponding);
	visual_layout->addWidget(uncorresponding_vert);

	/// connects
	connect(patch_conner, SIGNAL(toggled(bool)), this, SLOT(SetPatchConnerDisplay(bool)));
    connect(patch_edge, SIGNAL(toggled(bool)), this, SLOT(SetPatchEdgeDisplay(bool)));
    connect(patch_face, SIGNAL(toggled(bool)), this, SLOT(SetPatchFaceDisplay(bool)));
	connect(outrange_vert, SIGNAL(toggled(bool)), this, SLOT(SetOutRangeVertDisplay(bool)));
	connect(select_patch, SIGNAL(toggled(bool)), this, SLOT(SetSelectedPatchDisplay(bool)));
	connect(flipped_triangle, SIGNAL(toggled(bool)), this, SLOT(SetFlippedTriangleDisplay(bool)));
	connect(corresponding, SIGNAL(toggled(bool)), this, SLOT(SetCorrespondingDisplay(bool)));
	connect(uncorresponding_vert, SIGNAL(toggled(bool)), this, SLOT(SetUnCorrespondingDisplay(bool)));

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

	QPushButton* load_corresponding = new QPushButton(tr("Load Corresponding"));
	QPushButton* corresponding = new QPushButton(tr("Corresponding"));
	QPushButton* select_patch = new QPushButton(tr("Select Patch"));

	connect(load_corresponding, SIGNAL(clicked()), this, SLOT(LoadCorrespondingFile()));
	connect(select_patch, SIGNAL(clicked()), m_gl_viewer_1, SLOT(SetParamDrawerSelectPatchMode()));
	connect(corresponding, SIGNAL(clicked()), m_gl_viewer_1, SLOT(SetParamDrawerSelectVertMode()));
	connect(corresponding, SIGNAL(clicked()), m_gl_viewer_2, SLOT(SetParamDrawerCorrespondMode()));
    
	/// layout
	QVBoxLayout* corresponding_layout = new QVBoxLayout(corresponding_group);
	corresponding_layout->addWidget(load_corresponding);
	corresponding_layout->addWidget(corresponding);
	corresponding_layout->addWidget(select_patch);

	return corresponding_group;
}

void CrossParamControl::LoadCorrespondingFile()
{
	std::string prev_file_path;
	ifstream fin ("open_file_path.txt");
	if(!fin.fail()){
		fin >> prev_file_path; 
	}
	fin.close();

	QString fileName = QFileDialog::getOpenFileName(
		this,
		tr("Open Corresponding File"),
		tr(prev_file_path.c_str()),
		tr("Corr file(*.corr);;"
		"All files(*)"));
	std::string f = std::string((const char *) fileName.toLocal8Bit());
	if(f.size() != 0)
	{
		p_cross_parameter = p_cross_parameter = boost::shared_ptr<PARAM::CrossParameter> (new
			PARAM::CrossParameter(*(m_gl_viewer_1->p_param), *(m_gl_viewer_2->p_param)));

		p_cross_parameter->LoadCorrespondingFile(f);
		
		std::string file_path, file_title, file_ext;
		Utility util;
		util.ResolveFileName(f, file_path, file_title, file_ext);
		ofstream fout("open_file_path.txt");
		if(!fout.fail()){
			fout << file_path << std::endl;
		}
		fout.close();
		
	} 
	
}

void CrossParamControl::FindCorrespondingOnA()
{
	int vid = m_gl_viewer_2->p_param_drawer->GetSelectedVertID();
	if(vid == -1) return;
	if(m_gl_viewer_2->p_param == NULL) return;
	if(m_gl_viewer_1->p_param_drawer == NULL) return;
	if(p_cross_parameter == NULL) return;
	int chart_id = m_gl_viewer_2->p_param->GetVertexChartID(vid);
	PARAM::ParamCoord param_coord = m_gl_viewer_2->p_param->GetVertexParamCoord(vid);
	PARAM::ChartParamCoord chart_param_coord(param_coord, chart_id);
	PARAM::SurfaceCoord surface_coord;
	p_cross_parameter->GetSurfaceCoordOnA(chart_param_coord, surface_coord);
	m_gl_viewer_1->p_param_drawer->SetSelectedVertCoord(surface_coord);
	m_gl_viewer_1->updateGL();
}

void CrossParamControl::FindCorrespondingOnB()
{
//	int vid = m_gl_viewer_1->p_param_drawer->GetSelectedVertID();
	PARAM::SurfaceCoord surf_coord = m_gl_viewer_1->p_param_drawer->GetSelectedSurfaceCorod();
	if(surf_coord.face_index == -1) return;
	if(m_gl_viewer_1->p_param == NULL) return;
	if(m_gl_viewer_2->p_param_drawer == NULL) return;
	if(p_cross_parameter == NULL) return;

	int fid = surf_coord.face_index;
	PARAM::Barycentrc baryc = surf_coord.barycentric;

	const PolyIndexArray& fIndex1 = m_gl_viewer_1->p_mesh->m_Kernel.GetFaceInfo().GetIndex();
	const PolyIndexArray& fIndex2 = m_gl_viewer_2->p_mesh->m_Kernel.GetFaceInfo().GetIndex();
	const CoordArray& vCoord2 = m_gl_viewer_2->p_mesh->m_Kernel.GetVertexInfo().GetCoord();

	const IndexArray& face1 = fIndex1[fid];
	std::vector<PARAM::SurfaceCoord> face_vert_sf_vec(3);
	std::vector<Coord> face_vert_coord_vec(3); 
	for(int i=0; i<3; ++i)
	{
		int vid = face1[i];   
		PARAM::SurfaceCoord surf_coord = p_cross_parameter->m_corresponding_AB[vid];
		int _fid = surf_coord.face_index;
		if(_fid == -1) { std::cerr <<" Fid == -1" << std::endl; return;}
		PARAM::Barycentrc _baryc = surf_coord.barycentric;
		const IndexArray& face2 = fIndex2[_fid];	
		face_vert_coord_vec[i] = vCoord2[face2[0]]*_baryc[0] + vCoord2[face2[1]]*_baryc[1] + vCoord2[face2[2]]*_baryc[2];

	//	int chart_id = m_gl_viewer_1->p_param->GetVertexChartID(vid);
	//	PARAM::ParamCoord param_coord = m_gl_viewer_1->p_param->GetVertexParamCoord(vid);
		//std::cout << "Surface A --- chart id : " << chart_id << ", param_coord " << param_coord.s_coord <<" " << param_coord.t_coord << std::endl; 
	//	PARAM::ChartParamCoord chart_param_coord(param_coord, chart_id);

		 face_vert_sf_vec[i]= p_cross_parameter->m_corresponding_AB[vid];
	}
	Coord corr_coord = face_vert_coord_vec[0]*baryc[0] + face_vert_coord_vec[1]*baryc[1] + face_vert_coord_vec[2]*baryc[2];

//	p_cross_parameter->GetSurfaceCoordOnB(chart_param_coord, surface_coord);
//	m_gl_viewer_2->p_param_drawer->SetSelectedVertCoord(surface_coord);
	m_gl_viewer_2->p_param_drawer->SetSelectedVertCoord(corr_coord);

	int vid_2 = m_gl_viewer_2->p_param_drawer->GetSelectedVertID();
//	std::cout << "Surface B --- chart id : " << chart_id << ", param_coord " << param_coord.s_coord <<" "<<param_coord.t_coord << std::endl;
	m_gl_viewer_2->updateGL();
}

void CrossParamControl::ComputeCrossParam()
{
	if(m_gl_viewer_1 == 0 || m_gl_viewer_2 == 0) return;
	if(m_gl_viewer_1->p_param == 0 || m_gl_viewer_2->p_param == 0) return;

	p_cross_parameter->FindCorrespondingAB();

// 	p_cross_parameter->FindCorrespondingBA();
// 	p_cross_parameter->VertTextureTransferBA();
// 
// 	const std::vector<TexCoord>& vert_tex_array_B = p_cross_parameter->GetTransferedVertTexArrayB();
// 
// 	boost::shared_ptr<MeshModel> p_mesh_B = m_gl_viewer_2->p_mesh;
// 	p_mesh_B->m_Kernel.GetVertexInfo().GetTexCoord() = vert_tex_array_B;
// 	p_mesh_B->m_Kernel.GetFaceInfo().GetTexIndex() = (p_mesh_B->m_Kernel.GetFaceInfo().GetIndex());

// 	p_cross_parameter->FaceTextureTransferBA();
// 
// 	const std::vector<TexCoordArray>& face_tex_array_B = p_cross_parameter->GetTransferedFaceTexArrayB();
// 
// 	boost::shared_ptr<MeshModel> p_mesh_B = m_gl_viewer_2->p_mesh;
// 	p_mesh_B->m_Kernel.GetFaceInfo().GetTexCoord() = face_tex_array_B;
	
	m_gl_viewer_2->p_param_drawer->SetUnCorrespondingVertArray(p_cross_parameter->GetUnCorrespondingVertArrayOnB());
}

void CrossParamControl::OptimizeCrossParam()
{
//	if(p_cross_param == NULL) return;
// 	boost::shared_ptr<PARAM::CrossParamHarmonicOptimizer> p_cp_optimizer(new
// 		PARAM::CrossParamHarmonicOptimizer(*p_cross_param.get()));
//
// 	p_cp_optimizer->Optimize();
}

void CrossParamControl::SetPatchConnerDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer) {
		m_gl_viewer_1->p_param_drawer->SetDrawPatchConner(toggled);
		m_gl_viewer_1->updateGL();
	}
	if(m_gl_viewer_2 && m_gl_viewer_2->p_param_drawer) {
		m_gl_viewer_2->p_param_drawer->SetDrawPatchConner(toggled);
		m_gl_viewer_2->updateGL();
	}
}

void CrossParamControl::SetPatchEdgeDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawPatchEdge(toggled);
		m_gl_viewer_1->updateGL();
	}
	if(m_gl_viewer_2 && m_gl_viewer_2->p_param_drawer){
		m_gl_viewer_2->p_param_drawer->SetDrawPatchEdge(toggled);
		m_gl_viewer_2->updateGL();
	}
}

void CrossParamControl::SetPatchFaceDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawPatchFace(toggled);
		m_gl_viewer_1->updateGL();
	}
	if(m_gl_viewer_2 && m_gl_viewer_2->p_param_drawer){
		m_gl_viewer_2->p_param_drawer->SetDrawPatchFace(toggled);
		m_gl_viewer_2->updateGL();
	}
}

void CrossParamControl::SetOutRangeVertDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawOutRangeVertices(toggled);
		m_gl_viewer_1->updateGL();
	}

	if(m_gl_viewer_2 && m_gl_viewer_2->p_param_drawer){
		m_gl_viewer_2->p_param_drawer->SetDrawOutRangeVertices(toggled);
		m_gl_viewer_2->updateGL();
	}
}

void CrossParamControl::SetSelectedPatchDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawSelectedPatch(toggled);
		m_gl_viewer_1->updateGL();
	}
}

void CrossParamControl::SetFlippedTriangleDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawFlipFace(toggled);
		m_gl_viewer_1->updateGL();
	}
}

void CrossParamControl::SetCorrespondingDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawSelectedVertex(toggled);
		m_gl_viewer_1->updateGL();
	}

	if(m_gl_viewer_2 && m_gl_viewer_2->p_param_drawer){
		m_gl_viewer_2->p_param_drawer->SetDrawSelectedVertex(toggled);
		m_gl_viewer_2->updateGL();
	}

}

void CrossParamControl::SetUnCorrespondingDisplay(bool toggled)
{
	if(m_gl_viewer_1 && m_gl_viewer_1->p_param_drawer){
		m_gl_viewer_1->p_param_drawer->SetDrawUnCorresponding(toggled);
		m_gl_viewer_1->updateGL();
	}

	if(m_gl_viewer_2 && m_gl_viewer_2->p_param_drawer){
		m_gl_viewer_2->p_param_drawer->SetDrawUnCorresponding(toggled);
		m_gl_viewer_2->updateGL();
	}
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
   main_group_layout->addWidget(m_corresponding_group);
   main_group_layout->addWidget(m_cross_param_group);
   main_group_layout->addWidget(m_texture_setting_group);
   main_group_layout->addWidget(m_visualization_group);
   main_group_layout->addStretch(1);

}
