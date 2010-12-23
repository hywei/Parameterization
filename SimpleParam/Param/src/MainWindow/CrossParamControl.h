#ifndef CROSSPARAMCONTROL_H_
#define CROSSPARAMCONTROL_H_

#include "QGLViewer.h"

namespace PARAM
{
	class CrossParameter;
}

class CrossParamControl : public QWidget
{
	Q_OBJECT

public:
	CrossParamControl(QGLViewer* _gl_viewer_1, QGLViewer* _gl_viewer_2, QWidget* parent = 0);

private:
	QGroupBox* CreateSurface1Group(QWidget* parent = 0);
	QGroupBox* CreateSurface2Group(QWidget* parent = 0);
	QGroupBox* CreateCrossParamGroup(QWidget* parent = 0);
	QGroupBox* CreateTextureGroup(QWidget* parent = 0);
	QGroupBox* CreateVisualizationGroup(QWidget* parent = 0);
	QGroupBox* CreateOptimizerGroup(QWidget* parent = 0);
	QGroupBox* CreateCorrespondingGroup(QWidget* parent = 0);

	void CreateMainLayout();

private slots:
	void ComputeCrossParam();
	void OptimizeCrossParam();

private:
	QGLViewer* m_gl_viewer_1;
	QGLViewer* m_gl_viewer_2;

	QGroupBox* m_surface_1_group;
	QGroupBox* m_surface_2_group;
	QGroupBox* m_cross_param_group;
	QGroupBox* m_texture_setting_group;
	QGroupBox* m_visualization_group;
	QGroupBox* m_optimizer_group;
	QGroupBox* m_corresponding_group;

	boost::shared_ptr<PARAM::CrossParameter> p_cross_parameter;

};

#endif