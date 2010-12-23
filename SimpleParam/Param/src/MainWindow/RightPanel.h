#ifndef RIGHTPANEL_H
#define RIGHTPANEL_H

#include <QWidget>

class QTabWidget;
class QGLViewer;
class MainWindow;

class RightPanel : public QWidget
{
	Q_OBJECT

public:
	RightPanel(MainWindow* _main_window, QWidget* parent = 0);

public:
	QTabWidget* tabWidget;
};

class SimpleCrossParamTab : public QWidget
{
	Q_OBJECT

public:
	SimpleCrossParamTab(MainWindow* _main_window, QWidget* parent=0);

private slots:
	int LoadMesh1();
	int LoadMesh2();
	int LoadQuadPatch1();
	int LoadQuadPatch2();

	void SelectTexture();

	void ComputeCrossParam();

private:
	MainWindow* main_window;

};

#endif // RIGHTPANEL_H
