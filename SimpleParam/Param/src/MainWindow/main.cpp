#include <vld.h>
#include "MainWindow.h"
#include <QMainWindow>
#include <QtOpenGL>

#include "arthurstyle.h"

int main(int argc, char *argv[])
{

	QApplication::setColorSpec( QApplication::CustomColor );
	QApplication app(argc,argv);

	if ( !QGLFormat::hasOpenGL() ) {
		QString msg = "System has no OpenGL support!";
		QMessageBox::critical( 0, QString("OpenGL"), msg + QString(argv[1]) );
		return -1;
	}
	// create widget
	QStyle *arthurStyle = new ArthurStyle();
	MainWindow mainWin;
	mainWin.setStyle(arthurStyle);

	mainWin.show();

	return app.exec();
}
