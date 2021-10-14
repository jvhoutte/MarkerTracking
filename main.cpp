#include "mainwindow.h"
#include <QApplication>

#include "QVTKOpenGLWidget.h"
#include <QSurfaceFormat>

#include "vtkInteractorStyleImage.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkJPEGReader.h"
#include "vtkImageViewer.h"

#include "vtkComputeDistanceField.h"
#include "GeometricROI.h"


int main(int argc, char *argv[])
{
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    
    return a.exec();
}
