#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageViewer2.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSphereSource.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkImageActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkSTLReader.h>

#include <sstream>
#include <iostream>
#include <string>

#include "Defs.h"
#include "vtkAnnotateMarker.h"
#include "vtkMotionFilter.h"
#include "vtkComputeDistanceField.h"
#include "GeometricROI.h"

#include "KneeAPI.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:
    
    void drawRealSurfaceModel(int time=0, bool keepviewdirection = true);

    
    void on_btnOpenFemur_clicked();
    void on_btnOpenTibia_clicked();
    void on_btnSaveLandmarksFemur_clicked();    
    void on_btnSaveLandmarksTibia_clicked();
    
    void on_btnOpenRealSurfaceFemur_clicked();
    void on_btnOpenRealSurfaceTibia_clicked();
    void on_btnOpenRealLandmarksFemur_clicked();
    void on_btnOpenRealLandmarksTibia_clicked();
    
    void on_btnCalibrateFemur_clicked();
    void on_btnCalibrateTibia_clicked();
    
    void on_btnShowROI_clicked();
    void on_btnOpenMotion_clicked();
    
    void on_btnCalcDistances_clicked();
    void on_hSliderTIME_sliderMoved(int position);
    
    void ShowGraph(std::vector<std::vector<double>> data);
    
    void drawSingleSurfaceModel(int modelid);
    
    std::string GetSelectedFile();
    
    void SaveLandmarks(vtkSmartPointer<vtkPoints> points, std::string filename);
    
    void ShowVirtualModel(int blockid);
    
    void drawClosestPoint(int time);
    
    vtkSmartPointer<vtkPolyData> GetMarker(double* pos);
    
private:
    Ui::MainWindow *ui;
    vtkSmartPointer<vtkImageSliceMapper> ImageSliceMapper;
    vtkSmartPointer<vtkImageSlice> ImageSliceActor;
    vtkSmartPointer<vtkRenderer> renderer;
    
    KneeAPI* knee_api;
    
    int mMinSliderX;
    int mMaxSliderX;
    
    vtkSmartPointer<AnnotationStyle> style_femur;
    vtkSmartPointer<AnnotationStyle> style_tibia;
    
    bool showROIonly = false;
    
    std::vector<std::vector<int>> closestpoints;
    
};

#endif // MAINWINDOW_H
