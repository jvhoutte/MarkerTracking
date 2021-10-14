#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
    ImageSliceActor = vtkSmartPointer<vtkImageSlice>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    
    knee_api = new KneeAPI();
    
    //multiblock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    //motionfilter = vtkSmartPointer<vtkMotionFilter>::New();
    
    // Annotation style
    style_femur = vtkSmartPointer<AnnotationStyle>::New();
    style_tibia = vtkSmartPointer<AnnotationStyle>::New();
    
    
    ui->hSliderTIME->setEnabled(false);
    
    ui->btnSaveLandmarksFemur->setEnabled(false);
    ui->btnSaveLandmarksTibia->setEnabled(false);
    
    closestpoints = std::vector<std::vector<int>>(2);
    
    ///////////////////////// TEMPORARY CODE  /////////////////////////////////////////////////////////////
    // The functionalities used here are accesible through the Qt buttons                                //
    // This code should be removed in final version and is only meant for convenience for the developper //
    
    knee_api->openFile(0, "/home/jeroen/Documents/KNEE-project/InputData/STL/femur_remesh.stl");
    knee_api->openFile(1, "/home/jeroen/Documents/KNEE-project/InputData/STL/tibia_remesh.stl");
    
    knee_api->OpenMarkerFile("/home/jeroen/Documents/KNEE-project/InputData/OpticalData/markers/femur_summary.txt",0,"landmarks");
    knee_api->OpenMarkerFile("/home/jeroen/Documents/KNEE-project/InputData/OpticalData/markers/tibia_summary.txt",1,"landmarks");
    knee_api->OpenMarkerFile("/home/jeroen/Documents/KNEE-project/InputData/OpticalData/surfaceprobe/femur_surface.txt",0,"surfprobe");
    knee_api->OpenMarkerFile("/home/jeroen/Documents/KNEE-project/InputData/OpticalData/surfaceprobe/tibia_surface.txt",1,"surfprobe");
        
    knee_api->OpenMotionFile("/home/jeroen/Documents/KNEE-project/InputData/OpticalData/motion/motiondata.txt");
    
    knee_api->SetVirtualLandmarks(0, "virtlandmarks_femur.vtk");
    knee_api->SetVirtualLandmarks(1, "virtlandmarks_tibia.vtk");
    
    ShowVirtualModel(0);
    ShowVirtualModel(1);
    
    knee_api->Calibrate(0);
    knee_api->Calibrate(1);
    
}

MainWindow::~MainWindow()
{
    delete ui;
}

std::string MainWindow::GetSelectedFile(){
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open STL file"), "", tr("STL file (*.stl);;All Files (*)"));
    return fileName.toUtf8().constData();
}

void MainWindow::ShowVirtualModel(int blockid)
{    
    // Draw the model in the single viewer
    drawSingleSurfaceModel(blockid);

    // Draw the model in the motion display window
    if(knee_api->GetRealModel()->GetNumberOfBlocks()==2){
        drawRealSurfaceModel(0, false);
    }
    
    // Enable button to save landmarks 
    if(blockid==0){
        ui->btnSaveLandmarksFemur->setEnabled(true);
    }
    else if(blockid==1){
        ui->btnSaveLandmarksTibia->setEnabled(true);
    }
         
}


void MainWindow::drawSingleSurfaceModel(int modelid){
    vtkSmartPointer<vtkPolyDataMapper> polymapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    polymapper->SetInputDataObject(knee_api->GetVirtualModel(modelid));
    polymapper->ScalarVisibilityOn();
    polymapper->SetScalarRange(0,0);
    
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(polymapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(.3, .6, .3); // Background color green
    
    if(modelid == 0){
        style_femur->SetDefaultRenderer(renderer);
        style_femur->SetData( knee_api->GetVirtualModel(modelid) );
        
        ui->femurRenderer->GetRenderWindow()->AddRenderer(renderer);
        ui->femurRenderer->GetRenderWindow()->Render();
        ui->femurRenderer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style_femur);
    }
    else if(modelid == 1){
        style_tibia->SetDefaultRenderer(renderer);
        style_tibia->SetData( knee_api->GetVirtualModel(modelid) );
        
        ui->tibiaRenderer->GetRenderWindow()->AddRenderer(renderer);
        ui->tibiaRenderer->GetRenderWindow()->Render();
        ui->tibiaRenderer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style_tibia);
    }
}

vtkSmartPointer<vtkTable> CvrtDataToTable(std::vector<double> data){
    
        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
        table->AddColumn(vtkSmartPointer<vtkDoubleArray>::New());
        table->AddColumn(vtkSmartPointer<vtkDoubleArray>::New());
        table->GetColumn(0)->SetName("X");
        table->GetColumn(1)->SetName("Y");

        table->SetNumberOfRows(data.size());
        for(int i = 0; i < data.size(); ++i)
        {
            table->SetValue(i, 0, i);
            table->SetValue(i, 1, data[i]);
        }
        
        return table;
}

void MainWindow::ShowGraph(std::vector<std::vector<double>> data){
        
    // Prepare plot data
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    vtkSmartPointer<vtkPlot> line1 = chart->AddPlot(vtkChart::LINE);
    line1->SetInputData(CvrtDataToTable(data[0]), 0, 1);
    line1->SetColor(0,1,0);
    vtkSmartPointer<vtkPlot> line2 = chart->AddPlot(vtkChart::LINE);
    line2->SetInputData(CvrtDataToTable(data[1]), 0, 1);
    line2->SetColor(1,0,0);
    
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->SetRenderWindow(ui->vtkRenderChart->GetRenderWindow());
    view->GetScene()->AddItem(chart);
       
    ui->vtkRenderChart->SetRenderWindow( view->GetRenderWindow() );
}

void MainWindow::on_btnCalcDistances_clicked(){
    std::vector<std::vector<double>> distvalues = knee_api->CalculateDistances(&closestpoints);
    
    ShowGraph( distvalues );
}

vtkSmartPointer<vtkPolyData> MainWindow::GetMarker(double* pos){
    // Get the bounds of the object such that we know how large we can make the marker
    double bounds[6];
    knee_api->GetVirtualModel(1)->GetBounds(bounds);
    double diagonal = std::pow(std::pow(bounds[1]-bounds[0],2)+std::pow(bounds[3]-bounds[2],2)+std::pow(bounds[5]-bounds[4],2),0.5);

    // Create a sphere
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(pos[0], pos[1], pos[2]);
    sphereSource->SetRadius(diagonal/200.);
    sphereSource->Update();

    return sphereSource->GetOutput();
}

void MainWindow::drawClosestPoint(int time)
{
    // Keep camera direction 
    vtkCamera* camera = NULL;
    if(ui->vtkRenderTibia->GetRenderWindow()->GetRenderers()->GetFirstRenderer() != NULL){    
        camera = ui->vtkRenderTibia->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera();
    }
    
    // Get the virtual model 
    vtkSmartPointer<vtkPolyData> poly = knee_api->GetVirtualModel(1, true);
    
    // Get point coordinates
    double p1[3]; double p2[3];
    poly->GetPoint(closestpoints[0][time], p1);
    poly->GetPoint(closestpoints[1][time], p2);
    
    vtkSmartPointer<vtkPolyDataMapper> markermapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    markermapper1->SetInputData(GetMarker(p1));
    vtkSmartPointer<vtkPolyDataMapper> markermapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    markermapper2->SetInputData(GetMarker(p2));
    vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
    actor1->SetMapper(markermapper1);
    actor1->GetProperty()->SetColor(0,1,0);
    vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
    actor2->SetMapper(markermapper2);
    actor2->GetProperty()->SetColor(1,0,0);
    
    
    vtkSmartPointer<vtkDataArray> arr = vtkPolyData::SafeDownCast(knee_api->GetRealModel(time,true)->GetBlock(1))->GetPointData()->GetArray("SignedDistances");
    if(arr!=NULL){ poly->GetPointData()->SetScalars(arr); }
    
    vtkSmartPointer<vtkPolyDataMapper> polymapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    polymapper->SetInputDataObject(poly);
    polymapper->SetScalarVisibility(true);
    polymapper->SetScalarRange(0,10);
  
    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(polymapper->GetLookupTable());
    scalarBar->SetTitle("Distance");
    scalarBar->SetNumberOfLabels(2);
        
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(polymapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    ui->vtkRenderTibia->GetRenderWindow()->AddRenderer(renderer);
    
    if(camera != NULL){ renderer->SetActiveCamera(camera); }
    
    renderer->AddActor(actor);
    renderer->AddActor(actor1);
    renderer->AddActor(actor2);
    renderer->AddActor2D(scalarBar);
    renderer->SetBackground(.3, .6, .3); // Background color green
    
    ui->vtkRenderTibia->GetRenderWindow()->Render();
        
}


void MainWindow::drawRealSurfaceModel(int time, bool keepviewdirection)
{
    vtkCamera* camera; 
    if(keepviewdirection){
        camera = ui->vtkRenderModel->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera();
    }
    
    vtkSmartPointer<vtkCompositePolyDataMapper2> polymapper = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
    polymapper->SetInputDataObject(knee_api->GetRealModel(time, showROIonly));
    polymapper->SetScalarVisibility(true);
    polymapper->SetScalarRange(0,10);
  
    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(polymapper->GetLookupTable());
    scalarBar->SetTitle("Distance");
    scalarBar->SetNumberOfLabels(2);
    
    if(showROIonly==true){polymapper->SetBlockOpacity(0, 0);}
    
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(polymapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    if(keepviewdirection){
        renderer->SetActiveCamera(camera);
    }
    ui->vtkRenderModel->GetRenderWindow()->AddRenderer(renderer);

    renderer->AddActor(actor);
    renderer->AddActor2D(scalarBar);
    renderer->SetBackground(.3, .6, .3); // Background color green
    
    ui->vtkRenderModel->GetRenderWindow()->Render();
}

void MainWindow::on_btnShowROI_clicked()
{
    showROIonly = !showROIonly;
    drawRealSurfaceModel( ui->hSliderTIME->sliderPosition(), true);
}

void MainWindow::on_btnOpenFemur_clicked()
{
    std::string filename = GetSelectedFile();
    knee_api->openFile(0, filename);
    ShowVirtualModel(0);
}

void MainWindow::on_btnOpenTibia_clicked()
{
    std::string filename = GetSelectedFile();
    knee_api->openFile(1, filename);
    ShowVirtualModel(1);
}

void MainWindow::SaveLandmarks(vtkSmartPointer<vtkPoints> points, std::string filename){
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(filename.c_str());
    wr->SetInputData(poly);
    wr->Update();
}

void MainWindow::on_btnSaveLandmarksFemur_clicked()
{
    knee_api->SetVirtualLandmarks(0, style_femur->GetLandmarks() );    
    SaveLandmarks(style_femur->GetLandmarks(), "virtlandmarks_femur.vtk");
}

void MainWindow::on_btnSaveLandmarksTibia_clicked()
{
    knee_api->SetVirtualLandmarks(1, style_tibia->GetLandmarks() );    
    SaveLandmarks(style_tibia->GetLandmarks(), "virtlandmarks_tibia.vtk");
}

void MainWindow::on_btnOpenRealSurfaceFemur_clicked(){ knee_api->OpenMarkerFile(GetSelectedFile(),0,"surfprobe"); }
void MainWindow::on_btnOpenRealSurfaceTibia_clicked(){ knee_api->OpenMarkerFile(GetSelectedFile(),1,"surfprobe"); }
void MainWindow::on_btnOpenRealLandmarksFemur_clicked(){ knee_api->OpenMarkerFile(GetSelectedFile(),0,"landmarks"); }
void MainWindow::on_btnOpenRealLandmarksTibia_clicked(){ knee_api->OpenMarkerFile(GetSelectedFile(),1,"landmarks"); }

void MainWindow::on_btnCalibrateFemur_clicked(){ knee_api->Calibrate(0); }
void MainWindow::on_btnCalibrateTibia_clicked(){ knee_api->Calibrate(1); }


void MainWindow::on_btnOpenMotion_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open motion file"), "", tr("TXT file (*.txt);;All Files (*)"));
    std::string filename = fileName.toUtf8().constData();
    knee_api->OpenMotionFile(filename);
    ui->hSliderTIME->setEnabled(true);
    ui->hSliderTIME->setMinimum(0);
    ui->hSliderTIME->setMaximum(knee_api->GetNumberOfTimeFrames()-1);
}

void MainWindow::on_hSliderTIME_sliderMoved(int position)
{
    drawRealSurfaceModel(position, true);   
    if(closestpoints[0].size()>position && closestpoints[1].size()>position){
        drawClosestPoint(position);
    }
}
