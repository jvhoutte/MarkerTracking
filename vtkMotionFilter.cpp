#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"

#include "vtkMotionFilter.h"


vtkStandardNewMacro(vtkMotionFilter);

vtkMotionFilter::vtkMotionFilter()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  
  motiondata.resize(1); // 1 time frame
  calibrationdata.resize(2); // 2 compounds
}
 
vtkMotionFilter::~vtkMotionFilter()
{
}
 
int vtkMotionFilter::FillInputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    return 1;
    }
  return 0;
}
 
 
int vtkMotionFilter::FillOutputPortInformation(int vtkNotUsed(portNumber), vtkInformation *info)
{

  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");

  return 1;
}
 
 
void vtkMotionFilter::SetVirtualLandmarks(int modelidx, vtkSmartPointer<vtkPoints> lm_virt){calibrationdata[modelidx].landmarks_virtual->DeepCopy(lm_virt);
std::cout << "Saving " << calibrationdata[modelidx].landmarks_virtual->GetNumberOfPoints() << " virtual landmarks" << std::endl;
}
void vtkMotionFilter::SetRealLandmarks(int modelidx, vtkSmartPointer<vtkPolyData> lm_real){calibrationdata[modelidx].landmarks_real->DeepCopy(lm_real);}
void vtkMotionFilter::SetVirtualSurface(int modelidx, vtkSmartPointer<vtkPolyData> surf_virt){calibrationdata[modelidx].model_virtual->DeepCopy(surf_virt);}
void vtkMotionFilter::SetRealSurface(int modelidx, vtkSmartPointer<vtkPolyData> surf_real){calibrationdata[modelidx].model_real->DeepCopy(surf_real);}
 
int vtkMotionFilter::RequestData(vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the input and output
    this->surfaceRest = vtkMultiBlockDataSet::GetData(inputVector[0],0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    this->surfacePose = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    this->surfacePose->DeepCopy(surfaceRest);

    // Transform 
    
    vtkSmartPointer<vtkTransform> totaltransf = vtkSmartPointer<vtkTransform>::New();
    totaltransf->Identity();
    totaltransf->PreMultiply();
    totaltransf->Concatenate(calibrationdata[0].GetInverseTransform());
    totaltransf->Concatenate(motiondata[timeframe][0].GetInverseTransform());
    totaltransf->Concatenate(motiondata[timeframe][1].GetTransform());
    totaltransf->Concatenate(calibrationdata[1].GetTransform());
    
    // only transform the tibia 
    vtkSmartPointer<vtkTransformPolyDataFilter> transformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformer->SetTransform(totaltransf);
    transformer->SetInputData(vtkPolyData::SafeDownCast(surfaceRest->GetBlock(1)));
    transformer->Update();
    surfacePose->SetBlock(1, transformer->GetOutput());
    
    return 1;
}
 
bool vtkMotionFilter::Calibrate(int modelidx){
    return calibrationdata[modelidx].Calibrate(modelidx);
}
 
void vtkMotionFilter::AddTimeFrame(int time, int model, double q1, double q2, double q3, double q4, double t1, double t2, double t3){
    if(time>=motiondata.size()){ // Add time frame to the motiondata
        MultiCompoundTransform* multicompound = new MultiCompoundTransform();
        motiondata.push_back(*multicompound);
    }
    
    motiondata[time][model].SetParameters(q1,q2,q3,q4,t1,t2,t3);  
    
    
    
}
 
//----------------------------------------------------------------------------
void vtkMotionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}














////////////////////////////////////////////////////////////////////////////////////////////////////////



vtkStandardNewMacro(vtkGeneratePointCloudFilter);

vtkGeneratePointCloudFilter::vtkGeneratePointCloudFilter()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  framedata.resize(1);
}
 
vtkGeneratePointCloudFilter::~vtkGeneratePointCloudFilter()
{
}
 
int vtkGeneratePointCloudFilter::FillOutputPortInformation(int vtkNotUsed(portNumber), vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
 
int vtkGeneratePointCloudFilter::RequestData(vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the input and output
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    this->pointcloud = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    
    //pointcloud = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    
    // Transform
    for(int markeridx=0; markeridx<framedata.size(); markeridx++){
        // Transform of the bone 
        vtkSmartPointer<vtkTransform> transform1 = framedata[markeridx][1].GetTransform();
        transform1->Inverse();
        
        // Transform of the bone 
        vtkSmartPointer<vtkTransform> transform = framedata[markeridx][0].GetTransform();
        transform->PostMultiply();
        transform->Concatenate( transform1 );
        
        double point[3] = {0,0,0};
        double out[3] = {0,0,0};
        transform->TransformPoint(point, out);
        points->InsertNextPoint(out);
    }
    
    poly->SetPoints(points);
    
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(poly);
    vertexFilter->Update();

    pointcloud->DeepCopy(vertexFilter->GetOutput());
    
    return 1;
}
 
 
void vtkGeneratePointCloudFilter::AddTimeFrame(int time, int model, double q1, double q2, double q3, double q4, double t1, double t2, double t3){
    if(time>=framedata.size()){ // Add time frame to the framedata
        MultiCompoundTransform* multicompound = new MultiCompoundTransform();
        framedata.push_back(*multicompound);
    }
    
    framedata[time][model].SetParameters(q1,q2,q3,q4,t1,t2,t3);  
}
 
//----------------------------------------------------------------------------
void vtkGeneratePointCloudFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}



