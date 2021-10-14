#include "vtkComputeDistanceField.h"



vtkStandardNewMacro(vtkComputeDistanceField);

vtkComputeDistanceField::vtkComputeDistanceField()
{
  this->SetNumberOfInputPorts(1);
  //this->SetNumberOfOutputPorts(1);
  
}
 
vtkComputeDistanceField::~vtkComputeDistanceField()
{
}
 
int vtkComputeDistanceField::FillInputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    return 1;
    }
  return 0;
}
 
 
int vtkComputeDistanceField::FillOutputPortInformation(int vtkNotUsed(portNumber), vtkInformation *info)
{

  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");

  return 1;
}
 
int vtkComputeDistanceField::RequestData(vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the input and output
    
    vtkSmartPointer<vtkPolyData> inputsurf = vtkPolyData::GetData(inputVector[0],0);
    
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkPolyData> outputsurf  = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkPolyData> outs = EvaluateDistanceFunction(inputsurf);
    
    outputsurf->DeepCopy(outs);
    
    return 1;
}
 

void vtkComputeDistanceField::SetRefModel(vtkSmartPointer<vtkMultiBlockDataSet> multiblock, int refid){
    distfunc = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    
    // Create normals for the reference model
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(vtkPolyData::SafeDownCast(multiblock->GetBlock(refid)));
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOff();
    normalGenerator->FlipNormalsOn();
    normalGenerator->Update();
    
    // Set the referencemodel 
    distfunc->SetInput(normalGenerator->GetOutput());
    
} 

vtkSmartPointer<vtkPolyData> vtkComputeDistanceField::EvaluateDistanceFunction(vtkSmartPointer<vtkPolyData> object){
    
    // Add distances to each point
    vtkSmartPointer<vtkDoubleArray> signedDistances = vtkSmartPointer<vtkDoubleArray>::New();
    signedDistances->SetNumberOfTuples(object->GetNumberOfPoints());
    signedDistances->SetNumberOfComponents(1);
    signedDistances->SetName("SignedDistances");
    
    // Evaluate the signed distance function at all of the grid points
    for(int pointid=0; pointid<object->GetNumberOfPoints(); pointid++){
        double point[3];
        object->GetPoint(pointid, point);
        
        double dist = - distfunc->EvaluateFunction(point);
        
        signedDistances->SetComponent(pointid, 0, dist);
    }
        
    object->GetPointData()->SetScalars(signedDistances);
        
    return object;
}
 
 
//----------------------------------------------------------------------------
void vtkComputeDistanceField::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
