#ifndef __vtkDistanceField_h
#define __vtkDistanceField_h
 
#include <numeric>
#include "Defs.h"
 
class vtkComputeDistanceField : public vtkDataObjectAlgorithm 
{
public:
  vtkTypeMacro(vtkComputeDistanceField,vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkComputeDistanceField *New();
  
  void SetRefModel(vtkSmartPointer<vtkMultiBlockDataSet> multiblock, int refid);
  
protected:
  vtkComputeDistanceField();
  ~vtkComputeDistanceField();
 
  int FillInputPortInformation( int port, vtkInformation* info ) VTK_OVERRIDE;
  int FillOutputPortInformation(int portNumber, vtkInformation *info) VTK_OVERRIDE;
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  vtkSmartPointer<vtkPolyData> EvaluateDistanceFunction(vtkSmartPointer<vtkPolyData> object);
  
private:
  vtkComputeDistanceField(const vtkComputeDistanceField&);  // Not implemented.
  void operator=(const vtkComputeDistanceField&);  // Not implemented.
    
  vtkSmartPointer<vtkImplicitPolyDataDistance> distfunc;
  
 
};

#endif
