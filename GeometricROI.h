#ifndef GEOMROI_H
#define GEOMROI_H


#include "Defs.h"

class GeometricROI{
public:
    
    GeometricROI(){}
    
    vtkSmartPointer<vtkMultiBlockDataSet> GetJointSurface( vtkSmartPointer<vtkMultiBlockDataSet> multiblock );
    
    vtkSmartPointer<vtkMultiBlockDataSet> SplitBySymmetryPlane( vtkSmartPointer<vtkMultiBlockDataSet> multiblock );
    
    
    
    
protected:
    
    void GetCOM(vtkSmartPointer<vtkPolyData> object, double (*center)[3]);
    void GetElongationAxis(vtkSmartPointer<vtkPolyData> object, double (*axis)[3]);
    
    vtkSmartPointer<vtkPolyData> CutModel(vtkSmartPointer<vtkPolyData> object, double (*point1)[3], double (*point2)[3]);
    void GetEndPoints(vtkSmartPointer<vtkPolyData> object, double (*p1)[3], double (*p2)[3]);
    
    vtkSmartPointer<vtkPolyData> SplitBySymmetryPlane( vtkSmartPointer<vtkPolyData> object );
    vtkQuaterniond RotationFromReferenceAxis(double referenceAxis[3], double axis[3]);
    void FitPlane(vtkSmartPointer<vtkPolyData> source, vtkSmartPointer<vtkPolyData> mirror, double (*origin)[3], double (*evec)[3]);
    vtkSmartPointer<vtkTransform> AlignSymmetryPlane(vtkSmartPointer<vtkPolyData> source);
private:
    
};


#endif
