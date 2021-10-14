#ifndef VPICK_PUBLIC_H
#define VPICK_PUBLIC_H

#include <Defs.h>
#include <vtkFieldData.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkObjectFactory.h>
#include <vtkSphereSource.h>
#include <vtkDoubleArray.h>
#include <vtkSelectPolyData.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataGeometryFilter.h>
#include "vtkFloatArray.h"
//Code based on: https://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/CellPicking
 
// Define interaction style
class AnnotationStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static AnnotationStyle* New();
 
    AnnotationStyle();
 
    vtkSmartPointer<vtkPoints> GetLandmarks(){ return Landmark_points; }
    
    void SetData(vtkSmartPointer<vtkPolyData> data){this->Data->DeepCopy(data);}
   
  protected:
 
    void Init(vtkSmartPointer<vtkRenderer> renderer, vtkSmartPointer<vtkPolyDataMapper> mapper, vtkSmartPointer<vtkPolyData> object);
 
    vtkSmartPointer<vtkPolyData> GetMarker(double* pos);
 
    
	virtual void OnKeyPress();
	
	void UpdateVisualiser();
 
    
    private:
    vtkSmartPointer<vtkPoints> Landmark_points;
        
    vtkSmartPointer<vtkPolyData> Data;
    vtkSmartPointer<vtkPolyDataMapper> DataMapper;
    vtkSmartPointer<vtkDataSetMapper> selectedMapper;
    vtkSmartPointer<vtkActor> selectedActor;
    vtkSmartPointer<vtkDataSetMapper> landmarkMapper;
    vtkSmartPointer<vtkActor> landmarkActor;
    
    vtkSmartPointer<vtkMultiBlockDataSet> Landmarkmarkers;
    
    
};
 


#endif
