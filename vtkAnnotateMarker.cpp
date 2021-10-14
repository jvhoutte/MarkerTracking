#include "vtkAnnotateMarker.h"


vtkStandardNewMacro(AnnotationStyle);

AnnotationStyle::AnnotationStyle(){

    Data = vtkSmartPointer<vtkPolyData>::New();

    // vtkMultiBlockDataSets for marker visualisation
    Landmarkmarkers = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // Mapper and actor for new marker position
    landmarkMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    landmarkActor = vtkSmartPointer<vtkActor>::New();
    landmarkActor->SetMapper(landmarkMapper);
    landmarkActor->GetProperty()->SetColor(0,1,0);
    
    Landmark_points = vtkSmartPointer<vtkPoints>::New();

}


void AnnotationStyle::Init(vtkSmartPointer<vtkRenderer> renderer, vtkSmartPointer<vtkPolyDataMapper> mapper, vtkSmartPointer<vtkPolyData> object){
    this->SetDefaultRenderer(renderer);
    this->Data->DeepCopy(object);
    
    this->DataMapper = mapper;
    
    // Read arrays if available
    if(object->GetFieldData()->GetArray("Landmarks") != NULL){
      vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast( object->GetFieldData()->GetArray("Landmarks") );
      for(int i=0; i<array->GetNumberOfTuples();i++){
        double point[3]; array->GetTuple(i, point);
        Landmarkmarkers->SetBlock(Landmark_points->GetNumberOfPoints(), GetMarker(point));
        Landmark_points->InsertNextPoint(point);
      }
    }
    
    // Connect the actors of this class to the renderer
     UpdateVisualiser();
  }
  
  
  vtkSmartPointer<vtkPolyData> AnnotationStyle::GetMarker(double* pos){
        // Get the bounds of the object such that we know how large we can make the marker
        double bounds[6];
        Data->GetBounds(bounds);
        double diagonal = std::pow(std::pow(bounds[1]-bounds[0],2)+std::pow(bounds[3]-bounds[2],2)+std::pow(bounds[5]-bounds[4],2),0.5);

        // Create a sphere
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(pos[0], pos[1], pos[2]);
        sphereSource->SetRadius(diagonal/300.);
        sphereSource->Update();

        return sphereSource->GetOutput();
    }
    
    void AnnotationStyle::OnKeyPress()
	{
	  std::string key = this->GetInteractor()->GetKeySym();

	  if(key == "s"){ // add point to the loop
		// Get the location of the click (in window coordinates)
		int* pos = this->GetInteractor()->GetEventPosition();

		vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
		picker->SetTolerance(0.0005);
		picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

		if(picker->GetCellId() != -1){
		    Landmarkmarkers->SetBlock(Landmark_points->GetNumberOfPoints(), GetMarker(picker->GetPickPosition()));
		    Landmark_points->InsertNextPoint(picker->GetPickPosition());
		}
	  }
	  else if(key == "d"){
        if(Landmarkmarkers->GetNumberOfBlocks()>0){
            Landmarkmarkers->RemoveBlock(Landmarkmarkers->GetNumberOfBlocks()-1);
            
            vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
            for(int idx=0; idx<Landmark_points->GetNumberOfPoints()-1; idx++){
                newpoints->InsertNextPoint(Landmark_points->GetPoint(idx));
            }
            Landmark_points->DeepCopy(newpoints);
            
        }
	  }

	  UpdateVisualiser();
	}
  
  void AnnotationStyle::UpdateVisualiser(){
      
        // Update renderer with landmark markers
        vtkSmartPointer<vtkCompositeDataGeometryFilter> aggregatefilter2 = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
        aggregatefilter2->SetInputData(Landmarkmarkers);
        aggregatefilter2->Update();
        landmarkMapper->SetInputData(aggregatefilter2->GetOutput()); 
        landmarkActor->SetMapper(landmarkMapper);
        landmarkActor->GetProperty()->SetColor(1,0,0);
        this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(landmarkActor);
        this->Interactor->GetRenderWindow()->Render();	

        // Forward events
        vtkInteractorStyleTrackballCamera::OnKeyPress();
	}
 
