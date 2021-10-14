
#include "GeometricROI.h"





void GeometricROI::GetCOM(vtkSmartPointer<vtkPolyData> object, double (*center)[3]){
    vtkSmartPointer<vtkCenterOfMass> comfilter = vtkSmartPointer<vtkCenterOfMass>::New();
    comfilter->SetInputData(object);
    comfilter->SetUseScalarsAsWeights(false);
    comfilter->Update();

    comfilter->GetCenter(*center);
}


vtkSmartPointer<vtkMultiBlockDataSet> GeometricROI::GetJointSurface( vtkSmartPointer<vtkMultiBlockDataSet> multiblock ){
    
    // Get end points of model 1 
    double p11[3]; double p12[3];
    GetEndPoints(vtkPolyData::SafeDownCast(multiblock->GetBlock(0)), &p11, &p12);
    std::cout << p11[0] << " " << p11[1] << " " << p11[2] << std::endl;
    std::cout << p12[0] << " " << p12[1] << " " << p12[2] << std::endl; 
    
    // Get end points of model 2
    double p21[3]; double p22[3];
    GetEndPoints(vtkPolyData::SafeDownCast(multiblock->GetBlock(1)), &p21, &p22);
    
    std::cout << p21[0] << " " << p21[1] << " " << p21[2] << std::endl;
    std::cout << p22[0] << " " << p22[1] << " " << p22[2] << std::endl; 
    
    
    
    // Find set of closest points and use this to cut the polydatas at the correct side
    double d1 = vtkMath::Distance2BetweenPoints(p11,p21);
    double d2 = vtkMath::Distance2BetweenPoints(p11,p22);
    double d3 = vtkMath::Distance2BetweenPoints(p12,p21);
    double d4 = vtkMath::Distance2BetweenPoints(p12,p22);
    std::vector<double> v{d1, d2, d3, d4};
    std::vector<double>::iterator result = std::min_element(std::begin(v), std::end(v));
    int minelement = std::distance(std::begin(v), result);
   
    vtkSmartPointer<vtkMultiBlockDataSet> cutmodel = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    
    if(minelement==0 || minelement==1){cutmodel->SetBlock(0, CutModel(vtkPolyData::SafeDownCast(multiblock->GetBlock(0)), &p11, &p12));}
    else{cutmodel->SetBlock(0, CutModel(vtkPolyData::SafeDownCast(multiblock->GetBlock(0)), &p12, &p11));}
    
    if(minelement==0 || minelement==2){cutmodel->SetBlock(1, CutModel(vtkPolyData::SafeDownCast(multiblock->GetBlock(1)), &p21, &p22));}
    else{cutmodel->SetBlock(1, CutModel(vtkPolyData::SafeDownCast(multiblock->GetBlock(1)), &p22, &p21));}
    
    return cutmodel;
}

vtkSmartPointer<vtkPolyData> GeometricROI::CutModel(vtkSmartPointer<vtkPolyData> object, double (*point1)[3], double (*point2)[3]){
    double axis[3]; vtkMath::Subtract(*point1, *point2, axis);
    vtkMath::MultiplyScalar(axis,0.9);
    double planecenter[3];
    vtkMath::Add(*point2,axis,planecenter);
    vtkMath::Normalize(axis);
    
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(planecenter);
    plane->SetNormal(axis);
    
    vtkSmartPointer<vtkPlaneCollection> planes = vtkSmartPointer<vtkPlaneCollection>::New();
    planes->AddItem(plane);
    
    vtkSmartPointer<vtkClipClosedSurface> clipper = vtkSmartPointer<vtkClipClosedSurface>::New();
    clipper->SetClippingPlanes(planes);
    clipper->SetInputData(object);
    clipper->Update();
    
    return clipper->GetOutput();
}


void GeometricROI::GetEndPoints(vtkSmartPointer<vtkPolyData> object, double (*p1)[3], double (*p2)[3]){
    
    // Calculate COM 
    double com[3];
    GetCOM(object, &com);
    
    // Calculate elongation axis 
    double elongaxis[3];
    GetElongationAxis(object, &elongaxis);
    
    
    // Intersect the polydata with a line 
    vtkSmartPointer<vtkOBBTree> tree = vtkSmartPointer<vtkOBBTree>::New();
    tree->SetDataSet(object);
    tree->BuildLocator();

    double bounds[6];
    object->GetBounds(bounds);
    double length = pow( pow(bounds[1]-bounds[0],2) + pow(bounds[3]-bounds[2],2) + pow(bounds[5]-bounds[4],2), 0.5);
    
    // Intersect the locator with the line
	double lineP1[3] = { com[0] + length*elongaxis[0],  com[1] + length*elongaxis[1],  com[2] + length*elongaxis[2]};
	double lineP0[3] = { com[0] - length*elongaxis[0],  com[1] - length*elongaxis[1],  com[2] - length*elongaxis[2]};
	vtkSmartPointer<vtkPoints> intersectPoints = vtkSmartPointer<vtkPoints>::New();
    
	tree->IntersectWithLine(lineP0, lineP1, intersectPoints, NULL);
		
	// Get first and last intersection point
	intersectPoints->GetPoint(0, *p1);
    intersectPoints->GetPoint(intersectPoints->GetNumberOfPoints()-1, *p2);
}

void GeometricROI::GetElongationAxis(vtkSmartPointer<vtkPolyData> object, double (*axis)[3]){
    
    vtkSmartPointer<vtkDoubleArray> xArray = vtkSmartPointer<vtkDoubleArray>::New();
    xArray->SetNumberOfComponents(1);
    xArray->SetName("x");

    vtkSmartPointer<vtkDoubleArray> yArray = vtkSmartPointer<vtkDoubleArray>::New();
    yArray->SetNumberOfComponents(1);
    yArray->SetName("y");

    vtkSmartPointer<vtkDoubleArray> zArray = vtkSmartPointer<vtkDoubleArray>::New();
    zArray->SetNumberOfComponents(1);
    zArray->SetName("z");
    
    for(vtkIdType i = 0; i < object->GetNumberOfPoints(); i++){
        double p[3];
        object->GetPoint(i,p);
        xArray->InsertNextValue(p[0]);
        yArray->InsertNextValue(p[1]);
        zArray->InsertNextValue(p[2]);
    }
    
    vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
    datasetTable->AddColumn(xArray);
    datasetTable->AddColumn(yArray);
    datasetTable->AddColumn(zArray);  

    vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
    pcaStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );

    pcaStatistics->SetColumnStatus("x", 1 );
    pcaStatistics->SetColumnStatus("y", 1 );
    pcaStatistics->SetColumnStatus("z", 1 );
    
    pcaStatistics->RequestSelectedColumns();
    pcaStatistics->SetDeriveOption(true);
    pcaStatistics->Update();

    ///////// Eigenvectors ////////////
    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvectors(eigenvectors);

    vtkSmartPointer<vtkDoubleArray> evec1 = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvector(0, evec1);
    
    (*axis)[0]=evec1->GetValue(0);
    (*axis)[1]=evec1->GetValue(1);
    (*axis)[2]=evec1->GetValue(2);
    
    std::cout << evec1->GetValue(0) << " " << evec1->GetValue(1) << " " << evec1->GetValue(2) << std::endl;
}

vtkSmartPointer<vtkMultiBlockDataSet> GeometricROI::SplitBySymmetryPlane( vtkSmartPointer<vtkMultiBlockDataSet> object ){
    vtkSmartPointer<vtkMultiBlockDataSet> multiblock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    for(int blockid=0; blockid<object->GetNumberOfBlocks();blockid++){
        multiblock->SetBlock(blockid, SplitBySymmetryPlane(vtkPolyData::SafeDownCast(object->GetBlock(blockid))));
    }   
    return multiblock; 
}

vtkSmartPointer<vtkPolyData> GeometricROI::SplitBySymmetryPlane( vtkSmartPointer<vtkPolyData> object ){
    // Calculate COM 
    double com[3];
    GetCOM(object, &com);
    
    // Calculate elongation axis 
    double elongaxis[3];
    GetElongationAxis(object, &elongaxis);
    
    // Transform object such that the PC-axis aligns with world x-axis
    double refaxis[3] = {1,0,0};
    vtkQuaterniond quat = RotationFromReferenceAxis(elongaxis, refaxis);
	double rotaxis[3]; double rotangle = quat.GetRotationAngleAndAxis(rotaxis);
	vtkSmartPointer<vtkTransform> rotation = vtkSmartPointer<vtkTransform>::New();
	rotation->PostMultiply();
	rotation->Translate(-com[0], -com[1], -com[2]);
	rotation->RotateWXYZ(vtkMath::DegreesFromRadians(rotangle), rotaxis);
    
    vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformfilter->SetInputData(object);
    transformfilter->SetTransform(rotation);
    transformfilter->Update();
    
    vtkSmartPointer<vtkTransform> symmplanetransform = AlignSymmetryPlane(transformfilter->GetOutput());
    
    rotation->Concatenate( symmplanetransform );
    transformfilter->Update();
    
    vtkSmartPointer<vtkDoubleArray> symmarray = vtkSmartPointer<vtkDoubleArray>::New();
    symmarray->SetNumberOfTuples(object->GetNumberOfPoints());
    symmarray->SetNumberOfComponents(1);
    symmarray->SetName("SymmPlane");
    for(int i=0; i<object->GetNumberOfPoints(); i++){
        double p[3];
        transformfilter->GetOutput()->GetPoint(i,p);
        if(p[2]>0){symmarray->SetComponent(i,0,1);}
        else{symmarray->SetComponent(i,0,0);}
    }
    
    object->GetPointData()->AddArray(symmarray);
    
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(object);
    writer->SetFileName("/home/jeroen/Documents/KNEE-project/Output/symmobj.vtk");
    writer->Update();
    
    return object;
}


vtkQuaterniond GeometricROI::RotationFromReferenceAxis(double referenceAxis[3], double axis[3])
{
  vtkQuaterniond newOrientation;
  // Code greatly inspired by: http://www.fastgraph.com/makegames/3drotation/ .

  // Normalize. This is the unit vector in the "new Z" direction.
  const double epsilon = 1e-6;
  if (vtkMath::Norm(referenceAxis) < epsilon
    || vtkMath::Normalize(axis) < epsilon)
    {
    return newOrientation;
    }

  // The dot product of axis and the referenceAxis gives projection of
  // of axis on referenceAxis.
  double projection = vtkMath::Dot(axis, referenceAxis);

  // First try at making a View Up vector: use World Up.
  double viewUp[3];
  viewUp[0] = referenceAxis[0] - projection*axis[0];
  viewUp[1] = referenceAxis[1] - projection*axis[1];
  viewUp[2] = referenceAxis[2] - projection*axis[2];

  // Check for validity:
  double magnitude = vtkMath::Normalize(viewUp);
  if (magnitude < epsilon)
    {
    // Second try: Use Y axis default  (0,1,0).
    viewUp[0] = -axis[1]*axis[0];
    viewUp[1] = 1-axis[1]*axis[1];
    viewUp[2] = -axis[1]*axis[2];

    // Check for validity:
    magnitude = vtkMath::Normalize(viewUp);

    if (magnitude < epsilon)
      {
      // Final try: Use Z axis default  (0,0,1).
      viewUp[0] = -axis[2]*axis[0];
      viewUp[1] = -axis[2]*axis[1];
      viewUp[2] = 1-axis[2]*axis[2];

      // Check for validity:
      magnitude = vtkMath::Normalize(viewUp);

      if (magnitude < epsilon)
        {
        return newOrientation;
        }
      }
    }

  // Calculate the Right vector. Use cross product of axis and Up.
  double viewRight[3];
  vtkMath::Cross(viewUp, axis, viewRight);
  vtkMath::Normalize(viewRight); //Let's be paranoid about the normalization.

  // Get the rest transform matrix.
  newOrientation.SetRotationAngleAndAxis(acos(projection), viewRight);
  return newOrientation.Normalized();
}

void GeometricROI::FitPlane(vtkSmartPointer<vtkPolyData> source, vtkSmartPointer<vtkPolyData> mirror, double (*origin)[3], double (*evec)[3])
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	
	for(unsigned int i=0; i<source->GetNumberOfPoints(); i++){
		//Declare a variable to store the index of the point that gets added. This behaves just like an unsigned int.
		vtkIdType pid[1];
        
		//Add a point to the polydata and save its index, which we will use to create the vertex on that point.
		double p[3];
    	source->GetPoint(i,p);
    	
    	double q[3];
    	mirror->GetPoint(i,q);
    	
		pid[0] = points->InsertNextPoint((p[0]+q[0])/2, (p[1]+q[1])/2, (p[2]+q[2])/2);
		
		//create a vertex cell on the point that was just added.
		vertices->InsertNextCell ( 1,pid );
	}
	
	//create a polydata object
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	
	//set the points and vertices we created as the geometry and topology of the polydata
	polydata->SetPoints ( points );
	polydata->SetVerts ( vertices );


	/////////////////////////////
	// Fit plane to point set  //
	/////////////////////////////
	
	int npts = points->GetNumberOfPoints();
	// todo - error if no points or nPoints < 3
	double sumx=0.0, sumy=0.0, sumz=0.0;
	for (int i=0; i<npts; i++)
	{
		double *pt = points->GetPoint(i);
		sumx += pt[0];
		sumy += pt[1];
		sumz += pt[2];
	}
	double ox = sumx/npts; double oy = sumy/npts; double oz = sumz/npts;
	double sxx=0.0, syy=0.0, szz=0.0, sxy=0.0, sxz=0.0, syz=0.0;
	for (int i=0; i<npts; i++)
	{
		double* pt = points->GetPoint(i);
		double dx = ox-pt[0], dy = oy-pt[1], dz = oz-pt[2];
		sxx+=dx*dx; syy+=dy*dy; szz+=dz*dz;
		sxy+=dx*dy; sxz+=dx*dz; syz+=dy*dz;
	}
	double *mat[3],evals[3],*evecs[3];
	double mat0[3],mat1[3],mat2[3];
	double evec0[3],evec1[3],evec2[3];
	//setup
	mat[0] = mat0; mat[1] = mat1; mat[2] = mat2;
	evecs[0] = evec0; evecs[1] = evec1; evecs[2] = evec2;
	mat[0][0] = sxx; mat[0][1] = sxy; mat[0][2] = sxz;
	mat[1][0] = sxy; mat[1][1] = syy; mat[1][2] = syz;
	mat[2][0] = sxz; mat[2][1] = syz; mat[2][2] = szz;
	vtkMath::Jacobi(mat,evals,evecs);
	// vtkMath::Jacobi sorts results in dec order, vecs normalized
	// min evec is normal to plane - derive dip, dipdir
	// NB: evecs in columns
	(*evec)[0] = evecs[0][2]; (*evec)[1] = evecs[1][2]; (*evec)[2] = evecs[2][2];
	
    // want plane normal pointing UP
	/*
	if ((*evec)[2] < 0.0)
	{
		(*evec)[0] *= -1;
		(*evec)[1] *= -1;
		(*evec)[2] *= -1;
	}
	*/
    
    (*origin)[0] = ox;
    (*origin)[1] = oy;
    (*origin)[2] = oz;
}

vtkSmartPointer<vtkTransform> GeometricROI::AlignSymmetryPlane(vtkSmartPointer<vtkPolyData> source){

	////////////////////////
	// Mirror the mesh    //
	////////////////////////
	
	vtkSmartPointer<vtkTransform> mtransform = vtkSmartPointer<vtkTransform>::New();  
	mtransform->Scale(1,1,-1);

	std::cout << "mirrortransform... "  << std::endl;
	
	vtkSmartPointer<vtkTransformPolyDataFilter> mirrortransformfilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	mirrortransformfilter->SetInputData(source);
	mirrortransformfilter->SetTransform(mtransform);
	mirrortransformfilter->Update();
	
	std::cout << "reversefilter... "  << std::endl;
	
	vtkSmartPointer<vtkReverseSense> reverseSense = vtkSmartPointer<vtkReverseSense>::New();
	reverseSense->SetInputData(mirrortransformfilter->GetOutput());
	reverseSense->ReverseNormalsOn();
	reverseSense->Update();

	// scale+vtkreversesense or reflectionfilter	
	std::cout << "icp MATCHING... "  << std::endl;
	// Do ICP matching
	
	vtkSmartPointer<vtkIterativeClosestPointTransform> icp = vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
	icp->SetSource(source);
	icp->SetTarget(reverseSense->GetOutput());
	//	icp->DebugOn();
	icp->SetMaximumNumberOfIterations(100);
	icp->SetMaximumNumberOfLandmarks(source->GetNumberOfPoints());
	icp->SetCheckMeanDistance(1);
	icp->SetMaximumMeanDistance(0.0000001);
	icp->GetLandmarkTransform()->SetModeToRigidBody(); // rotation and translation only
	icp->Inverse();
	icp->Modified();
	icp->Update();
	
	// Transform the transformed mesh back by the ICP transformation solution
	vtkSmartPointer<vtkTransformPolyDataFilter> icpTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	icpTransformFilter->SetInputData(reverseSense->GetOutput());
	icpTransformFilter->SetTransform(icp); // Move 
	icpTransformFilter->Update();
	
	vtkSmartPointer<vtkPolyData> mirror = vtkSmartPointer<vtkPolyData>::New();
   	mirror->DeepCopy(icpTransformFilter->GetOutput());	
	
	// Calculate average distance between the two meshes
	vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter = vtkSmartPointer<vtkDistancePolyDataFilter>::New();
	distanceFilter->SetInputConnection( 0, icpTransformFilter->GetOutputPort() );
	distanceFilter->SetInputData( 1, source);
	distanceFilter->SignedDistanceOff(); 	
	distanceFilter->Update();
	
	
	double avgdist =0;
	int N = distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetNumberOfTuples();
	for(int i=0; i<N; i++){
		avgdist += *(distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetTuple(i));
	}		
	avgdist/=N;
	std::cout << "Avg dist " << avgdist << std::endl;
	
	// Check different starting conditions to avoid ending up in a local minimum during ICP.
	// The following code can be written shorter by using modified-function. Not sure if there is a way using this with the vtkIterativeClosestPointTransform object
	double dangle=90; // degrees
	double bestangle = 0;
	double minavgdist = avgdist;
	for(double it=dangle; it<360; it+=dangle){  
		mtransform->Identity();
		mtransform->RotateY(it); // At each iteration, an additional rotation around y axis (ie elongation axis) has been applied, before the mirroring by xy plane.
		mtransform->Scale(1,1,-1);
		mirrortransformfilter->SetTransform(mtransform);
		mirrortransformfilter->Update();
		reverseSense->Update();
		icp->Update();
		icpTransformFilter->Update();
		distanceFilter->Update();
		
		avgdist = 0;
		for(int i=0; i<N; i++){
			avgdist += *(distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetTuple(i));
		}		
		avgdist/=N;
		
		std::cout << it << "\t" << avgdist << std::endl; 
		if(avgdist < minavgdist){bestangle = it; minavgdist = avgdist; mirror->DeepCopy(icpTransformFilter->GetOutput());	}
	}

	std::cout << "Best starting rotation " << bestangle << std::endl;
	
	///////////////////////////////
	// Calculate symmetry plane  //
	///////////////////////////////
	
	double origin[3];
	double evec[3];
    
	FitPlane(source, mirror, &origin, &evec);

	double refaxis[3] = {0,0,1};
	vtkQuaterniond quat = RotationFromReferenceAxis(evec, refaxis);

	double rotaxis[3];
    double rotangle = quat.GetRotationAngleAndAxis(rotaxis);

	vtkSmartPointer<vtkTransform> rotation = vtkSmartPointer<vtkTransform>::New();
	rotation->RotateWXYZ(vtkMath::DegreesFromRadians(rotangle), rotaxis);
	vtkSmartPointer<vtkTransform> transl = vtkSmartPointer<vtkTransform>::New();
	transl->Translate(-origin[0],-origin[1],-origin[2]);
	transl->Concatenate(rotation);

	return transl;
}

