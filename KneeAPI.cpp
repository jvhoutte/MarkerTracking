#include "KneeAPI.h"


void split(const std::string &s, char delim, std::vector<std::string> &elems){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) { elems.push_back(item); }
}



KneeAPI::KneeAPI(){
    multiblock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    motionfilter = vtkSmartPointer<vtkMotionFilter>::New();
    
}


void KneeAPI::openFile(int blockid, std::string filename){
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    multiblock->SetBlock(blockid, reader->GetOutput());
    
    // Save the surface in the calibration class
    motionfilter->SetVirtualSurface(blockid, reader->GetOutput());
    
    // Determine joint region and mirror plane
    DetermineROI();
}

void KneeAPI::DetermineROI(){
    if(multiblock->GetNumberOfBlocks()==2){
        motionfilter->SetInputData(0, multiblock);
        motionfilter->Update();
        
        // Determine the ROI around the joint
        GeometricROI* roi = new GeometricROI();
        multiblock_roi = roi->GetJointSurface(multiblock);
        
        // Determine the symmetry plane 
        multiblock_roi = roi->SplitBySymmetryPlane(multiblock_roi);
        
        delete roi;
        
        std::cout << "Calculating distance function of femur ... " << std::endl;
        distfilter = vtkSmartPointer<vtkComputeDistanceField>::New();
        distfilter->SetRefModel(multiblock_roi, 0);
        
    }
}

void KneeAPI::SetVirtualLandmarks(int modelid, std::string filename){
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str()); 
    reader->Update();
    
    if(reader->GetOutput() != NULL){ SetVirtualLandmarks(modelid, reader->GetOutput()->GetPoints()); }
}

void KneeAPI::SetVirtualLandmarks(int modelid, vtkSmartPointer<vtkPoints> lmarks){
    motionfilter->SetVirtualLandmarks(modelid, lmarks );
}

void KneeAPI::Calibrate(int modelid){
    motionfilter->Calibrate(modelid);
}


void KneeAPI::OpenMarkerFile(std::string filename, int modelid, std::string type){
    std::ifstream infile (filename);
    std::string line;

    int time=0;
    
    vtkSmartPointer<vtkGeneratePointCloudFilter> pcgenerator = vtkSmartPointer<vtkGeneratePointCloudFilter>::New();
    
    while (std::getline(infile, line))
    {
        std::vector<std::string> row_values;

        split(line, '\t', row_values);

        
        if(row_values.size() > 55 && time>0 ){
            // Force '.' as the radix point in the conversion string to double. 
            // Otherwise: const auto oldLocale=std::setlocale(LC_NUMERIC,nullptr);
            std::setlocale(LC_NUMERIC,"C");
            
            if(row_values[14]=="OK" &&  row_values[18]=="OK" && row_values[22]=="OK" && (row_values[26]=="OK"|| row_values[26]=="Used OOV") && row_values[33]=="OK" && row_values[43]=="OK" && row_values[47]=="OK" && row_values[51]=="OK" && row_values[55]=="OK" ){
            
            // Bone 1
            pcgenerator->AddTimeFrame(time-1, 0, std::stod(row_values[5]), std::stod(row_values[6]), std::stod(row_values[7]),std::stod(row_values[8]),std::stod(row_values[9]),std::stod(row_values[10]),std::stod(row_values[11]));

            // Bone 2
            pcgenerator->AddTimeFrame(time-1, 1, std::stod(row_values[34]), std::stod(row_values[35]), std::stod(row_values[36]), std::stod(row_values[37]), std::stod(row_values[38]),std::stod(row_values[39]),std::stod(row_values[40]));
                
            time++;
            }
        }
        
        if(time==0){time++;}
        
    }
     
    pcgenerator->Update();
    
    vtkSmartPointer<vtkPolyData> polydata = vtkPolyData::SafeDownCast( pcgenerator->GetOutput() );
       
    
    // Store the marker data in the calibration class object 
    
    if(type=="landmarks"){
        motionfilter->SetRealLandmarks(modelid, polydata);
    }
    else if(type=="surfprobe"){
        motionfilter->SetRealSurface(modelid, polydata);
    }
    
}

int KneeAPI::GetNumberOfTimeFrames(){
    return motionfilter->GetNumberOfTimeFrames();
}

void KneeAPI::OpenMotionFile(std::string filename){
    std::ifstream infile (filename);
    std::string line;

    int time=0;
    
    while (std::getline(infile, line))
    {
        std::vector<std::string> row_values;

        split(line, '\t', row_values);

        if(row_values.size() > 43 && time>0 ){
            // Force '.' as the radix point in the conversion string to double. 
            // Otherwise: const auto oldLocale=std::setlocale(LC_NUMERIC,nullptr);
            std::setlocale(LC_NUMERIC,"C");
            
            if(row_values[14]=="OK" &&  row_values[18]=="OK" && row_values[22]=="OK" && (row_values[26]=="OK"|| row_values[26]=="Used OOV") && row_values[33]=="OK" && row_values[43]=="OK" && row_values[47]=="OK" && row_values[51]=="OK" && row_values[55]=="OK" ){
            
            // Bone 1
            motionfilter->AddTimeFrame(time-1, 1, std::stod(row_values[5]), std::stod(row_values[6]), std::stod(row_values[7]),std::stod(row_values[8]),std::stod(row_values[9]),std::stod(row_values[10]),std::stod(row_values[11]));

            // Bone 2
            motionfilter->AddTimeFrame(time-1, 0, std::stod(row_values[34]), std::stod(row_values[35]), std::stod(row_values[36]), std::stod(row_values[37]), std::stod(row_values[38]),std::stod(row_values[39]),std::stod(row_values[40]));
            
            time++;
            
            }
        }
        if(time==0){time++;}  
    }
}

vtkSmartPointer<vtkMultiBlockDataSet> KneeAPI::GetRealModel(int time, bool roi_only){
    if(roi_only==true){ motionfilter->SetInputData(0, multiblock_roi); }
    else{ motionfilter->SetInputData(0, multiblock); }
    
    motionfilter->SetTimeFrame(time);
    motionfilter->Update();
    
    vtkSmartPointer<vtkMultiBlockDataSet> newshape = vtkMultiBlockDataSet::SafeDownCast(motionfilter->GetOutput());
    
    if(roi_only==true){
        distfilter->SetInputData(0, newshape->GetBlock(1));
        distfilter->Update();
        newshape->SetBlock(1, vtkPolyData::SafeDownCast(distfilter->GetOutput()));
    }
    
    return newshape;
}


vtkSmartPointer<vtkPolyData> KneeAPI::GetVirtualModel(int modelidx, bool roi_only){
    if(roi_only){
        return vtkPolyData::SafeDownCast(multiblock_roi->GetBlock(modelidx));
    }
    else{
        return vtkPolyData::SafeDownCast(multiblock->GetBlock(modelidx));
    }
}

std::vector<std::vector<double>> KneeAPI::CalculateDistances(std::vector<std::vector<int>>* points){
    
    // container for distance values
    std::vector<std::vector<double>> result(2);
    result[0] = std::vector<double>(0);
    result[1] = std::vector<double>(0);
    
    // container for point ids
    if(points != NULL){
        (*points)[0] = std::vector<int>(motionfilter->GetNumberOfTimeFrames());
        (*points)[1] = std::vector<int>(motionfilter->GetNumberOfTimeFrames());
    }
    
    for(int time=0; time<motionfilter->GetNumberOfTimeFrames(); time++){
        // Get the tibia model for this time frame with the distance data-array and the lateral-medial data-array
        vtkSmartPointer<vtkPolyData> tibia = vtkPolyData::SafeDownCast( vtkMultiBlockDataSet::SafeDownCast( GetRealModel( time, true) )->GetBlock(1));
        
        vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr->SetInputData(tibia);
        wr->SetFileName(("/home/jeroen/Documents/KNEE-project/Output/TibiaDistField_"+std::to_string(time)+".vtk").c_str());
        wr->Update();
        
        // Find minimum distance
        vtkSmartPointer<vtkDataArray> distarray =  tibia->GetPointData()->GetArray("SignedDistances");
        vtkSmartPointer<vtkDataArray> symmarray =  tibia->GetPointData()->GetArray("SymmPlane");
        
        std::vector<double> minvalues = {distarray->GetRange()[1], distarray->GetRange()[1]};
                
        for(int p=0; p<tibia->GetNumberOfPoints(); p++){
            double point[3]; tibia->GetPoint(p, point);
            int side = symmarray->GetComponent(p,0);
            if(distarray->GetComponent(p,0)<minvalues[side]){ 
                minvalues[side] = distarray->GetComponent(p,0);
                if(points != NULL){ (*points)[side][time] = p; }
            }
        }
        
        std::cout << "Min distance " << minvalues[0] << " " << minvalues[1] << std::endl;
        
        result[0].push_back(minvalues[0]);
        result[1].push_back(minvalues[1]);
    }
    
    return result;
}
