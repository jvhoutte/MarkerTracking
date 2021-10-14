#ifndef __vtkTransformFilter_h
#define __vtkTransformFilter_h
 
#include <numeric>
#include "vtkDataObjectAlgorithm.h"
#include "vtkTrivialProducer.h"
#include "vtkPolyDataNormals.h"
#include "vtkMultiBlockDataSet.h"
#include "Defs.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include "vtkMatrix4x4.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkIterativeClosestPointTransform.h"

class CompoundTransform{
    public:
        CompoundTransform(){
            totaltransformmatrix = vtkSmartPointer<vtkMatrix4x4>::New();
            totaltransformmatrix->Identity();
        }
        
        void SetParameters(double q1, double q2, double q3, double q4, double t1, double t2, double t3){
            double quatarray[4] = {q1, q2, q3, q4};
            vtkQuaterniond quat(quatarray);
            double rotation[3][3]; quat.ToMatrix3x3(rotation);
            
            totaltransformmatrix->Identity();
            totaltransformmatrix->SetElement(0,0,rotation[0][0]);
            totaltransformmatrix->SetElement(0,1,rotation[0][1]);
            totaltransformmatrix->SetElement(0,2,rotation[0][2]);
            totaltransformmatrix->SetElement(1,0,rotation[1][0]);
            totaltransformmatrix->SetElement(1,1,rotation[1][1]);
            totaltransformmatrix->SetElement(1,2,rotation[1][2]);
            totaltransformmatrix->SetElement(2,0,rotation[2][0]);
            totaltransformmatrix->SetElement(2,1,rotation[2][1]);
            totaltransformmatrix->SetElement(2,2,rotation[2][2]);
            totaltransformmatrix->SetElement(0,3,t1);
            totaltransformmatrix->SetElement(1,3,t2);
            totaltransformmatrix->SetElement(2,3,t3);
        }
    
        vtkSmartPointer<vtkTransform> GetTransform(){
            vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
            transform->SetMatrix(totaltransformmatrix);
            return transform;
        }
        
        vtkSmartPointer<vtkTransform> GetInverseTransform(){
            vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
            transform->SetMatrix(totaltransformmatrix);
            transform->Inverse();
            return transform;
        }
    
    private:
        vtkSmartPointer<vtkMatrix4x4> totaltransformmatrix;
    
};
 
class MultiCompoundTransform{
  public:
    MultiCompoundTransform(int NUMBER_OF_COMPOUNDS=2){
        transforms.resize(NUMBER_OF_COMPOUNDS);
    }
    
    CompoundTransform& operator[](int idx) { return transforms[idx]; }
    const CompoundTransform& operator[](int idx) const { return transforms[idx]; }
    
  private:
    std::vector<CompoundTransform> transforms;
};


class CalibrationObject
{
public:
    CalibrationObject(){
        landmarks_virtual = vtkSmartPointer<vtkPoints>::New();
        landmarks_real = vtkSmartPointer<vtkPolyData>::New();
        model_real = vtkSmartPointer<vtkPolyData>::New();
        model_virtual = vtkSmartPointer<vtkPolyData>::New();
        
        calib = vtkSmartPointer<vtkTransform>::New();
        calib->Identity();
        
        invcalib = vtkSmartPointer<vtkTransform>::New();
        invcalib->Identity();
    }
    
    vtkSmartPointer<vtkPoints> landmarks_virtual;
    vtkSmartPointer<vtkPolyData> landmarks_real;
    vtkSmartPointer<vtkPolyData> model_real;
    vtkSmartPointer<vtkPolyData> model_virtual;
    
    vtkSmartPointer<vtkTransform> GetTransform(){return calib;}
    vtkSmartPointer<vtkTransform> GetInverseTransform(){return invcalib;}
    
    bool Calibrate(int modelidx){
        // Check if everything is ready for calibration 
        
        if(landmarks_real->GetNumberOfPoints()!=landmarks_virtual->GetNumberOfPoints() ){std::cout << "error 1..." << std::endl;}
        if(landmarks_real->GetNumberOfPoints()<3){std::cout << "error 2..." << std::endl;}
        if(model_real->GetNumberOfPoints()<100){std::cout << "error 3..." << std::endl;}
        if(model_virtual->GetNumberOfPoints()<100){std::cout << "error 4..." << std::endl;}
            
        if(landmarks_real->GetNumberOfPoints()==landmarks_virtual->GetNumberOfPoints() && landmarks_real->GetNumberOfPoints()>3 && model_real->GetNumberOfPoints()>100 && model_virtual->GetNumberOfPoints()>100){
            std::cout << "Calibrating..." << std::endl;
            Registration(modelidx);
            return true;
        }
        else{
            return false;
        }
    }
    
protected:
    void Registration(int modelidx){
        vtkSmartPointer<vtkLandmarkTransform> icp_lm = vtkSmartPointer<vtkLandmarkTransform>::New();
        icp_lm->SetSourceLandmarks(landmarks_virtual);
        icp_lm->SetTargetLandmarks(landmarks_real->GetPoints());
        icp_lm->SetModeToRigidBody(); // rotation and translation only
        icp_lm->Update();
        
        // Transform the source points by the ICP solution
        vtkSmartPointer<vtkTransformPolyDataFilter> ppmTransformFilter =  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        ppmTransformFilter->SetInputData(model_virtual);
        ppmTransformFilter->SetTransform(icp_lm);
        ppmTransformFilter->Update();
        
        
        vtkSmartPointer<vtkIterativeClosestPointTransform> icp = vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
        icp->SetTarget(ppmTransformFilter->GetOutput());
        icp->SetSource(model_real);
        icp->SetMaximumNumberOfLandmarks(model_real->GetNumberOfPoints());
        icp->GetLandmarkTransform()->SetModeToRigidBody();
        icp->SetMaximumNumberOfIterations(1000);
        icp->StartByMatchingCentroidsOff();
        icp->Update();
        
        // ICP->Inverse does not return the inverse of the calculated transformtion, but switches source and target. Instead, extract the transform-matrix and calculate inverse afterwards
        vtkSmartPointer<vtkMatrix4x4> icpmatrix = icp->GetMatrix();
        vtkSmartPointer<vtkTransform> icptransform = vtkSmartPointer<vtkTransform>::New();
        icptransform->SetMatrix(icpmatrix);
        icptransform->Inverse();
        
        //icp->Inverse();
        
        calib->Identity();
        calib->PostMultiply();
        calib->Concatenate(icp_lm);
        calib->Concatenate(icptransform);
        
        invcalib->DeepCopy(calib);
        invcalib->Inverse();
        
        vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformfilter->SetInputData(model_virtual);
        transformfilter->SetTransform(calib);
        transformfilter->Update();
        
        vtkSmartPointer<vtkPolyDataWriter> wr1 = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr1->SetFileName(("/home/jeroen/Documents/KNEE-project/Output/Target_reg_"+std::to_string(modelidx)+".vtk").c_str());
        wr1->SetInputData(model_real);
        wr1->Update();
        
        vtkSmartPointer<vtkPolyDataWriter> wr2 = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr2->SetFileName(("/home/jeroen/Documents/KNEE-project/Output/Result_reg_"+std::to_string(modelidx)+".vtk").c_str());
        wr2->SetInputData(transformfilter->GetOutput());
        wr2->Update();
        
        vtkSmartPointer<vtkPolyDataWriter> wr3 = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr3->SetFileName(("/home/jeroen/Documents/KNEE-project/Output/Init_reg_"+std::to_string(modelidx)+".vtk").c_str());
        wr3->SetInputData(ppmTransformFilter->GetOutput());
        wr3->Update();
        
        std::cout << "done " << std::endl;
    }
    
    
    
private:
    
    vtkSmartPointer<vtkTransform> invcalib;
    vtkSmartPointer<vtkTransform> calib;
};

class vtkMotionFilter : public vtkDataObjectAlgorithm 
{
public:
  vtkTypeMacro(vtkMotionFilter,vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkMotionFilter *New();
  
  vtkSmartPointer<vtkTransform> GetBlockTransform(int blockid);
  
  void SetTimeFrame(int frameidx){timeframe = frameidx; this->Modified();}
  
  void AddTimeFrame(int time, int model, double q1, double q2, double q3, double q4, double t1, double t2, double t3);
  
  void SetVirtualLandmarks(int modelidx, vtkSmartPointer<vtkPoints> lm_virt);
  void SetRealLandmarks(int modelidx, vtkSmartPointer<vtkPolyData> lm_real);
  void SetVirtualSurface(int modelidx, vtkSmartPointer<vtkPolyData> surf_virt);
  void SetRealSurface(int modelidx, vtkSmartPointer<vtkPolyData> surf_real);
  
  bool Calibrate(int modelidx);
  
  int GetNumberOfTimeFrames(){return motiondata.size();}
  
protected:
  vtkMotionFilter();
  ~vtkMotionFilter();
 
  int FillInputPortInformation( int port, vtkInformation* info ) VTK_OVERRIDE;
  int FillOutputPortInformation(int portNumber, vtkInformation *info) VTK_OVERRIDE;
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
private:
  vtkMotionFilter(const vtkMotionFilter&);  // Not implemented.
  void operator=(const vtkMotionFilter&);  // Not implemented.

  int timeframe = 0;
  
  std::vector<MultiCompoundTransform> motiondata;
  
  std::vector<CalibrationObject> calibrationdata;
  
  vtkSmartPointer<vtkMultiBlockDataSet> surfaceRest;
  vtkSmartPointer<vtkMultiBlockDataSet> surfacePose;
};



class vtkGeneratePointCloudFilter : public vtkDataObjectAlgorithm 
{
public:
  vtkTypeMacro(vtkGeneratePointCloudFilter,vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkGeneratePointCloudFilter *New();
  
  void SetTimeFrame(int frameidx){timeframe = frameidx; this->Modified();}
  
  void AddTimeFrame(int time, int model, double q1, double q2, double q3, double q4, double t1, double t2, double t3);
  
  
  
protected:
  vtkGeneratePointCloudFilter();
  ~vtkGeneratePointCloudFilter();
 
  int FillOutputPortInformation(int portNumber, vtkInformation *info) VTK_OVERRIDE;
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
private:
  vtkGeneratePointCloudFilter(const vtkGeneratePointCloudFilter&);  // Not implemented.
  void operator=(const vtkGeneratePointCloudFilter&);  // Not implemented.

  int timeframe = 0;
  
  std::vector<MultiCompoundTransform> framedata;
  
  vtkSmartPointer<vtkPolyData> pointcloud;
};




#endif
