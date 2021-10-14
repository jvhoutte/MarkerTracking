#ifndef KNEEAPI_H
#define KNEEAPI_H


#include "Defs.h"
#include "vtkMotionFilter.h"
#include "vtkComputeDistanceField.h"
#include "GeometricROI.h"
#include <clocale>

class KneeAPI
{
    public:
    
        /* Constructor */
        KneeAPI();
        
        /* Open virtual model file (.stl) */
        void openFile(int blockid, std::string filename);
    
        /* Set landmarks annotated on the virtual stl model. Can either be a vtkpoints or a vtk file */
        void SetVirtualLandmarks(int modelid, vtkSmartPointer<vtkPoints> lmarks);
        void SetVirtualLandmarks(int modelid, std::string filename);
        
        /* Open optical marker file */
        void OpenMarkerFile(std::string filename, int modelid, std::string type);
        
        /* Open optical motion data */
        void OpenMotionFile(std::string filename);
    
        /* Get virtual model */
        vtkSmartPointer<vtkPolyData> GetVirtualModel(int modelidx, bool roi_only=false);
        
        /* Get articulated real model at time t */
        vtkSmartPointer<vtkMultiBlockDataSet> GetRealModel(int time=0, bool roi_only=false);
        
        /* Peform calibration */
        void Calibrate(int modelid);        
        
        std::vector<std::vector<double>> CalculateDistances(std::vector<std::vector<int>>* points = NULL);
        
        int GetNumberOfTimeFrames();
        
    protected:
        
        void DetermineROI();
        
    private:
    
        vtkSmartPointer<vtkMotionFilter> motionfilter;
        vtkSmartPointer<vtkMultiBlockDataSet> multiblock;
        vtkSmartPointer<vtkMultiBlockDataSet> multiblock_roi;
        
        vtkSmartPointer<vtkComputeDistanceField> distfilter;
    
};



#endif



