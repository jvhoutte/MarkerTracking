# MarkerTracking
GUI to register a surface model to optically tracked markers and compute the joint space width during articulation.


![alt text](https://github.com/jvhoutte/MarkerTracking/blob/main/figures/gui.png?raw=true)

## Requirements

* VTK: 9.0
* Qt: 5
* CMake: >2.8

## The cmake way 
    
The commands to build the project would be:
``` 
    mkdir build_dir
    cd build_dir
    cmake .. 
    make
```

## Step-by-step instructions

* Load Femur/Tibia: load the digital surface models of the bones (.stl).
* The surface models are displayed in the window. Start annotating the same landmarks as were annotated on the real bones by the digital annotator. Annotate a landmark by hoovering over the rendered model and pressing "s" on the keyboard when the cursor is above the required position. To remove the last landmark press "d". Save the annotated landmarks by pressing "Save landmarks".
* Read surface probe: Load the optical data acquired by moving the digital annotator over the real bone's surface (.txt).
* Real landmarks: Load the optical data representing the landmarks on the real bones (.txt). 
* Calibrate: when all files are loaded succesfully, the calibration can be performed. 
* Load motion data: Load the optical data acquired during articulating motion (.txt). Only timeframes without occlusions are considered. 
* Calculate distances: Calculate the joint space width for all the timeframes. The minimal
distance for medial and lateral side are shown in the graph. The small render window
in the bottom right corner shows the distance field on the tibia joint surface with the
point of minimal distance highlighted by two spheres. Note that the colorscale is cut off
at zero and ten millimeters.
* Show ROI only: Alter between the full bone model visualisation and the region of
interest visualisation.

## Contact 

For any question, please contact:
Jeroen Van Houtte (jeroen.vanhoutte@uantwerpen.be)
