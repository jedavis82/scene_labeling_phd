# scene_labeling
This codebase contains the implementation of the Scene To Text system that incorporates spatial relationship information developed as PhD dissertation research in CITE PHD. 

This code is free to use but the authors ask that if you make use of any of the code during research you cite the work using PAPER CITATION. 

# Installation 
## Anaconda Environment
The codebase makes use of an Anaconda enviroment. This environment can be installed by running the following command from the conda_env directory in an Anaconda prompt: 

`conda env create -f environment.yml` 

## MatLab Engine for Python
The Histogram of Forces (HOF) code is used to compute the spatial relationships between object two tuples in an image. The HOF code is implemented in MatLab and the MatLab engine for Python is required to run the HOF code. Installation instructions for the MatLab engine for Python can be found at this [link](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).  You will need to follow the instructions for intalling the engine API at the system command prompt and this will need to be done inside of the Anaconda environment. 

For example, from an Anaconda prompt, run the following commands: 

`conda activate sensitivity_analysis`

`cd <matlabroot>\extern\engines\python`

`python setup.py install`
## YOLOv3 Object Detection Model 
The code base uses the YOLOv3 object detection model. Due to size constraints on the repository, this model could not be uploaded. The model can be downloaded from the [YOLOv3 site](https://pjreddie.com/darknet/yolo/). 

The files required are: 
- yolov3.cfg
- yolov3.weights

The code base originally stored these files in the input/models/ directory, as can be seen in the object_detection.py script. 

# Usage 

# Example System Output 

