# EASE: EM-Assisted Sources Extraction. 

EASE is a toolbox for fusing calcium imaging data and densely reconstructed Electron Microscopy (EM) segments. 

## Requirements 
**Running environment**
1. <img src="https://upload.wikimedia.org/wikipedia/commons/2/21/Matlab_Logo.png" height="20" /> MATLAB: all computations. The earliest MATLAB version passing my test is R2015a. However, you need versions after R2015b if you have to use DataJoint to load EM meshes. 

2. <img src="https://upload.wikimedia.org/wikipedia/en/6/62/MySQL.svg" height="20"/> MySQL: store data. The other option is to save data into *.mat file directly.

**Pakages** 

1. DataJoint: access data from MySQL. It's not required if you don't need accessing data.  

**Data**

1. EM meshes 
2. calcium imaging videos. 
3. high-resolution 2-photon scanning of the imaged volume. 

details of the data information will be explained in the section of [Run EASE](run-ease); 

## Installation
1. clone the package (bash environment) 
```bash 
git clone --recurse-submodules  https://github.com/zhoupc/ease.git
```

2. setup the path (MATLAB environment)
```matlab
>> run ease_setup.m
```

## Run EASE 
We created a GUI for running EASE. Once we set up the data paths and configurations, just run 
```matlab 
>> ease.startGUI()
```
Then you can start processing data with the GUI.  Easy-peasy! 

We also provided many scripts for processing data automatically. More details of running EASE will be added soon. 



