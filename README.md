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
>> run_ease; 
```
Then you can start processing data with the GUI.  Easy-peasy! 

We also provided many scripts for processing data automatically. More details of running EASE will be added soon. 

## Package structure
* **@EM2P**: define a class object for organzing data, options and high level functions. 
* **@MF3D**: define a class object (matrix factorization 3D) for running matrix factorization on 3D data. 
* **+erun**: (Ease run) a collection of scripts for processing data. 
* **functions**: matlab functions used by the package. 
* **GUI**:  GUI callbacks and layouts. 
* **scripts**: a collection of scripts. (I'm going to slowly move all scripts into the folder +erun). 
* **config_ease.yaml**: an example yaml file for storing configurations 
* **ease_setup.m**: setup ease. 
* **run_ease.m**: a demo script 

## Copyright 
[Pengcheng Zhou](https://zhoupc.github.io) @Columbia University, 2019



