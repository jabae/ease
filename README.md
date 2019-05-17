# EASE: EM-Assisted Sources Extraction. 
EASE is a toolbox for fusing calcium imaging data and densely reconstructed Electron Microscopy (EM) segments. 
![](graph_abstract.png)

## Requirements 
**Running environment**
1. <img src="https://upload.wikimedia.org/wikipedia/commons/2/21/Matlab_Logo.png" height="20" /> MATLAB: all computations. The earliest MATLAB version passing my test is R2015a. However, you need versions after R2015b if you have to use DataJoint to load EM meshes. 

2. <img src="https://upload.wikimedia.org/wikipedia/en/6/62/MySQL.svg" height="20"/> MySQL: store data. We mainly use MySQL to store EM meshes and other large intermediate results relating to preprocessing EM segments. We can also skip using MySQL if we have EM footprints on each 2P scanning plane.   

**Pakages** 

1. DataJoint: access data from MySQL. It's not required if you don't need accessing data.  

**Data**

1. EM meshes 
2. calcium imaging videos. 
3. high-resolution 2-photon scanning of the imaged volume. 

details of the data information will be explained in the section of [Run EASE](run-ease); 

## Installation
1. install [F-image](https://github.com/zhoupc/F-image) and then install EASE with F-image in MATLAB 
   ```matlab 
    >> fi.install('ease')
   ```
2. [install mysql](https://dev.mysql.com/doc/refman/8.0/en/installing.html) and configure the database. You also need to test your [datajoint](https://docs.datajoint.io/matlab/admin/Admin.html) connection with MySQL
## Use EASE 
* We described a typical pipeline of running EASE. Please check it [here](./how_to_use.md).
* We also uploaded our running scripts for producing results in our [EASE paper](xxx). Please check it [here](https://github.com/zhoupc/ease_project)

## References 
## Copyright 
[Pengcheng Zhou](https://zhoupc.github.io) @Columbia University, 2019



