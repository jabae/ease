### step 0: import EASE 
```matlab
>> fi.usepkg('ease'); 
```

### step 1: (first time use) create a project folder to store data/code/results 
EASE requires 
```matlab 
>> ease.init_project();  
```

<details><summary>CLICK FOR details</summary>
<p>

EASE will ask you to interactively choose a project folder, and then create the following subfolders and files. 
* /data/
* /scripts/
* /results/
* /Figures/
* /Videos/
* /metainfo.yaml: a yaml file storing a list of datasets, database hosts and some info relating to each dataset. 

Here is an example **metainfo.yaml**: 

```yaml 
datasets_list: [pinky100]
databases_list: ['127.0.0.1:3306']
example_data: {datajoint_name: ta3p100, database_name: microns_ta3p100, rel_mesh: ta3p100.Mesh}
```

This metafile tells EASE that this project has one dataset named 'pinky100' and one database server at 127.0.0.1:3306. **'database_name'** is the name of this dataset in database and **'datajoint_name'** is an alias used in MATLAB to access this database. **'rel_mesh'** tells the table name of the EM meshes.  

You can either manually edit this **metainfo.yaml** or do it in MATLAB command window. 
```matlab 
>> ease.add_database(); 
>> ease.add_dataset(); # see step 2 for details.  
```
**When you want to add new datasets or database hosts, don't forget to update the metainfo.yaml file.**
</p>
</details>


### step 2: (first time use) add a dataset 

```matlab 
>> ease.add_dataset()
name of the datasets: pinky100
database name (e.g., microns_ta3): microns_ta3p100
database alias (e.g. ta3): ta3p100
table of the EM meshes (e.g., ta3.Mesh): ta3p100.Mesh
pinky100: dataset info added. Next we need the actual data for this dataset.
	1: registration.csv. A csv file for aligning EM space and the CI space.
	2: calcium imaging videos.
	3: 2p structural imaging.
You can find the detailed description of these datasets on xxx.
Once you get everything ready, run 
	>>ease.import_data();
to import them.
```
This command requires you to provide some basic info of a dataset. 


The above command only requires some basic info of the dataset. As shown above, we also need to import the actural data. Here is the summary of the data required by EASE
#### required datasets 
1. EM meshes.
   - The data are saved into a SQL database.
   - each segment  is uniquely determined by a **SEGMENT_ID** & **SEGMENTATION** (segmentation version) 
   - The mesh of each segment is given by **VERTICES** & **FACES** (check [WIKI](https://en.wikipedia.org/wiki/Triangle_mesh) for details of trimesh)
2. coordinates of the registed points in both 2p space and EM space. 
   - the file is **registration.csv** 
   - each row is the coordinates of the one point $(x_{em}, y_{em}, z_{em},x_{2p},y_{2p},z_{2p},)$. The unit of the EM coordinates is **nm**, while the 2P unit is **um**.  
    
3. calcium imaging videos. 
   - a set of videos named like **scan1_slice2_block_3.mat**. Here **1,2,3** correspond to scan ID, slice ID & block ID respectively. 
   - Each mat file contains a video with the dimension **d1 x d2 x T**. 
   - Each scan imaged multiple slices at different z planes. The data were temporally divided into multiple blocks.   
4. 2p structural imaging stack. 
   - the file is **stack_2p.mat** (you can use other names)
   - the file contains a 3D high-resolution structural imaging of a large volume (k*d1 x k*d2 x d3, k is the spatial up-sampling factor). All calcium imaging vedios were collected within the volume.  

There are also some configurations you need to modify for adapting your dataset. You can run 
```matlab 
>> ease.modify_configs(); 
```
to edit them directly. 
<details><summary>CLICK for details </summary>
<p>

**important configurations to modify**
1. video_Fs: frame rate 
2. use_denoise: boolean; use the denoised result or not. (default: false)
3. num_scans: number of scans 
4. num_slices: number of planes 
5. num_blocks: number of blocks for one plane. 
6. dims_stack: 1*3 vector. dimension of the stack data. 
7. dims_video: 1*2 vector. dimension of each 2D video. 
8. range_2p: 1*3 vector. spatial range of the stack volume in physical space (unit: um). 

```yaml
yaml_path: /data/home/zhoupc/ease_test/pinky_config.yaml
data_folder: /data/home/zhoupc/ease_test/data/pinky
output_folder: /data/home/zhoupc/ease_test/results/pinky
denoised_folder: cropped_denoised_video
raw_folder: cropped_raw_video
matfile_video: functional_data.mat
matfile_stack: stack_2p.mat
registration_csv: registration.csv
matfile_transformation: coor_convert.mat
matfile_em: em_2.mat
video_shifts:
  ii: []
  jj: []
video_zvals: []
video_zvals_updated: []
video_Fs: 15.0
video_T: 0.0
use_denoise: false
stack_shifts: []
d1: []
d2: []
d3: []
extra_margin: 5.0
FOV: []
FOV_stack: []
align_max_zshift: 2.0
num_scans: 8.0
num_slices: 3.0
num_blocks: 3.0
dims_stack: [512.0, 512.0, 310.0]
dims_video: [256.0, 256.0]
range_2p: [400.0, 400.0, 310.0]
scan_id: 1.0
slice_id: 1.0
block_id: 1.0
options_init: {init_method: tf, order_statistics: l3norm, min_similarity: 0.6, clear_results: false,
  save_fig: false, K_candidate: 3000.0, show_fig: true, K_new: 50.0}
show_em_only: false
em_shifts: []
em_segmentation: 2.0
em_load_flag: true
em_zblur: 8.0
em_scale_factor: 0.001
score_method: corr
nam_show: cn

```
 </p>
</details>

**Once you have all the datasets, simply run the following command to import them**. 
```matlab 
>> ease.import_data()
step 1: choose the folder containing all calcium imaging video data.
step 2: choose the high resolution structural data
step 3: choose the csv file for aligning 2P space and EM space.
```

 
### step 3: run EASE by selecting a dataset and a database 
```matlab 
>> run_ease; 
```
Then you can follow the prompts to choose datasets and database hosts.
<details><summary>CLICK FOR EXAMPLE INPUTS</summary>
<p>

```matlab
>> run_ease

**********choose the data to use**********
1: pinky100
2: pinky40
********************************************
data ID: 1
you selected data pinky100

The configuration of the EASE environment has been updated from 
/data/home/zhoupc/github/new_ease_project/pinky100_config.yaml

************ SELECT A DATABASE ************
1: 127.0.0.1:1111
2: xxxx.us-east-1.rds.amazonaws.com:1111

database ID: 1

**************** Done ********************
You are going to connect to a database
	127.0.0.1:1111.
Now type your login information
username: root
password: *******

  connection_id() 
 +---------------+
  21              


ans = 

  Connection with properties:

             host: '127.0.0.1:1111'
             user: 'root'
        initQuery: ''
    inTransaction: 0
           connId: 0
         packages: [0×1 containers.Map]
      foreignKeys: [0×0 struct]
      isConnected: 1

Database connected
```
</p>
</details>

### step 4: process data of one scan 
```matlab
ease.scan_id = 1; 
ease.block_id = 1; 
neuron = ease.get_MF3D();   % create a wrapper to run 3D matrix factorization 
K_new = [50, 50, 50];       % add 50 neurons in each initialization (3 iterations in total)
K_candidate = [3000, 4000, 5000];  % choose the top-k candidate EM components 
nb = [1, 2, 3];    % number of background components in each iteration 
min_ranks = [2, 2, 2];    % rematch EM component if the current EM is not among the top-2 match
configs = ease.create_running_configs(K_new, K_candidate, nb, min_ranks);  % create a configuration to run EASE automatically 

neuron.at_ease(configs);   % your computer is working very hard, so it's your coffee time. 
```

### step 5: visualize results and do manual interventions
```matlab 
ease.startGUI(); 
```
There are short keyboard shortcuts for the GUI. 
* n: show next neuron 
* b: show previous neuron 
* s: label the neuron as **soma**
* d: label the neuron as **dendrite**
* 0: mark this neuron to be deleted later. 