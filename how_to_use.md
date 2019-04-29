### step 0: import EASE 
```matlab
>> fi.usepkg('ease'); 
```

### step 1: (first time use) create a project folder to store data/code/results 
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
datasets_list: [pinky40, pinky100]
databases_list: ['xxxxxx.us-east-1.rds.amazonaws.com:1111', '127.0.0.1:1111']
pinky40: {datajoint_name: ta3, rel_mesh: ta3.MeshFragment}
pinky100: {datajoint_name: ta3p100, rel_mesh: ta3p100.Mesh}
```

This metafile tells EASE that this project has two datasets named 'pinky40' and 'pinky100' and two available database hosts. **'datajoint_name'** provides the database name of the selected data and **'rel_mesh'** tells the table name of the EM meshes.  
</p>
</details>

**When you want to add new datasets or database hosts, don't forget to update the metainfo.yaml file.**

### step 2: add a dataset 
```matlab 
>> ease.add_dataset(); 
```
<details><summary>CLICK FOR details </summary>
<p>

```matlab 
>> ease.add_dataset()
name of the datasets: test
database name (e.g. ta3): ta3test
table of the EM meshes (e.g., ta3.Mesh): ta3test.Mesh 
test: dataset added. Here are things you need to do: 
	1. add data files to folder: /data/home/zhoupc/github/new_ease_project/data/test
	2. modify data options: /data/home/zhoupc/github/new_ease_project/test_config.yaml
	3. add database schema to access ta3test
```
**I. required data**:

0. **EM meshes**: It's saved into a database. By default, it contains at least 3 tables 
   (1). **Segmentation**: specifies the version of the EM segmentation. there are two columns: segmentation; segmentation_description
   (2). **Segment**: specify a unique key value for each EM segment. It has two fields: segmentation, segment_id
   (3) **Mesh**: the meshes for all EM segments. Beside the unique key value from **Segment**, it has the following fields:  n_vertices (bigint), n_triangles (bigint), vertices (longblob) and triangles(longblob).
 

1. **stack_2p.mat**: there is a variable stack_2p with a dimension of [d1, d2, d3] to represent the high-resolution 2p stack data. 

2. **functional_data.mat**: info relating to the calcium imaging video data. it has a cell array ([num_scans, num_planes, num_blocks]) of data loaders for each video data: 

3. registration.csv: xyz coordinates of a list of registered points in both 2p space and EM space. [x_2p, y_2p, z_2p, x_em, y_em, z_em] 

**II. important configurations to modify**
1. video_Fs: frame rate 
2. use_denoise: boolean; use the denoised result or not. (default: false)
3. num_scans: number of scans 
4. num_slices: number of planes 
5. num_blocks: number of blocks for one plane. 
6. dims_stack: 1*3 vector. dimension of the stack data. 
7. dims_video: 1*2 vector. dimension of each 2D video. 
8. range_2p: 1*3 vector. spatial range of the stack volume in physical space (unit: um). 
 </p>
</details>

 
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