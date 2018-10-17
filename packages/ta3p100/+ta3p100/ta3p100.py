import datajoint as dj
import numpy as np
from collections import Counter
from cloudvolume import CloudVolume

schema = dj.schema('microns_ta3p100')
ta3p100 = dj.create_virtual_module('ta3p100', 'microns_ta3p100')

@schema
class Ingest(dj.Manual):
    definition = """
    # Pinky100 ingest table using Princeton's naming conventions
    cleft_segid           : bigint   # maps to synapse_id in ta3p100.Synapses()
    ---
    postsyn_basin=null    : bigint   # not used in the other tables
    postsyn_segid=null    : bigint   # maps to postsyn in ta3p100.Synapses()
    postsyn_x=null        : bigint   # not used in the other tables
    postsyn_y=null        : bigint   # not used in the other tables
    postsyn_z=null        : bigint   # not used in the other tables
    presyn_segid=null     : bigint   # maps to presyn in ta3p100.Synapses()
    presyn_basin=null     : bigint   # not used in the other tables
    presyn_x=null         : bigint   # not used in the other tables
    presyn_y=null         : bigint   # not used in the other tables
    presyn_z=null         : bigint   # not used in the other tables
    centroid_x            : bigint   # maps to synapse_x in ta3p100.Synapses()
    centroid_y            : bigint   # maps to synapse_y in ta3p100.Synapses()
    centroid_z            : bigint   # maps to synapse_z in ta3p100.Synapses()
    bbox_bx               : bigint   # maps to syn_bbox_x1 in ta3p100.Synapses()
    bbox_by               : bigint   # maps to syn_bbox_y1 in ta3p100.Synapses()
    bbox_bz               : bigint   # maps to syn_bbox_z1 in ta3p100.Synapses()
    bbox_ex               : bigint   # maps to syn_bbox_x2 in ta3p100.Synapses()
    bbox_ey               : bigint   # maps to syn_bbox_y2 in ta3p100.Synapses()
    bbox_ez               : bigint   # maps to syn_bbox_z2 in ta3p100.Synapses()
    size                  : bigint   # maps to size in ta3p100.Synapses()
    """

@schema
class Segment(dj.Manual):
    definition = """
    # Segment: a volumetric segmented object
    segment_id   : bigint   # unique segment id
    """
    
@schema
class Synapse(dj.Manual):
    definition = """
    # Anatomically localized synapse between two Segments from Pinky100
    synapse_id        : bigint   # synapse index within the segmentation 
    ---
    (presyn) -> ta3p100.Segment(segment_id)
    (postsyn) -> ta3p100.Segment
    synapse_x         : bigint   # (EM Voxels)
    synapse_y         : bigint   # (EM Voxels)
    synapse_z         : bigint   # (EM Voxels)
    syn_bbox_x1       : bigint   # (EM Voxels) - bounding box
    syn_bbox_y1       : bigint   # (EM Voxels) - bounding box
    syn_bbox_z1       : bigint   # (EM Voxels) - bounding box
    syn_bbox_x2       : bigint   # (EM Voxels) - bounding box
    syn_bbox_y2       : bigint   # (EM Voxels) - bounding box
    syn_bbox_z2       : bigint   # (EM Voxels) - bounding boxta3p100 = dj.create_virtual_module('ta3p100', 'microns_ta3p100')
    size              : bigint   #
    """

@schema    
class Mesh(dj.Imported):
    volume = CloudVolume('https://storage.googleapis.com/neuroglancer/pinky100_v0/seg/lost_no-random/bbox1_0', progress=False)
    
    definition = """
    # Mesh
    -> Segment
    ---
    n_vertices     : int        # number of vertices in this mesh
    n_triangles    : int        # number of triangles in this mesh
    vertices       : longblob   # x,y,z coordinates of vertices
    triangles      : longblob   # triangles (triplets of vertices)
    """
    
    def make(self, key):
        
        mesh = self.volume.mesh.get(key['segment_id'])
        vertices = np.asarray(mesh['vertices'])
        faces = mesh['faces']
        num_triangles = int(len(faces) / 3)
        triangles = np.reshape(np.array(faces), (num_triangles, 3))
        
        key['n_vertices'] = mesh['num_vertices']
        key['n_triangles'] = num_triangles
        key['vertices'] = vertices
        key['triangles'] = triangles
        
        self.insert1(key)

@schema
class QualityScore(dj.Computed):
    definition = """
    -> Mesh
    ---
    n_non_manifold_edges    : int   # Number of non-manifold edges the mesh contains
    """

    def non_manifold_count(self, triangles):
        mesh_edges = []
        for triangle in triangles:
            mesh_edges.extend([(triangle[0], triangle[1]),
                               (triangle[1], triangle[2]),
                               (triangle[2], triangle[0])])

        ordered_mesh_edges = [(a, b) if a < b else (b, a) for a, b in mesh_edges]

        count_of_edges = Counter(ordered_mesh_edges)

        score = 0
        for edge, count in count_of_edges.most_common():
            if count != 2:
                score += 1
        return score
        
    def make(self, key):
        
        triangles = (Mesh() & key).fetch1('triangles')
        key['n_non_manifold_edges'] = self.non_manifold_count(triangles)
        self.insert1(key)
