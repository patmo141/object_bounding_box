bl_info = {
    "name": "Object Bounding Box",
    "author": "Patrick R. Moore",
    "version": (0, 1),
    "blender": (2, 7, 3),
    "location": "View3D > Add > Mesh > New Object",
    "description": "Adds new cube which is minumum bounding box",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh"}

import bpy
import bmesh
import math
import random
import time
from mathutils import Vector, Matrix
from bpy.props import BoolProperty, FloatProperty, IntProperty
import numpy as np

def bbox_orient(bme_verts, mx):
    '''
    takes a lsit of BMverts
    '''
    verts = [mx * v.co for v in bme_verts]
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    zs = [v[2] for v in verts]
    
    return (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs))

def bbox_vol(box):
    
    V = (box[1]-box[0]) * (box[3]-box[2]) * (box[5]-box[4])
    
    return V

def box_cords(box):
    '''
    returns vertices in same configuration as default cube in blender
    easy to asign v.co of a cube primitive
    '''
    cords = [Vector((box[0],box[2],box[4])),
             Vector((box[0],box[2],box[5])),
             Vector((box[0],box[3],box[4])),
             Vector((box[0],box[3],box[5])),
             Vector((box[1],box[2],box[4])),
             Vector((box[1],box[2],box[5])),
             Vector((box[1],box[3],box[4])),
             Vector((box[1],box[3],box[5])),
             ]
    
    return cords

def main(context, rand_sample, spin_res, make_sphere):
    start = time.time()
    #rand_sample = 400  #randomly select this many directions on a solid hemisphere to measure from
    #spin_res = 180   #180 steps is 0.5 degrees
    
    world_mx = context.object.matrix_world
    me = context.object.data
    bme = bmesh.new()
    bme.from_mesh(me)
    
    convex_hull  = bmesh.ops.convex_hull(bme, input = bme.verts, use_existing_faces = True)
    total_hull = convex_hull['geom']
    
    hull_verts = [item for item in total_hull if hasattr(item, 'co')]
    hull_faces = [item for item in total_hull if hasattr(item, 'no')]
    hull_bme = bmesh.new()
    #hull_bme.verts = hull_verts
    #hull_bme.faces = hull_faces
    
    min_mx = Matrix.Identity(4)
    min_box = bbox_orient(hull_verts, min_mx)
    min_V = bbox_vol(min_box)
    print('initial volume %f' % min_V)
    min_axis = Vector((0,0,1))
    min_angle = 0
    axes = []
    for i in range(0,rand_sample):
        u = random.random()
        v = random.random()
        
        theta = math.pi * u  
        phi =  math.acos(2 * v - 1)
        
        x = math.cos(theta) * math.sin(phi)
        y = math.sin(theta) * math.sin(phi)  
        z = math.cos(phi)
        
        axis = Vector((x,y,z))
        axes.append(axis)
        for n in range(0, spin_res):
            angle = math.pi/2 * n/spin_res
            rot_mx = Matrix.Rotation(angle,4,axis)
            
            box = bbox_orient(hull_verts, rot_mx)
            test_V = bbox_vol(box)
            
            if test_V < min_V:
                min_V = test_V
                min_axis = axis
                min_angle = angle
                min_box = box
                min_mx = rot_mx
    elapsed_time = time.time() - start
    print('did %i iterations in %f seconds' % (rand_sample*spin_res, elapsed_time))
    print("final volume %f" % bbox_vol(min_box))         
    box_verts = box_cords(min_box)
    bpy.ops.mesh.primitive_cube_add()
    context.object.matrix_world = min_mx.inverted()*world_mx
    context.object.draw_type = 'BOUNDS'
    for i, v in enumerate(box_verts):
        context.object.data.vertices[i].co = v
    
    #visualize the sample vectors
    if make_sphere:
        sample_sphere = bmesh.new()
        for ax in axes:
            sample_sphere.verts.new(ax)
        
        sphere_me = bpy.data.meshes.new('Bound Samples')
        dest_ob = bpy.data.objects.new('Bound Samples',sphere_me)
        sample_sphere.to_mesh(sphere_me)
        context.scene.objects.link(dest_ob)
        sample_sphere.free()
    
    bme.free() 

def main_SVD(context, down_sample, method, spin_res):
    start = time.time()
    
    world_mx = context.object.matrix_world
    me = context.object.data
    bme = bmesh.new()
    bme.from_mesh(me)
    
    convex_hull  = bmesh.ops.convex_hull(bme, input = bme.verts, use_existing_faces = True)
    total_hull = convex_hull['geom']
    
    hull_verts = [item for item in total_hull if hasattr(item, 'co')]
    hull_faces = [item for item in total_hull if hasattr(item, 'no')]
    #hull_bme = bmesh.new()
    #hull_bme.verts = hull_verts
    #hull_bme.faces = hull_faces
    
    
    vert_data = [v.co for v in hull_verts]  #ToDo...world coords better?
    v0 = np.array(vert_data, dtype=np.float64, copy=True)
    ndims = v0.shape[0]
    
    #move data to origin
    t0 = -np.mean(v0, axis=1)
    M0 = np.identity(ndims+1)
    M0[:ndims, ndims] = t0
    v0 += t0.reshape(ndims, 1)
    
    #A = np.zeros(shape = [3,len(vert_data)])
    #for i in range(0,len(vert_data)):
    #        V1 = vert_data[i]
    #        A[0][i], A[1][i], A[2][i] = V1[0], V1[1], V1[2]
    
    
    U, s, V = np.linalg.svd(v0, full_matrices=True)
    
    print(V) #these should be eigenvectors which are the principle axes?
    
    rmx = Matrix.Identity(4)  #may need to transpose this stuff
    rmx[0][0], rmx[0][1], rmx[0][2] = V[0][0], V[0][1], V[0][2]
    rmx[1][0], rmx[1][1], rmx[1][2] = V[1][0], V[1][1], V[1][2]
    rmx[2][0], rmx[2][1], rmx[2][2] = V[2][0], V[2][1], V[2][2]
    
    
    min_box = bbox_orient(hull_verts, rmx)
    min_vol = bbox_vol(min_box)
    min_angle = 0
    min_mx = rmx
    
    if method == 1:  #PCAmax, spin around the biggest PCA component, good for linear objects
        axis = Vector((V[0][0], V[0][1], V[0][2]))
    elif method == 2:
        axis = Vector((V[1][0], V[1][1], V[1][2]))
    else:   
        axis = Vector((V[2][0], V[2][1], V[2][2]))
    
    for n in range(0, spin_res):
        print('did spin it')
        angle = math.pi/2 * n/spin_res
        rot_mx = Matrix.Rotation(angle,4,axis)
            
        box = bbox_orient(hull_verts, rot_mx)
        test_V = bbox_vol(box)
            
        if test_V < min_vol:
            print('found better volume with PCA')
            min_angle = angle
            min_box = box
            min_mx = rot_mx
    
    elapsed_time = time.time() - start
    
    print('found bbox of volume %f in %f seconds with SVD' % (min_vol, elapsed_time))
    
    box_verts = box_cords(min_box)
    bpy.ops.mesh.primitive_cube_add()
    context.object.matrix_world = min_mx.inverted() * world_mx
    context.object.draw_type = 'BOUNDS'
    for i, v in enumerate(box_verts):
        context.object.data.vertices[i].co = v
    
    bme.free()     
    
class ObjectMinBoundBox(bpy.types.Operator):
    """Find approximate minimum bounding box of object"""
    bl_idname = "object.min_bounds"
    bl_label = "Min Bounding Box"

    # generic transform props
    sample_vis = BoolProperty(
            name="Visualize Sample",
            description = 'add a sphere to the scene showing random direction sample',
            default=False,
            )
    area_sample = IntProperty(
            name="Direction Samples",
            description = 'number of random directions to test calipers in',
            default = 200)
    angular_sample = IntProperty(
            name="Direction samples",
            description = 'angular step to rotate calipers 90 = 1 degree steps, 180 = 1/2 degree steps',
            default = 90)
    
    method = IntProperty(
            name="Method",
            description = 'dummy prop.  0 is brute force, 1 is PCAmax, 2 is PCAmin',
            default = 0)
    
    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.active_object.type == 'MESH'

    
    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


    def draw(self, context):
        layout = self.layout
        row =layout.row()
        row.prop(self, "sample_vis")
        
        row =layout.row()
        row.prop(self, "area_sample")
        
        row =layout.row()
        row.prop(self, "angular_sample")
        
        #row =layout.row()
        #row.prop(self, "method")
        pass
    
    
    def execute(self, context):
        if self.method == 0:
            main(context, self.area_sample, self.angular_sample, self.sample_vis)
        
        else:
            main_SVD(context, 1, self.method, self.angular_sample)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(ObjectMinBoundBox)


def unregister():
    bpy.utils.unregister_class(ObjectMinBoundBox)


if __name__ == "__main__":
    register()
