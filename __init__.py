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
from bpy.props import BoolProperty, FloatProperty, IntProperty, EnumProperty
import numpy as np

def bbox_orient(bme_verts, mx):
    '''
    takes a lsit of BMverts ora  list of vectors
    '''
    if hasattr(bme_verts[0], 'co'):
        verts = [mx * v.co for v in bme_verts]
    else:
        verts = [mx * v for v in bme_verts]
        
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
    scale = world_mx.to_scale()
    trans = world_mx.to_translation()
    
    tr_mx = Matrix.Identity(4)
    sc_mx = Matrix.Identity(4)
    
    tr_mx[0][3], tr_mx[1][3], tr_mx[2][3] = trans[0], trans[1], trans[2]
    sc_mx[0][0], sc_mx[1][1], sc_mx[2][2] = scale[0], scale[1], scale[2]
    r_mx = world_mx.to_quaternion().to_matrix().to_4x4()
    
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
    
    fmx = tr_mx * r_mx * min_mx.inverted() * sc_mx
    context.object.matrix_world = fmx
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

def main_SVD(context, down_sample, method, spin_res, make_box):
    start = time.time()
    
    world_mx = context.object.matrix_world
    scale = world_mx.to_scale()
    trans = world_mx.to_translation()
    
    tr_mx = Matrix.Identity(4)
    sc_mx = Matrix.Identity(4)
    
    tr_mx[0][3], tr_mx[1][3], tr_mx[2][3] = trans[0], trans[1], trans[2]
    sc_mx[0][0], sc_mx[1][1], sc_mx[2][2] = scale[0], scale[1], scale[2]
    r_mx = world_mx.to_quaternion().to_matrix().to_4x4()
    
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
    
    U, s, V = np.linalg.svd(v0, full_matrices=True)
    
    #make a rotation matrix from eigenvectors (easy)
    rmx = Matrix.Identity(4)
    rmx[0][0], rmx[0][1], rmx[0][2] = V[0][0], V[0][1], V[0][2]
    rmx[1][0], rmx[1][1], rmx[1][2] = V[1][0], V[1][1], V[1][2]
    rmx[2][0], rmx[2][1], rmx[2][2] = V[2][0], V[2][1], V[2][2]
    
    
    min_box = bbox_orient(vert_data, rmx)
    min_vol = bbox_vol(min_box)
    min_angle = 0
    min_mx = rmx
      
    #these are our PCA directions
    X = Vector(( V[0][0], V[0][1], V[0][2]))
    Y = Vector(( V[1][0], V[1][1], V[1][2]))
    Z = Vector(( V[2][0], V[2][1], V[2][2]))
    
    for n in range(0, 2 * spin_res):
        angle = math.pi * n/(2 * spin_res)
        
        rmx = Matrix.Identity(4)
            
        if method == 'pca_x':  #keep x axis and rotate around it
            rmx[0][0], rmx[0][1], rmx[0][2] = X[0], X[1], X[2]
            
            y = math.cos(angle) * Y + math.sin(angle) * Z
            y.normalize()
            rmx[1][0], rmx[1][1], rmx[1][2] = y[0], y[1], y[2]
            
            z = -math.sin(angle) * Y + math.cos(angle) * Z
            z.normalize()
            rmx[2][0], rmx[2][1], rmx[2][2] = z[0], z[1], z[2]
            
        elif method == 'pca_y': #keep y axis and rotate around it
            x = math.cos(angle) * X - math.sin(angle) * Z
            x.normalize()
            rmx[0][0], rmx[0][1], rmx[0][2] = x[0], x[1], x[2]
            #keep y
            rmx[1][0], rmx[1][1], rmx[1][2] = Y[0], Y[1], Y[2]
        
            z = math.sin(angle) * X + math.cos(angle) * Z
            z.normalize()
            rmx[2][0], rmx[2][1], rmx[2][2] = z[0], z[1], z[2]
            
        else:
            x = math.cos(angle) * X + math.sin(angle) * Y
            x.normalize()
            rmx[0][0], rmx[0][1], rmx[0][2] = x[0], x[1], x[2]
            
            y = -math.sin(angle) * X + math.cos(angle) * Y
            y.normalize()
            rmx[1][0], rmx[1][1], rmx[1][2] = y[0], y[1], y[2]
            
            #Keep Z
            rmx[2][0], rmx[2][1], rmx[2][2] = Z[0], Z[1], Z[2]

        box = bbox_orient(vert_data, rmx)
        test_V = bbox_vol(box)
        if test_V < min_vol:
            min_angle = angle
            min_box = box
            min_mx = rmx
            min_vol = test_V
            
        if make_box:
            box_verts = box_cords(box)
            bpy.ops.mesh.primitive_cube_add()
        
            context.object.matrix_world =rmx.transposed().inverted() * world_mx
            context.object.draw_type = 'BOUNDS'
            for i, v in enumerate(box_verts):
                context.object.data.vertices[i].co = v    
        
    
    elapsed_time = time.time() - start
    
    print('found bbox of volume %f in %f seconds with SVD followed by rotating calipers' % (min_vol, elapsed_time))
    
    box_verts = box_cords(min_box)
    bpy.ops.mesh.primitive_cube_add()
    #FinalMatrix = TranslationMatrix * RotationMatrix * ScaleMatrix
    fmx = tr_mx * r_mx * min_mx.inverted() * sc_mx
    
    context.object.matrix_world =  fmx
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
    
    make_box = BoolProperty(
            name="Visualize Boxes",
            description = 'add a cube for all bounding boxes tried.  VERY MESS!',
            default=False,
            )
    area_sample = IntProperty(
            name="Direction Samples",
            description = 'number of random directions to test calipers in',
            default = 200)
    angular_sample = IntProperty(
            name="Direction samples",
            description = 'angular step to rotate calipers 90 = 1 degree steps, 180 = 1/2 degree steps',
            default = 50)
    
    #(identifier, name, description, icon, number)
    method_enum = [('brute_force', "BRUTE FORCE", 'Checks a bunch of random boxes'),
                   ('pca_y', 'PCAY','Good for linear objects'),
                   ('pca_x', 'PCAX', 'Good for flat objects'),
                   ('pca_z', 'PCAZ', 'Good for some things')]
        
    method = bpy.props.EnumProperty(
        name="Method", 
        description="Min BBox method to use", 
        items=method_enum, 
        default='brute_force',
        options={'ANIMATABLE'})
    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.active_object.type == 'MESH'

    
    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


    def draw(self, context):
        layout = self.layout
        row =layout.row()
        row.prop(self, "sample_vis")
        row.prop(self, "make_box")
        
        row =layout.row()
        row.prop(self, "area_sample")
        
        row =layout.row()
        row.prop(self, "angular_sample")
        
        row =layout.row()
        row.prop(self, "method")
        pass
    
    
    def execute(self, context):
        if self.method == 'brute_force':
            main(context, self.area_sample, self.angular_sample, self.sample_vis)
        
        else:
            main_SVD(context, 1, self.method, self.angular_sample, self.make_box)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(ObjectMinBoundBox)


def unregister():
    bpy.utils.unregister_class(ObjectMinBoundBox)


if __name__ == "__main__":
    register()
