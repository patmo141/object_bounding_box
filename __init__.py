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
    rand_sample = 400  #randomly select this many directions on a solid hemisphere to measure from
    spin_res = 180   #180 steps is 0.5 degrees
    
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
        pass
    
    
    def execute(self, context):
        main(context, self.area_sample, self.angular_sample, self.sample_vis)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(ObjectMinBoundBox)


def unregister():
    bpy.utils.unregister_class(ObjectMinBoundBox)


if __name__ == "__main__":
    register()
