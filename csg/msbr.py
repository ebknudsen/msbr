import openmc
import pandas as pd
import numpy as np
import math

import materials

in2cm=2.54
cm2in=1.0/2.54


class msbr(openmc.Model):
  def __init__(self):
    self.materials=mats=materials.define_materials()
    self.md={m.name:m for m in self.materials}

  def geom(self):
    zoneIA_stringer=self.zoneIA_stringer()
    zoneIIA_stringer=self.zoneIIA_stringer()
    zoneA_lattice = openmc.RectLattice()
    zoneA_lattice.lower_left=(-4*in2cm*23/2.0,-4*in2cm*23/2.0)
    zoneA_lattice.pitch=(4.0*in2cm,4.0*in2cm)
    zoneA_lattice.universes=np.tile(zoneIA_stringer,(23,23))
    zoneA_lattice.universes[5,5]=zoneIIA_stringer
    zoneA_lattice.universes[4,5]=zoneIIA_stringer
    zoneA_lattice.universes[5,4]=zoneIIA_stringer
    zoneA_lattice.universes[6,5]=zoneIIA_stringer
    zoneA_lattice.universes[5,6]=zoneIIA_stringer

    zoneI_IIB_boundary=self.stringer_boundary()
    
    bound=openmc.ZCylinder(r=50*in2cm,boundary_type='vacuum')

    core=openmc.Cell(fill=zoneA_lattice,region=zoneI_IIB_boundary)
    self.geometry=openmc.Geometry([core])
    self.geometry.export_to_xml()

  def stringer_boundary(self):
    diagonal_cut_left_back = openmc.Plane(a=1,b=-1, c=0, d= -4*in2cm*math.sqrt(2.0)*5)
    diagonal_cut_left_front = openmc.Plane(a=1,b=1, c=0, d= -4*in2cm*math.sqrt(2.0)*5)
    diagonal_cut_right_back = openmc.Plane(a=-1,b=-1, c=0, d= -4*in2cm*math.sqrt(2.0)*5)
    diagonal_cut_right_front = openmc.Plane(a=-1,b=1, c=0, d= -4*in2cm*math.sqrt(2.0)*5)

    straight_cut_left =   openmc.XPlane(x0=-4*in2cm*5)
    straight_cut_right =  openmc.XPlane(x0=4*in2cm*5)
    straight_cut_front =  openmc.YPlane(y0=-4*in2cm*5)
    straight_cut_back =   openmc.YPlane(y0=4*in2cm*5)

    top = openmc.ZPlane(z0= 171*in2cm/2.0)
    bot = openmc.ZPlane(z0=-171*in2cm/2.0)

    boundary = -top & +bot & +diagonal_cut_upper_left & +diagonal_cut_lower_left & +diagonal_cut_upper_right & +diagonal_cut_lower_right & +straight_cut_left & -straight_cut_right & +straight_cut_front & -straight_cut_back
    return boundary

  def stringer_sq(self, bore_radius=1.34/2.0):
    #return a zoneIA stringer as a universe or as a list of cells
    str_zIA_length=171*in2cm
    str_zI_side=3.9*in2cm
    str_zI_dot_r=0.1*in2cm/2.0
    bore_r=bore_radius*in2cm

    str_zIA_top = openmc.ZPlane(z0= str_zIA_length/2.0)
    str_zIA_bot = openmc.ZPlane(z0=-str_zIA_length/2.0)
    str_zIA_sq_l = openmc.XPlane(x0=-str_zI_side/2.0)
    str_zIA_sq_r = openmc.XPlane(x0= str_zI_side/2.0)
    str_zIA_sq_f = openmc.YPlane(y0=-str_zI_side/2.0)
    str_zIA_sq_b = openmc.YPlane(y0= str_zI_side/2.0)

    str_zIA_dot_l = openmc.ZCylinder(x0=-str_zI_side/2.0,y0=-str_zI_side/2.0+str_zI_dot_r, r=str_zI_dot_r)
    str_zIA_dot_r = openmc.ZCylinder(x0= str_zI_side/2.0,y0= str_zI_side/2.0-str_zI_dot_r, r=str_zI_dot_r)
    str_zIA_dot_f = openmc.ZCylinder(y0=-str_zI_side/2.0,x0= str_zI_side/2.0-str_zI_dot_r, r=str_zI_dot_r)
    str_zIA_dot_b = openmc.ZCylinder(y0= str_zI_side/2.0,x0=-str_zI_side/2.0+str_zI_dot_r, r=str_zI_dot_r)

    str_zIA_bore = openmc.ZCylinder(r=bore_r)

    str_zIA_body = -str_zIA_top & +str_zIA_bot & ( ( +str_zIA_sq_l & -str_zIA_sq_r & +str_zIA_sq_f & -str_zIA_sq_b ) | -str_zIA_dot_l | - str_zIA_dot_r | -str_zIA_dot_f | -str_zIA_dot_b ) & +str_zIA_bore
    str_zIA_central = -str_zIA_top & +str_zIA_bot & -str_zIA_bore
    str_zIA_outside = +str_zIA_top | -str_zIA_bot | ( (-str_zIA_sq_l | +str_zIA_sq_r | -str_zIA_sq_f | +str_zIA_sq_b) & (+str_zIA_dot_l | +str_zIA_dot_r | +str_zIA_dot_f | +str_zIA_dot_b) )

    stringer_zIA_body = openmc.Cell(fill=self.md['graphite'], region=str_zIA_body)
    stringer_zIA_fuel_bore = openmc.Cell(fill=self.md['salt'], region=str_zIA_central)
    stringer_zIA_fuel_outside = openmc.Cell(fill=self.md['salt'], region=str_zIA_outside)
  
    return [stringer_zIA_body,stringer_zIA_fuel_bore,stringer_zIA_fuel_outside]

  def zoneIA_stringer(self,univ=True):
    cs=self.stringer_sq()
    if(univ):
      u=openmc.Universe(cells=cs)
      return u
    else:
      return cs

  def zoneIIA_stringer(self,univ=True):
    cs=self.stringer_sq(bore_radius=2.6/2.0)
    if(univ):
      u=openmc.Universe(cells=cs)
      return u
    else:
      return cs

  def zoneIA_IIB(self):
    #this function builds the geometry based on stringers
    pass

  def zoneIIB(self):
    #this function build the region with the radial undermoderated region
    pass
  
  def plot(self):
    pixels=800
    plotxy = openmc.Plot.from_geometry(self.geometry)
    plotxy.width=(100,100)
    plotxy.pixels=(pixels,pixels)
    plotxy.color_by = 'material'
    plotxy.origin=(0,0,0)
    #plotxy.colors = matcolors
    plotxy.to_ipython_image()


if __name__=='__main__':
  mm=msbr()
  mm.geom()
  mm.plot()
