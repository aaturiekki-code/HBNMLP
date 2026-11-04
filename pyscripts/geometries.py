from __future__ import print_function,division
import numpy as np
from pyqula import geometry
from pyqula import sculpt
from pyqula import supercell
from pyqula import specialgeometry

def twisted_multilayer(m0=3,rotate=True,shift=None,
  sublattice=True,r=1,rot=[1,1,0,0],g=None,dz=None, rotate_center="atom", ab=0, mirror = 0):
  natoms = int(12*m0**2+12*m0*r+4*r**2)
  """Return the geometry for twisted multilayer graphene"""
  if g is None: g = geometry.honeycomb_lattice() # default is honeycomb
  if dz is None: # width + an (arbitrary) distance of 3
      dz = np.max(g.r[:,2]) - np.min(g.r[:,2]) + 3. # distance
  else: g = geometry.get_geometry(g) # get the geometry
#  g.has_sublattice = False # no sublattice
  from pyqula import supercell as supercelltk
# put the right angle between lattice vectors
  g = supercelltk.target_angle_volume(g,angle=1./3.,volume=1)
# old version
#  g = geometry.non_orthogonal_supercell(g,m=[[-1,0,0],[0,1,0],[0,0,1]])
#  print(g.a1,g.a2,g.a1.dot(g.a2)) ; exit()
  g0 = g.copy() # copy geometry
  theta = np.arccos((3.*m0**2+3*m0*r+r**2/2.)/(3.*m0**2+3*m0*r+r**2))
  if shift is None:
      shift = [[0.,0.] for r in rot] # initialize
  print("Theta",theta*180.0/np.pi)
  nsuper = [[m0,m0+r,0],[-m0-r,2*m0+r,0],[0,0,1]]
  g = geometry.non_orthogonal_supercell(g,m=nsuper,
           reducef=lambda x: 3*np.sqrt(x))

  g1 = g.copy().get_supercell(2)
  g2 = g.copy().get_supercell(2)

  g1.set_origin()
  g2.set_origin()

  dr = np.array([0, 0, 0])
  for i in range(0, natoms):
     if 1.01>np.linalg.norm((g1.r[0]-g1.r[i]))>0.99:
        dr = g1.r[0]-g1.r[i]
        print(dr)
        break
  
  if rotate_center == "hex":
     if ab == 0:
        g1.r = g1.r-dr
        g2.r = g2.r-dr
     elif ab == 1:
        g2.r = g2.r-dr
  elif rotate_center == "atom":
     if ab == 0:
        ()
     elif ab == 1:
        g1.r = g1.r+dr
  if mirror == 1:
     g2.sublattice = g2.sublattice*-1
    
  g1.r2xyz()
  g2.r2xyz()
  g1.z = g1.z - 1.18
  g2.z = g2.z + 1.18
  g1.xyz2r()
  g2.xyz2r()
  g2 = g2.rotate(theta*180.0/np.pi)

  g.r = np.concatenate([g1.r, g2.r]).copy() 
  g.r2xyz()
  if g.has_sublattice:
      g.sublattice = np.concatenate([g1.sublattice, g2.sublattice]).copy() 

  inuc = sculpt.sites_in_unit_cell(g.r,g.a1,g.a2,g.a3,dim=2)
  g.r = g.r[inuc==1.]
  g.sublattice = g.sublattice[inuc==1.]
  g.r2xyz() # update
  g.real2fractional() # update fractional coordinates 
  g = sculpt.rotate_a2b(g,g.a1,np.array([1.,0.,0.])) # rotate
  g.get_fractional() # get fractional coordinates 
  return g

def bilayer(l, mirror):
   gs=geometry.honeycomb_lattice()
   gs.center_in_atom(0)
   gs.xyz2r()
   gs = multilayer(gs, l=l)
   if mirror == 1:
         gs.sublattice = gs.sublattice*-1
         gs.r2xyz()
   return gs

def multilayer(g,l=[0],dr=[1.,0.,0.],dz=[0.,0.,1.]):
    g.center_in_atom(0)
    g.xyz2r()
    """Return a multilayer geometry"""
    if type(l)==str: # if it is a string, make the replacement
        l = l.replace("A","0")
        l = l.replace("B","1")
        l = l.replace("C","2")
        l = [int(il) for il in l] # transform to numbers
    dr = np.array(dr) # to array
    dz = np.array(dz) # to array
    ss = [] # list for the sublattice
    rs = [] # list for the positions
    ii = 0 # start
    for il in l: # loop over layers
        if g.has_sublattice:
            for (ri,si) in zip(g.r,g.sublattice): # loop over positions
                rs.append(ri + dr*il + (ii-1/2)*dz) # add position
                ss.append(si) # add sublattice
        else:
            for ri in g.r: # loop over positions
                rs.append(ri + dr*il + ii*dz) # add position
        ii += 1
    go = g.copy()
    go.r = np.array(rs)
    if g.has_sublattice: go.sublattice = np.array(ss)
    go.r2xyz()
    return go

def rotate_layer(g,theta,dr=[0.,0.,0.]):
    """Rotate one of the layers, retainign the same unit cell"""
    g1 = g.copy() # copy geometry
    g1 = g1.rotate(theta*180/np.pi)
    g1s = g1.supercell(2) # supercell
#    g1.a1 = g.a1.copy()
#    g1.a2 = g.a2.copy()
    inuc = sculpt.sites_in_unit_cell(g1s.r,g.a1,g.a2,g.a3,dim=2)
    g1.r = g1s.r[inuc==1.]
    if g1.has_sublattice:
        g1.sublattice = g1s.sublattice[inuc==1.]
#    g1.r = sculpt.retain_unit_cell(g1s.r,g.a1,g.a2,g.a3,dim=2) # get positions
    g1.r2xyz() # update
    g1.x += dr[0] # shift upper layer
    g1.y += dr[1] # shift upper layer
    g1.z += dr[2] # shift upper layer
    g1.xyz2r() # update
    return g1 # return the geometry

def change_sublattice(g):
   g.sublattice

def multilayer(g,l=[0],dr=[1.,0.,0.],dz=[0.,0.,1.]):
    """Return a multilayer geometry"""
    if type(l)==str: # if it is a string, make the replacement
        l = l.replace("A","0")
        l = l.replace("B","1")
        l = l.replace("C","2")
        l = [int(il) for il in l] # transform to numbers
    dr = np.array(dr) # to array
    dz = np.array(dz) # to array
    ss = [] # list for the sublattice
    rs = [] # list for the positions
    ii = 0 # start
    for il in l: # loop over layers
        if g.has_sublattice:
            for (ri,si) in zip(g.r,g.sublattice): # loop over positions
                rs.append(ri + dr*il + (ii-1/2)*dz) # add position
                ss.append(si) # add sublattice
        else:
            for ri in g.r: # loop over positions
                rs.append(ri + dr*il + ii*dz) # add position
        ii += 1
    go = g.copy()
    go.r = np.array(rs)
    if g.has_sublattice: go.sublattice = np.array(ss)
    go.r2xyz()
    return go
   

def supercell2d(g,n1=1,n2=1):
  """ Creates a supercell for a 2d system"""
  go = g.copy() # copy geometry
  if True: # brute force
      n = len(g.r)*n1*n2 # total number of positions
      rs = np.zeros((n,3)) # storage
      rs = supercell.supercell2d_jit(g.r.real,n1,n2,g.a1.real,g.a2.real,rs.real) # get the replicas
  else: # jit (not as fast for some reason)
      rs = replicate3d(g.r,g.a1,g.a2,np.array([0.,0.,1.]),n1,n2,1)
  go.r = np.array(rs) # store
  go.r2xyz()
  go.a1 = go.a1*n1
  go.a2 = go.a2*n2
  # shift to zero
  go.center()
  if g.has_sublattice: # supercell sublattice
    go.sublattice = np.concatenate([g.sublattice for i in range(n1*n2)])
  if g.atoms_have_names: # supercell sublattice
    go.atoms_names = g.atoms_names*n1*n2
  go.get_fractional() # get fractional coordinates
  supercell.record_diagonal_supercell(go,len(g.r),n1,n2,1) # unfolding bookkeeping
  return go
