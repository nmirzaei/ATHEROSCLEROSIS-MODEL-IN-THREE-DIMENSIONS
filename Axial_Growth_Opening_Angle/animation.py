from dolfin import *
import os
#This function creates an animation for us for the deformation
def animation(kk,u,t,vtkfile1,vtkfile2,vtkfile3,vtkfile4,vtkfile5):

    #Home Pc or HPC directory
    os.chdir("Home directory")
    #########################################################################################################################

    #Importing the mesh and its volume and facet
    Domain=Mesh('Mesh.xml')
    Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    #########################################################################################################################

    #Creating the required function spaces
    R = VectorFunctionSpace(Domain, "R", 0)
    V= VectorFunctionSpace(Domain, "P", 1)
    VV = VectorFunctionSpace(Domain, "P", 1)
    #########################################################################################################################

    #In case you have large outputs send them to a scratch directory on HPC
    os.chdir("/lustre/scratch/nmirzaei")
    #########################################################################################################################

    #Finding the center of mass for the reference domain
    position = Function(V)
    position.assign(Expression(["x[0]", "x[1]", "x[2]"], element=V.ufl_element()))
    c = TestFunction(R)
    volume = assemble(Constant(1.0)*dx(domain=Domain))
    centroid = assemble(dot(c, position)*dx)
    f = centroid / volume
    f_np = f.get_local()
    #########################################################################################################################

    #Moving the domain
    if kk>0:
        ALE.move(Domain,u)
    #########################################################################################################################

    #Finding the center of mass for the deformed mesh
    new_V = VectorFunctionSpace(Domain, 'P', 1)
    new_R = VectorFunctionSpace(Domain, "R", 0)
    new_position = Function(new_V)
    new_position.assign(Expression(["x[0]", "x[1]", "x[2]"], element=new_V.ufl_element()))
    new_c = TestFunction(new_R)
    new_volume = assemble(Constant(1.0)*dx(domain=Domain))
    new_centroid = assemble(dot(new_c, new_position)*dx)
    new_f = new_centroid / new_volume
    new_f_np = new_f.get_local() # numpy array
    #########################################################################################################################

    #deviation of the deformed center of mass from the original one
    dev = Constant(f_np - new_f_np)
    #########################################################################################################################

    #moving the deformed mesh back to the center of mass for the original mesh
    if kk>0:
        ALE.move(Domain , dev)
    #########################################################################################################################

    #Extracting the intima
    Sub_Mesh=SubMesh(Domain,Bulk,12)
    surface_marker = MeshFunction("size_t", Sub_Mesh, Sub_Mesh.topology().dim() - 1, 0)
    ncells = MeshFunction("size_t", Sub_Mesh, Sub_Mesh.topology().dim())
    volume_marker = MeshFunction("size_t", Sub_Mesh, Sub_Mesh.topology().dim())
    vmap = Sub_Mesh.data().array('parent_vertex_indices', 0)
    cmap = Sub_Mesh.data().array('parent_cell_indices', Sub_Mesh.topology().dim())
    #########################################################################################################################

    #########################################################################################################################
    n = 0
    for c in cells(Sub_Mesh):
         parent_cell = Cell(Domain, cmap[c.index()])
         volume_marker.array()[c.index()] = Bulk.array()[parent_cell.index()]
         for f in facets(parent_cell):
              for g in facets(c):
                  g_vertices = vmap[g.entities(0)]
                  if set(f.entities(0)) == set(g_vertices):
                       surface_marker.array()[g.index()] = Boundary.array()[f.index()]
              n=n+1
    #########################################################################################################################

    #Saving frames
    vtkfile1<<(Bulk,t)
    vtkfile2<<(Boundary,t)
    vtkfile3<<(Domain,t)
    vtkfile4<<(volume_marker,t)
    vtkfile5<<(surface_marker,t)
    #########################################################################################################################
