from dolfin import *
import numpy as np
import os

#This function creates our animation
def Reader(counter):
    os.chdir("/home/1925/Fenics/Opening_Angle_Cylindrical")
    #create the mesh and its volume and boundary
    Domain=Mesh('Mesh.xml')
    Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ##################################################################################################################

    #Function space
    S = FunctionSpace(Domain,'R',0)
    ##################################################################################################################

    #Reading the parameters from the home directory where it was saved
    a1 = Function(S,'home directory/alpha(%d).xml' %(counter))
    b1= Function(S,'home directory/beta(%d).xml' %(counter))
    c1 = Function(S,'home directory/gamma(%d).xml' %(counter))
    ##################################################################################################################

    #changing the directory back to the scratch memory in case of large outputs
    os.chdir("/lustre/scratch/nmirzaei")
    ##################################################################################################################
    
    return a1.vector().get_local()[0], b1.vector().get_local()[0], c1.vector().get_local()[0]
