from dolfin import *
import numpy as np
import os

#This function creates our animation
def Reader(counter):

    #Home PC or HPC home directory
    os.chdir("/home/1925/Fenics/Evolution_Cylindrical")
    ##################################################################################################################

    #create the mesh and its volume and boundary
    Domain=Mesh('Mesh.xml')
    Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ##################################################################################################################

    #Creating a function space
    S = FunctionSpace(Domain,'R',0)
    ##################################################################################################################

    #Reading growth parameters computed from before. Depending where they are saved the home directory should change
    a1 = Function(S,'home directory/alpha(%d).xml' %(counter))
    b1= Function(S,'home directory/beta(%d).xml' %(counter))
    c1 = Function(S,'home directory/gamma(%d).xml' %(counter))
    ##################################################################################################################

    #After extracting parameter from the home directory we can send the directory back to the scratch directory for saving large outputs.
    os.chdir("Scratch directory")
    ##################################################################################################################
        
    return a1.vector().get_local()[0], b1.vector().get_local()[0], c1.vector().get_local()[0]
