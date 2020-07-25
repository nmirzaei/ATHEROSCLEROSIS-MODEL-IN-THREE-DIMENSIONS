from dolfin import *
import numpy as np
import os

#This function creates our animation
def Reader(counter):
    Domain=Mesh('Mesh.xml')  #cut domain
    Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
    Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ##################################################################################################################
    #real function space for the parameters
    S = FunctionSpace(Domain,'R',0)
    ##################################################################################################################

    #resing the parameters from where they were saved
    a1 = Function(S,'directory/alpha(%d).xml' %(counter))
    b1= Function(S,'directory/beta(%d).xml' %(counter))
    c1 = Function(S,'directory/gamma(%d).xml' %(counter))
    ##################################################################################################################

    return a1.vector().get_local()[0], b1.vector().get_local()[0], c1.vector().get_local()[0]
