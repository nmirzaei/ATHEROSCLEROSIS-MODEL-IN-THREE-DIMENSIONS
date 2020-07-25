
from dolfin import *
import numpy as np

#this function calculates the strain
def strain(u1,t,C_i,C_m,C_a,vtkfile1,vtkfile2,vtkfile3):

    #reading the original configs
    Domain = Mesh('Mesh.xml') #cut domain
    bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the Domain
    boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ###################################################################################
    V = VectorFunctionSpace(Domain,'P',1)
    u2 = Function(V)
    u1.set_allow_extrapolation(True)
    u2.interpolate(u1)
    ####################################################################################

    #Splitting the tensors into single entries
    [C_11_i,C_12_i,C_13_i,C_21_i,C_22_i,C_23_i,C_31_i,C_32_i,C_33_i] = C_i.split(True)
    [C_11_m,C_12_m,C_13_m,C_21_m,C_22_m,C_23_m,C_31_m,C_32_m,C_33_m] = C_m.split(True)
    [C_11_a,C_12_a,C_13_a,C_21_a,C_22_a,C_23_a,C_31_a,C_32_a,C_33_a] = C_a.split(True)
    ####################################################################################

    #Function spaces
    WW = FunctionSpace(Domain,'DG',0)
    VV = VectorFunctionSpace(Domain,'P',1)
    ###################################################################################
    u = Function(VV)
    u.interpolate(u2)
    ###################################################################################
    C_11 = project(C_11_i+C_11_m+C_11_a,WW)
    C_22 = project(C_22_i+C_22_m+C_22_a,WW)
    C_33 = project(C_33_i+C_33_m+C_33_a,WW)
    ###################################################################################

    #Finding the center of mass
    R = VectorFunctionSpace(Domain, "R", 0)
    V= VectorFunctionSpace(Domain, "P", 1)
    position = Function(V)
    position.assign(Expression(["x[0]", "x[1]","x[2]"], element=V.ufl_element()))
    c = TestFunction(R)
    bulk1 = assemble(Constant(1.0)*dx(domain=Domain))
    centroid = assemble(dot(c, position)*dx)
    f = centroid / bulk1
    f_np = f.get_local()
    ###################################################################################

    #moving mesh
    ALE.move(Domain,u2)
    ###################################################################################

    #finding the new center of mass
    new_V = VectorFunctionSpace(Domain, 'P', 1)
    new_R = VectorFunctionSpace(Domain, "R", 0)
    new_position = Function(new_V)
    new_position.assign(Expression(["x[0]", "x[1]","x[2]"], element=new_V.ufl_element()))
    new_c = TestFunction(new_R)
    new_bulk = assemble(Constant(1.0)*dx(domain=Domain))
    new_centroid = assemble(dot(new_c, new_position)*dx)
    new_f = new_centroid / new_bulk#
    new_f_np = new_f.get_local()
    dev = Constant(f_np - new_f_np)
    ###################################################################################


    #moving the deformed Domain back to the center of mass for the original Domain
    ALE.move(Domain , dev)
    ###################################################################################




    #renaming for animation purposes
    C_11.rename('C_11','C_11')
    C_22.rename('C_22','C_22')
    C_33.rename('C_33','C_33')
    ###################################################################################


    #saving into file
    vtkfile1 << (C_11,t)
    vtkfile2 << (C_22,t)
    vtkfile3 << (C_33,t)
    print "All Done!"
    ###################################################################################
