from dolfin import *
import numpy as np

#this function calculates the maximum strqin
def MaxComp(U,t,F_i,F_m,F_a):

    #reading the original configs
    Domain = Mesh('Mesh.xml')  #cut domain
    bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the Domain
    boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
    ###################################################################################

    #allowing extrapolation in case of mesh refinement
    U.set_allow_extrapolation(True)
    ###################################################################################

    #Splitting the tensors into single entries
    [F_11_i,F_12_i,F_13_i,F_21_i,F_22_i,F_23_i,F_31_i,F_32_i,F_33_i] = F_i.split(True)
    [F_11_m,F_12_m,F_13_m,F_21_m,F_22_m,F_23_m,F_31_m,F_32_m,F_33_m] = F_m.split(True)
    [F_11_a,F_12_a,F_13_a,F_21_a,F_22_a,F_23_a,F_31_a,F_32_a,F_33_a] = F_a.split(True)
    ####################################################################################

    #defining function spaces
    VV = VectorFunctionSpace(Domain,'P',1)
    WW = FunctionSpace(Domain,'DG',0)
    TT = TensorFunctionSpace(Domain,'DG',0)
    u = Function(VV)
    u.interpolate(U)
    ###################################################################################

    #calculate the maximum
    def MaxComp(f, subdomains, subd_id):
        '''Minimum of f over subdomains cells marked with subd_id'''
        V = f.function_space()

        dm = V.dofmap()

        subd_dofs = np.unique(np.hstack(
                [dm.cell_dofs(c.index()) for c in SubsetIterator(subdomains, subd_id)]))

        Array= f.vector().get_local()[subd_dofs]
        Array[Array<0]=0
        return np.max(Array)
    ###################################################################################


    #calculating the center of mass for the reference domain
    R = VectorFunctionSpace(Domain, "R", 0)
    V= VectorFunctionSpace(Domain, "P", 1)
    position = Function(V)
    position.assign(Expression(["x[0]", "x[1]","x[2]"], element=V.ufl_element()))
    c = TestFunction(R)
    bulk0 = assemble(Constant(1.0)*dx(domain=Domain))
    centroid = assemble(dot(c, position)*dx)
    f = centroid / bulk0
    f_np = f.get_local()
    ###################################################################################

    #Moving the mesh
    ALE.move(Domain,u)
    ###################################################################################

    #finding the new center of mass
    new_V = VectorFunctionSpace(Domain, 'P', 1)
    new_R = VectorFunctionSpace(Domain, "R", 0)
    new_position = Function(new_V)
    new_position.assign(Expression(["x[0]", "x[1]","x[2]"], element=new_V.ufl_element()))
    new_c = TestFunction(new_R)
    new_bulk = assemble(Constant(1.0)*dx(domain=Domain))
    new_centroid = assemble(dot(new_c, new_position)*dx)
    new_f = new_centroid / new_bulk
    new_f_np = new_f.get_local()
    dev = Constant(f_np - new_f_np)
    ###################################################################################


    #moving the deformed Domain back to the center of mass for the original Domain
    ALE.move(Domain , dev)
    ###################################################################################



    #Send functions to MaxComp to compute their maxima on the respective subdomain
    Max_11_i = MaxComp(F_11_i,bulk,12)

    Max_22_i = MaxComp(F_22_i,bulk,12)

    Max_33_i = MaxComp(F_33_i,bulk,12)

    Max_11_m = MaxComp(F_11_m,bulk,13)

    Max_22_m = MaxComp(F_22_m,bulk,13)

    Max_33_m = MaxComp(F_33_m,bulk,13)

    Max_11_a = MaxComp(F_11_a,bulk,14)

    Max_22_a = MaxComp(F_22_a,bulk,14)

    Max_33_a = MaxComp(F_33_a,bulk,14)

    print "All Done!"

    return Max_11_i, Max_22_i, Max_33_i,Max_11_m,  Max_22_m,  Max_33_m ,Max_11_a, Max_22_a, Max_33_a
    ###################################################################################
