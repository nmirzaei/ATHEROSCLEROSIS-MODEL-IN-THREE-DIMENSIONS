from dolfin import *
import numpy as np
import os
import csv
#The following lines take an already converted Gmsh mesh and import it
mesh = Mesh('Mesh.xml')
######################################################################

# Facet functions
Volume = MeshFunction('size_t' , mesh , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh = MeshFunction('size_t', mesh , 'Mesh_facet_region.xml')  #saves the boundary info of the mesh
######################################################################


# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
ffc_options = {'optimize': True}
######################################################################

# define function space
V = VectorFunctionSpace(mesh, 'P', 1)
S = FunctionSpace(mesh,'R',0)
T1 = TensorFunctionSpace(mesh,'DG',0,shape=(3,3))
N = V.dim()
d = mesh.geometry().dim()
######################################################################


# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
######################################################################


# Define functions
U1 = Function(V)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
######################################################################



#defining spatial coordinates for changing to cylindrical coordinates
x = SpatialCoordinate(mesh)
######################################################################



#Defining a function that sets all the irrelevant layer values to zero.
def extract_values(u,cell_function,subdomain_id, V):
    dofmap=V.dofmap()
    mesh = V.mesh()
    for cell in cells(mesh):
        if cell_function[cell.index()]!=subdomain_id:
            dofs = dofmap.cell_dofs(cell.index())
            for dof in dofs:
                u.vector()[dof]=0.0
    return u
######################################################################


######################################################################
# Elasticity parameter for the Hopzafel cylinder (intima)
mu_i, nu_i, beta_i, eta_i, rho_i, phi_i = 27.9, 0.49, 170.88, 0, 0.51, 60.3*pi/180

# Elasticity parameter for the Hopzafel cylinder (Media)
mu_m, nu_m, beta_m, eta_m, rho_m, phi_m = 1.27, 0.49, 8.21, 0, 0.25, 20.61*pi/180

# Elasticity parameter for the Hopzafel cylinder (adventitia)
mu_a, nu_a, beta_a, eta_a, rho_a, phi_a = 7.56, 0.49, 85.03, 0, 0.55, 67*pi/180
######################################################################


#Collagen fibers direction in cylindrical coordinates  (before it was (0,cos,sin) I changed it to (0,sin,cos) which I believe is the right one )
b_i = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_i)
b_m = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_m)
b_a = Expression(('0','sin(phi)','cos(phi)'),degree=0,phi=phi_a)
######################################################################


#growth for the media
ginv_m = Identity(3)
G_m = Expression('1', degree=0)
######################################################################

#growth for the adventitia
ginv_a = Identity(3)
G_a = Expression('1', degree=0)
######################################################################


#facet Normals for blood pressure direction
n = Expression(('-1','0','0'),degree=0)
######################################################################


#User input for the desirable stress
response = raw_input('How many growth parameters do you want to extract?')
M = int(response)
response1 = raw_input('Calculate growth parameters every other -- term:')
Skip = int(response1)
######################################################################

#radius variable
radius = Expression("x[0]",degree=0)
######################################################################

#looping prameters
counter=0
epsilon, sigma, t = 0.0, 12.0, 0.0
deps, dsig, dt= 0.001, 2.0, 0.001
epsmax, sigmax = 0.001, 0.0
j=0
k=0
ID=0
Alpha = []
Beta=[]
Gamma=[]
######################################################################


while counter<=M:


 #Reading the desirable displacement
 if counter%Skip==0:
    alpha = Function(S,'Growth Directory/alpha(%d).xml' %(counter))
    beta = Function(S,'Growth Directory/beta(%d).xml' %(counter))
    gamma = Function(S,'Growth Directory/gamma(%d).xml' %(counter))
 ######################################################################



    Alpha.append(alpha.vector().get_local()[0])
    Beta.append(beta.vector().get_local()[0])
    Gamma.append(gamma.vector().get_local()[0])
######################################################################

#making sure we start saving results after the pressure is imposed (if there is any pressure). Also making sure we save results every 10 steps
 if sigma<=sigmax:
    sigma+=dsig
 if sigma>sigmax:
    epsilon+=deps
    t+=dt
    counter+=1
    ID+=1
 if ID!=0:
    sigma-=dsig
    ID=0
######################################################################

#printing the current loop variables
 print(epsilon,flush=True)
 print(t,flush=True)
 print(sigma,flush=True)
######################################################################

 k+=1
######################################################################



#saving in csv for matlab postprocessing
import csv

with open('alpha.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(Alpha)
with open('beta.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(Beta)
with open('gamma.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(Gamma)
######################################################################
