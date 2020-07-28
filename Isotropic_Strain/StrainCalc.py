from dolfin import *
import numpy as np
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
from strain import strain
from MaxComp import MaxComp
from Reader import Reader
import os
import csv
#The following lines take an already converted Gmsh mesh and import it
mesh = Mesh('Mesh.xml')
mesh1 = Mesh('Mesh1.xml')
######################################################################

# Facet functions
Volume = MeshFunction('size_t' , mesh , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh = MeshFunction('size_t', mesh , 'Mesh_facet_region.xml')  #saves the boundary info of the mesh
Volume1 = MeshFunction('size_t' , mesh , 'Mesh1_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh1 = MeshFunction('size_t', mesh , 'Mesh1_facet_region.xml')  #saves the boundary info of the mesh
######################################################################


# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
ffc_options = {'optimize': True}
######################################################################

# define function space
V = VectorFunctionSpace(mesh, 'P', 1)
V1 = VectorFunctionSpace(mesh1, 'P', 1)
S = FunctionSpace(mesh,'DG',0)
T = TensorFunctionSpace(mesh,'DG',0)
N = V.dim()
d = mesh.geometry().dim()
######################################################################


# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
######################################################################


# Define functions
u = Function(V)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
######################################################################



#defining spatial coordinates for changing to cylindrical coordinates
x = SpatialCoordinate(mesh)
######################################################################

#calculating the norm of a given vector
def Mag(u):
    MM = u[0]*u[0]+u[1]*u[1]+u[2]*u[2]
    M = MM**0.5
    return M
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
mu_i, nu_i, beta_i, eta_i, rho_i, phi_i = 27.9, 0.49, 170.88, 263.66, 0.51, 60.3*pi/180
######################################################################

# Elasticity parameter for the Hopzafel cylinder (Media)
mu_m, nu_m, beta_m, eta_m, rho_m, phi_m = 1.27, 0.49, 8.21, 21.60, 0.25, 20.61*pi/180
######################################################################

# Elasticity parameter for the Hopzafel cylinder (adventitia)
mu_a, nu_a, beta_a, eta_a, rho_a, phi_a = 7.56, 0.49, 85.03, 38.57, 0.55, 67*pi/180
######################################################################


#fibers for intitma
#in case of 3-D you should add x[2] for z
x, y, z= sym.symbols('x[0], x[1], x[2]')
theta = sym.atan2(y,x)
t = theta + pi
r_1 =  0.7512+0.0124*sym.cos(t)+0.0414*sym.cos(2*t)+0.0120*sym.cos(3*t)+0.0111*sym.cos(4*t)-0.0061*sym.cos(5*t)-0.0008*sym.cos(6*t)-0.00001*sym.cos(7*t)+0.0033*sym.cos(8*t)+ 0.0217*sym.sin(t)-0.0197*sym.sin(2*t)-0.0027*sym.sin(3*t)-0.0084*sym.sin(4*t)-0.0135*sym.sin(5*t)-0.0034*sym.sin(6*t)-0.0001*sym.sin(7*t)-0.0048*sym.sin(8*t)
r_2 = 1.0284+0.0362*sym.cos(t)-0.0557*sym.cos(2*t)+0.0190*sym.cos(3*t)-0.0185*sym.cos(4*t)+0.0038*sym.cos(5*t)+0.0045*sym.cos(6*t)-0.0002*sym.cos(7*t)+0.0006*sym.cos(8*t)-0.1984*sym.sin(t)-0.0173*sym.sin(2*t)-0.0008*sym.sin(3*t)-0.0133*sym.sin(4*t)-0.0026*sym.sin(5*t)-0.0066*sym.sin(6*t)-0.0029*sym.sin(7*t)+0.0064*sym.sin(8*t)
f_1_x = (r_1)*sym.cos(t)
f_1_y = (r_1)*sym.sin(t)
f_2_x = (r_2)*sym.cos(t)
f_2_y = (r_2)*sym.sin(t)
w_A_x = sym.diff(f_1_x,theta)
w_A_y = sym.diff(f_1_y,theta)
w_C_x = sym.diff(f_2_x,theta)
w_C_y = sym.diff(f_2_y,theta)
AE = ((x-f_1_x)**2+(y-f_1_y)**2)**0.5
EC = ((x-f_2_x)**2+(y-f_2_y)**2)**0.5
u_s = (AE*w_A_x+EC*w_C_x)/(AE+EC)
v_s = (AE*w_A_y+EC*w_C_y)/(AE+EC)
w_s = (sym.tan(phi_i)*(u_s**2+v_s**2)**0.5)
M_0 = 1/((u_s**2+v_s**2)*(1+sym.tan(phi_i)**2))**0.5
U_s = M_0*u_s
V_s = M_0*v_s
W_s = M_0*w_s
THETA_s = sym.atan2(V_s,U_s)
SinT_0 = sym.sin(THETA_s)
nSinT_0= -sym.sin(THETA_s)
CosT_0 = sym.cos(THETA_s)
u_s_code= sym.printing.ccode(U_s)
v_s_code= sym.printing.ccode(V_s)
w_s_code= sym.printing.ccode(W_s)
SinT_0_code=sym.printing.ccode(SinT_0)
CosT_0_code = sym.printing.ccode(CosT_0)
b_i=Expression((u_s_code,v_s_code,w_s_code), degree=0)
sint_i = Expression(SinT_0_code,degree=0)
cost_i = Expression(CosT_0_code ,degree=0)
###############################################################

#fibers for hopzafel (media)
r_3 = 1.2584+0.0362*sym.cos(t)-0.0557*sym.cos(2*t)+0.0190*sym.cos(3*t)-0.0185*sym.cos(4*t)+0.0038*sym.cos(5*t)+0.0045*sym.cos(6*t)-0.0002*sym.cos(7*t)+0.0006*sym.cos(8*t)-0.1984*sym.sin(t)-0.0173*sym.sin(2*t)-0.0008*sym.sin(3*t)-0.0133*sym.sin(4*t)-0.0026*sym.sin(5*t)-0.0066*sym.sin(6*t)-0.0029*sym.sin(7*t)+0.0064*sym.sin(8*t)
f_2_x = (r_2)*sym.cos(t)
f_2_y = (r_2)*sym.sin(t)
f_3_x = (r_3)*sym.cos(t)
f_3_y = (r_3)*sym.sin(t)
v_A_x = sym.diff(f_2_x,theta)
v_A_y = sym.diff(f_2_y,theta)
v_C_x = sym.diff(f_3_x,theta)
v_C_y = sym.diff(f_3_y,theta)
AB = ((x-f_2_x)**2+(y-f_2_y)**2)**0.5
BC = ((x-f_3_x)**2+(y-f_3_y)**2)**0.5
u_q = (AB*v_A_x+BC*v_C_x)/(AB+BC)
v_q = (AB*v_A_y+BC*v_C_y)/(AB+BC)
w_q = (sym.tan(phi_m)*(u_q**2+v_q**2)**0.5)
M_1 = 1/((u_q**2+v_q**2)*(1+sym.tan(phi_m)**2))**0.5
U_q = M_1*u_q
V_q = M_1*v_q
W_q = M_1*w_q
THETA_q = sym.atan2(V_q,U_q)
SinT_1 = sym.sin(THETA_q)
nSinT_1= -sym.sin(THETA_q)
CosT_1 = sym.cos(THETA_q)
u_q_code= sym.printing.ccode(U_q)
v_q_code= sym.printing.ccode(V_q)
w_q_code= sym.printing.ccode(W_q)
SinT_1_code=sym.printing.ccode(SinT_1)
CosT_1_code = sym.printing.ccode(CosT_1)
b_m=Expression((u_q_code,v_q_code,w_q_code), degree=0)
sint_m = Expression(SinT_1_code,degree=0)
cost_m = Expression(CosT_1_code ,degree=0)
###############################################################

#fibers for hopzafel (adventitia)
r_4 = 1.3584+0.0362*sym.cos(t)-0.0557*sym.cos(2*t)+0.0190*sym.cos(3*t)-0.0185*sym.cos(4*t)+0.0038*sym.cos(5*t)+0.0045*sym.cos(6*t)-0.0002*sym.cos(7*t)+0.0006*sym.cos(8*t)-0.1984*sym.sin(t)-0.0173*sym.sin(2*t)-0.0008*sym.sin(3*t)-0.0133*sym.sin(4*t)-0.0026*sym.sin(5*t)-0.0066*sym.sin(6*t)-0.0029*sym.sin(7*t)+0.0064*sym.sin(8*t)
f_4_x = (r_4)*sym.cos(t)
f_4_y = (r_4)*sym.sin(t)
u_A_x = sym.diff(f_3_x,theta)
u_A_y = sym.diff(f_3_y,theta)
u_C_x = sym.diff(f_4_x,theta)
u_C_y = sym.diff(f_4_y,theta)
AD = ((x-f_3_x)**2+(y-f_3_y)**2)**0.5
DC = ((x-f_4_x)**2+(y-f_4_y)**2)**0.5
u_p = (AD*u_A_x+DC*u_C_x)/(AD+DC)
v_p = (AD*u_A_y+DC*u_C_y)/(AD+DC)
w_p = (sym.tan(phi_a)*(u_p**2+v_p**2)**0.5)
M_2 = 1/((u_p**2+v_p**2)*(1+sym.tan(phi_a)**2))**0.5
U_p = M_2*u_p
V_p = M_2*v_p
W_p = M_2*w_p
THETA_p = sym.atan2(V_p,U_p)
SinT_2 = sym.sin(THETA_p)
nSinT_2= -sym.sin(THETA_p)
CosT_2 = sym.cos(THETA_p)
u_p_code= sym.printing.ccode(U_p)
v_p_code= sym.printing.ccode(V_p)
w_p_code= sym.printing.ccode(W_p)
SinT_2_code=sym.printing.ccode(SinT_2)
CosT_2_code = sym.printing.ccode(CosT_2)
b_a = Expression((u_p_code,v_p_code,w_p_code), degree=0)
sint_a = Expression(SinT_2_code,degree=0)
cost_a = Expression(CosT_2_code ,degree=0)
###############################################################


#growth for the media
ginv_m = Identity(3)
G_m = Expression('1', degree=0)
######################################################################

#growth for the adventitia
ginv_a = Identity(3)
G_a = Expression('1', degree=0)
######################################################################


# blood pressure direction
n = FacetNormal(mesh)
########################################################################################################################

#User input for the desirable stress
M4 = 0
Choose = raw_input('color map?(Yes=1, No =0)')
MC = int(Choose)
if MC == 1:
    response2 = raw_input('At what step do you want the color map:')
    M2 = int(response2)
response4 = raw_input('Do you also want to calculate the maximum strains?(Yes=1, No =0)')
M4 = int(response4)

response = raw_input('How many t do you want to calculate the stresses?')
M = int(response)

response1 = raw_input('Skip?:')
M1 = int(response1)
######################################################################



######################################################################
vtkfile1 = File('home directory/Residual/C_11.pvd')
vtkfile2 = File('home directory/Residual/C_22.pvd')
vtkfile3 = File('home directory/Residual/C_33.pvd')
#######################################################################

#looping prameters
counter=1
epsilon, sigma, t = 0.0, 16.0, 0.0
deps, dsig, dt= 0.001, 2.0, 0.001
epsmax, sigmax = 0.001, 16.0
j=0
k=0
ID=0
max11_i = []
max22_i = []
max33_i = []
max11_m = []
max22_m = []
max33_m = []
max11_a = []
max22_a = []
max33_a = []
#######################################################################

while counter<=M:
 #reading isotropic displacements from the directory they are saved in.
 if counter%M1==0:
    u = Function(V,'directory/displacement/u%d.xml' %(counter))
    u1 = Function(V,'directory/displacement/u%d.xml' %(counter))
 #######################################################################
    #Reading growth parameters by the Reader function and calculating the isotropic growth parameter
    if counter!=0:
        a1, b1, c1 = Reader(counter)
        alpha = (-1+((1+a1*epsilon)*(1+b1*epsilon)*(1+c1*epsilon))**(0.3333))/epsilon
    else:
        alpha = Constant(0)
    #######################################################################


    #Constructing the growth tensor from the growth parameters above
    EXP = Expression('eps*exp(-a*(x[2]-0.5)*(x[2]-0.5)-a*(x[1]+0.8)*(x[1]+0.8)-a*(x[0])*(x[0]))',degree=0,a=4.25,eps=epsilon)
    g_i = as_tensor([[(1+alpha*EXP)*(cost_i*cost_i)+(1+alpha*EXP)*(sint_i*sint_i),(1+alpha*EXP)*(cost_i*sint_i)-(1+alpha*EXP)*(cost_i*sint_i),0],[(1+alpha*EXP)*(cost_i*sint_i)-(1+alpha*EXP)*(cost_i*sint_i),(1+alpha*EXP)*(sint_i*sint_i)+(1+alpha*EXP)*(cost_i*cost_i),0],[0,0,1+alpha*EXP]])
    Transform = as_tensor([[cost_i,-1*sint_i,0],[sint_i,cost_i,0],[0,0,1]])
    G_i = (1+alpha*EXP)*(1+alpha*EXP)*(1+alpha*EXP)
    ginv_i = as_tensor([[(1/(1+alpha*EXP))*(cost_i*cost_i)+(1/(1+alpha*EXP))*(sint_i*sint_i),(1/(1+alpha*EXP))*(cost_i*sint_i)-(1/(1+alpha*EXP))*(cost_i*sint_i),0],[(1/(1+alpha*EXP))*(sint_i*cost_i)-(1/(1+alpha*EXP))*(sint_i*cost_i),(1/(1+alpha*EXP))*(sint_i*sint_i)+(1/(1+alpha*EXP))*(cost_i*cost_i),0],[0,0,1/(1+alpha*EXP)]])
    #######################################################################

    # Kinematics for initma
    #    U1.interpolate(u)
    II = Identity(3)            # Identity tensor
    F = II+grad(u)
    #project the fiber direction onto xy and then deform it
    BBi = as_vector([b_i[0],b_i[1],0])
    Bi = F*BBi
    bi = (1/Mag(Bi))*Bi
    Trans_i = as_tensor([[bi[0],-bi[1],0],[bi[1],bi[0],0],[0,0,1]])
    strain_i= F*ginv_i
    cc_i = (Trans_i.T)*strain_i*(Trans_i)
    c_i = project(cc_i,T)
    C_i = Function(T); C_i.vector()[:]=c_i.vector()
    C_i = extract_values(C_i,Volume,12,T)
    #######################################################################
    # Kinematics for media
    #project the fiber direction onto xy and then deform it
    BBm = as_vector([b_m[0],b_m[1],0])
    Bm = F*BBm
    bm = (1/Mag(Bm))*Bm
    Trans_m = as_tensor([[bm[0],-bm[1],0],[bm[1],bm[0],0],[0,0,1]])
    strain_m= F
    cc_m = (Trans_m.T)*strain_m*(Trans_m)
    c_m= project(cc_m,T)
    C_m = Function(T); C_m.vector()[:]=c_m.vector()
    C_m = extract_values(C_m,Volume,13,T)
    #######################################################################
    # Kinematics for adventitia
    #project the fiber direction onto xy and then deform it
    BBa = as_vector([b_a[0],b_a[1],0])
    Ba = F*BBa
    ba = (1/Mag(Ba))*Ba
    Trans_a = as_tensor([[ba[0],-ba[1],0],[ba[1],ba[0],0],[0,0,1]])
    strain_a= F
    cc_a = (Trans_a.T)*strain_a*(Trans_a)
    c_a = project(cc_a,T)
    C_a = Function(T); C_a.vector()[:]=c_a.vector()
    C_a = extract_values(C_a,Volume,14,T)
    #######################################################################
    #Calculating strain
    if MC==1 and counter%M2==0:
        strain(u1,t,C_i,C_m,C_a,vtkfile1,vtkfile2,vtkfile3)
    #######################################################################
    #calculating strain maxima
    if M4 == 1:
        i1, i2, i3, m1, m2, m3, a1, a2, a3 = MaxComp(u,t,C_i,C_m,C_a)
        max11_i.append(i1)
        max22_i.append(i2)
        max33_i.append(i3)
        max11_m.append(m1)
        max22_m.append(m2)
        max33_m.append(m3)
        max11_a.append(a1)
        max22_a.append(a2)
        max33_a.append(a3)
    #######################################################################


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
 print(epsilon)#,flush=True)
 print(t)#,flush=True)
 print(sigma)#,flush=True)
######################################################################

#updating
 k+=1
######################################################################

#Saving strain maxima in csv files
if M4==1:
    import csv

    with open('max11_i.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max11_i)

    with open('max22_i.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max22_i)


    with open('max33_i.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max33_i)


    with open('max11_m.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max11_m)


    with open('max22_m.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max22_m)


    with open('max33_m.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max33_m)


    with open('max11_a.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max11_a)


    with open('max22_a.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max22_a)


    with open('max33_a.csv', 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(max33_a)
else:
    print("No max strain was calculated!")
#######################################################################
