from dolfin import *
import numpy as np
import sympy as sym
from sympy import symbols
from sympy import atan2,Abs
from animation import animation
import os


#The following lines take an already converted Gmsh mesh and import it
mesh = Mesh('Mesh.xml')
########################################################################################################################

# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
ffc_options = {'optimize': True}
########################################################################################################################


#Auxillary function spaces for projection and plotting purposes
S = FunctionSpace(mesh,'R',0)
V = VectorFunctionSpace(mesh,'P',1)
######################################################################


d = mesh.geometry().dim()
########################################################################################################################

# Facet and volume functions
Volume = MeshFunction('size_t' , mesh , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh = MeshFunction('size_t', mesh , 'Mesh_facet_region.xml')  #saves the boundary info of the mesh
########################################################################################################################

#Home PC or HPC directory
os.chdir("/lustre/scratch/nmirzaei")
########################################################################################################################


#Saving the volume and boundary
file=File('Axial/Volume.pvd')
file<<Volume
file=File('Axial/bnd.pvd')
file<<bnd_mesh
########################################################################################################################


# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
########################################################################################################################


# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)               # Solution function
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
#SinT_0_code=sym.printing.ccode(SinT_0)
#CosT_0_code = sym.printing.ccode(CosT_0)
b_i=Expression((u_s_code,v_s_code,w_s_code), degree=0)
#sint_i = Expression(SinT_0_code,degree=0)
#cost_i = Expression(CosT_0_code ,degree=0)
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
#SinT_1_code=sym.printing.ccode(SinT_1)
#CosT_1_code = sym.printing.ccode(CosT_1)
b_m=Expression((u_q_code,v_q_code,w_q_code), degree=0)
#sint_m = Expression(SinT_1_code,degree=0)
#cost_m = Expression(CosT_1_code ,degree=0)
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
#SinT_2_code=sym.printing.ccode(SinT_2)
#CosT_2_code = sym.printing.ccode(CosT_2)
b_a = Expression((u_p_code,v_p_code,w_p_code), degree=0)
#sint_a = Expression(SinT_2_code,degree=0)
#cost_a = Expression(CosT_2_code ,degree=0)
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


#Files for saving animation
vtkfile1 = File('Axial/Vol_rec.pvd')
vtkfile2 = File('Axial/bnd_rec.pvd')
vtkfile3 = File('Axial/domain_rec.pvd')
vtkfile4 = File('Axial/Vol_int_rec.pvd')
vtkfile5 = File('Axial/Bnd_intima_rec.pvd')
########################################################################################################################

#looping prameters
counter=0
kk = 0.0
epsilon, sigma, t = 0.0, 16.0, 0.0
deps, dsig, dt= 0.001, 8.0, 0.001
epsmax, sigmax = 3.0, 16.0
ID = 0
######################################################################

while int(epsilon*1000)<=int(epsmax*1000):


    #Growth in the intima with the unknown parameters
    EXP = Expression('eps*exp(-a*(x[2]-0.5)*(x[2]-0.5)-a*(x[1]+0.8)*(x[1]+0.8)-a*(x[0])*(x[0]))',degree=0,a=4.25,eps=2.7*epsilon)
    g_i = as_tensor([[1,0,0],[0,1,0],[0,0,1+EXP]])
    G_i = (1+EXP)
    ginv_i = as_tensor([[1,0,0],[0,1,0],[0,0,1/(1+EXP)]])
    ########################################################################################################################

    # Kinematics for initma
    II = Identity(3)            # Identity tensor
    F =II + grad(u)             # Deformation gradient
    J = det(F)
    Fe_i= F*ginv_i
    Fe_i = variable(Fe_i)
    C_i = F.T*F                   # Right Cauchy-Green tensor
    Ce_i= Fe_i.T*Fe_i
    ####################################################################################################################

    # Kinematics for media
    Fe_m= F*ginv_m
    Fe_m = variable(Fe_m)
    C_m = F.T*F                   # Right Cauchy-Green tensor
    Ce_m= Fe_m.T*Fe_m
    ####################################################################################################################


    # Kinematics for adventitia
    Fe_a= F*ginv_a
    Fe_a = variable(Fe_a)
    C_a = F.T*F                   # Right Cauchy-Green tensor
    Ce_a= Fe_a.T*Fe_a
    ####################################################################################################################


    # Invariants of deformation tensors (intima)
    I_i = tr(C_i)
    Ie_i=tr(Ce_i)
    J_i = det(Fe_i)   #Jacobian of F or Fe?
    ####################################################################################################################


    # Invariants of deformation tensors (media)
    I_m = tr(C_m)
    Ie_m=tr(Ce_m)
    J_m = det(Fe_m)   #Jacobian of F or Fe?
    ####################################################################################################################


    # Invariants of deformation tensors (media)
    I_a = tr(C_a)
    Ie_a=tr(Ce_a)
    J_a = det(Fe_a)   #Jacobian of F or Fe?
    ####################################################################################################################

    #for hopzafel (intima)
    I4_i = conditional(lt(dot(b_i, Ce_i*b_i),Constant(1)),1,dot(b_i, Ce_i*b_i))
    #for hopzafel (media)
    I4_m = conditional(lt(dot(b_m, Ce_m*b_m),Constant(1)),1,dot(b_m, Ce_m*b_m))
    #for hopzafel (adventitia)
    I4_a = conditional(lt(dot(b_a, Ce_a*b_a),Constant(1)),1,dot(b_a, Ce_a*b_a))
    ######################################################################

    #printing the current loop variables
    print(epsilon,flush=True)
    print(t,flush=True)
    print(sigma,flush=True)
    ######################################################################

    # Stored strain energy density (neo-Hookean Hopzafel model)
    psi_i = G_i*((mu_i/2)*(Ie_i - 3)+(eta_i/beta_i)*(exp(beta_i*(rho_i*(I4_i-1)**2+(1-rho_i)*(Ie_i-3)**2))-1)+((mu_i*nu_i)/(1-2*nu_i))*(J_i-1)**2- mu_i*ln(J_i))
    psi_m = ((mu_m/2)*(Ie_m - 3)+((mu_m*nu_m)/(1-2*nu_m))*(J_m-1)**2- mu_m*ln(J_m)+(eta_m/beta_m)*(exp(beta_m*(rho_m*(I4_m-1)**2+(1-rho_m)*(Ie_m-3)**2))-1))
    psi_a = ((mu_a/2)*(Ie_a - 3)+((mu_a*nu_a)/(1-2*nu_a))*(J_a-1)**2- mu_a*ln(J_a)+(eta_a/beta_a)*(exp(beta_a*(rho_a*(I4_a-1)**2+(1-rho_a)*(Ie_a-3)**2))-1))
    ######################################################################

    #Stress variables for the weak form
    TT_i = G_i*diff(psi_i,Fe_i)*(ginv_i.T)
    TT_m = diff(psi_m,Fe_m)*(ginv_m.T)
    TT_a = diff(psi_a,Fe_a)*(ginv_a.T)
    ######################################################################

    #Pulling back the normals using Nanson's formula
    NansonOp = J*inv(F).T
    deformed_NN = dot(NansonOp,n)
    NormN = sqrt(dot(deformed_NN,deformed_NN))
    deformed_N = as_vector([deformed_NN[0],deformed_NN[1],deformed_NN[2]])
    ######################################################################

    #The weak form
    Pi = inner(TT_i,grad(v))*dx(12)+inner(TT_m,grad(v))*dx(13)+inner(TT_a,grad(v))*dx(14)#+(sigma)*dot(deformed_N,v)*ds(8)
    ########################################################################################################################

    #Solve
    solve(Pi == 0, u ,
          solver_parameters={"newton_solver": {"relative_tolerance": 9e-10,
                                               "absolute_tolerance": 9e-10,"maximum_iterations": 80}})
    ######################################################################

     #making sure we start saving results after the pressure is imposed (if there is any pressure). Also making sure we save results every 10 steps
    if sigma<=sigmax:
         sigma+=dsig
    if sigma>sigmax:
         #Save animation every 10 steps to save space and time
         if counter%10==0:
             animation(kk,u,t,vtkfile1,vtkfile2,vtkfile3,vtkfile4,vtkfile5)
         ######################################################################

         #Save the displacements for later stress computations
         File("Opening_Angle/displacement/u%d.xml" %(counter))<<project(u,V) #Save the displacements for later process
         ######################################################################

         #Update loop info
         epsilon+=deps
         counter+=1
         t+=dt
         ID+=1
         ######################################################################

     ######################################################################
    if ID!=0:
        sigma-=dsig
        ID=0
     ######################################################################

     ######################################################################
    kk+=1
     ######################################################################
