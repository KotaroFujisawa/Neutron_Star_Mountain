import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.integrate import quad
#Central density
rhoc = 1.0e17
#Constants
G = 6.674e-11
R = 1.19e4
c = 3e8
#Harmonics
l1 = 2
beta1 = l1*(l1+1)
#Range of the radial coordinate
r_min = 1e-12
r_max = R
n_points = 1000
r = np.linspace(r_min, r_max, n_points)
#Background density
def rhof(r): 
	x = (r/R)
	rhof = rhoc*(np.sinc(x))
	return rhof

#Background pressure
def pf(r):
    x = (r/R)
    pf = ((2*G*(rhoc**2)*(R**2))/np.pi)*((np.sinc(x))**2)
    return pf

#Square of the speed of sound
def csf(r): 
	x = (r/R)
	csf = 10*G*rhoc*(R**2)*((1)/(3*np.pi))*(np.sinc(x))
	return csf
	
#Gravitational acceleration
def gf(r): 
	x = (r/R)
	gf = 4*G*R*rhoc*((np.cos(x)/x) - (np.sinc(x)/x))
	return gf
	
#Gradient of background density
def rhogradf(r): 
	x = (r/R)
	rhogradf = rhoc*(np.pi/R)*((np.cos(x)/x) - (np.sinc(x)/x))
	return rhogradf

#Gradient of square of the speed of sound
def csgradf(r): 
	x = (r/R)
	csgradf = (10/3)*G*rhoc*R*((np.cos(x)/x) - (np.sinc(x)/x))
	return csgradf
	
rho = rhof(r)
cs = csf(r)
g = gf(r)
rhograd = rhogradf(r)
csgrad = csgradf(r)
#Adiabatic factor
Gamma = 2.2
Lambda1 = 3.2
Lambda2 = 4.2	
#The coupled ODEs
def eigenvalue_problem(r, y, om2):
	DPhi = y[0]
	Phi = y[1]
	W = y[2]
	V = y[3]
	dDPhidr = -(2/r)*DPhi + (((beta1/r)**2) - ((4*np.pi*G*rho)/cs))*Phi -	4*np.pi*G*(rhograd + ((rho*g)/(cs)))*W + ((4*np.pi*G*rho*(om2[0])*r)/cs)*V
	dPhidr = DPhi
	dWdr = (1/(cs))*Phi + ((g/cs) - (2/r))*W + (((beta1**2)/r) - (((om2[0])*r)/(cs)))*V
	dVdr = ((g/(cs*(om2[0])*r)) + ((rhograd)/(rho*(om2[0])*r)))*Phi + ((1/r) + ((Lambda2*g*rhograd)/((om2[0])*r*rho)) + (((g**2)*Lambda1)/(cs*r*(om2[0]))) + ((rhograd*csgrad)/(rho*r*(om2[0]))) + ((cs*(rhograd**2))/((rho**2)*(om2[0])*r)) + ((g*csgrad)/(cs*r*(om2[0]))))*W -  (((g)/(cs)) + (1/r) + (rhograd/rho))*V
	return np.vstack((dDPhidr, dPhidr, dWdr, dVdr))
	
#boundary conditions
def bc(y0, yR, om2):
	R = 1.19e4
	return np.array([y0[0], y0[1], y0[2], yR[0] + (3/R)*yR[1], yR[1] - om2[0]*R*yR[3]])

#initial solution guess
om2 = np.zeros(1)
y_guess = np.zeros((4, n_points))
om2[0] = 0.1
#Solution of the BVP
res_first = solve_bvp(eigenvalue_problem, bc, r, y_guess, om2)
#1st g mode
gom1 = res_first.p
y_guess = res_first.sol(r)
#Solution of the BVP, 2nd g-mode
res_second = solve_bvp(eigenvalue_problem, bc, r, y_guess, gom1)
#2nd g mode
gom2 = res_second.p
print("The first eigenfrequency is: ", np.sqrt(gom1))
print("The second eigenfrequency is: ", np.sqrt(gom2))
#Calculate the perturbed density
def Drho1(r): 
	Drho1 = -rho(r)*((res_first[1]/cs(r)) + ((g(r)*res_first[2])/cs(r)) - ((gom1*r*res_first[3])/cs(r)) + ((res_first[2]*rhograd(r))/rho(r))) 
	return Drho1

def Drho2(r): 
	Drho2 = -rho(r)*((res_second[1]/cs(r)) + ((g(r)*res_second[2])/cs(r)) - ((gom2*r*res_second[3])/cs(r)) + ((res_second[2]*rhograd(r))/rho(r)))
	return Drho2

#Multipole moments
I1, error = quad(Drho1(r)*(r**4), r_min, r_max)
I2, error = quad(Drho2(r)*(r**4), r_min, r_max)
#Orbital velocity
OM_min = 0.001
OM_max = 0.1
OM = np.linspace(OM_min, OM_max, n_points)
#Love numbers
def k1(OM): 
	k1 = ((2*np.pi*G)/(5*(R**5)))*(I1**2)/((gom1**2) - (2*(OM**2)))
	return k1

def k2(OM): 
	k2 = ((2*np.pi*G)/(5*(R**5)))*(I2**2)/((gom2**2) - (2*(OM**2)))
	return k2

k = k1 + k2 
plt.plot(OM, k)
plt.show()
