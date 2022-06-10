# %%
import numpy as np
from scipy import optimize, integrate, interpolate
from collections.abc import Iterable
# ***from collections import Iterable*** will be removed after Python 3.10
# and this feature is induced since Python 3.3 .
import os
import matplotlib.pyplot as plt
#%%
"""
    It should be known that in this programme, all the load is represented 
    by Characteristic Value and if you want to add up load, you need to
     use the correct combination yourself.
    And no impact factor is included in the value.
    

    Load system notes:
        I divide the load into two categories: dead loads and live loads. And they
        are further divided into XXXX concentrated loads and XXXX continuous distributed loads.

        Concentrated loads are represented by NumPy Ndarray objects in kN and continuous distributed loads
        are represented by Python functions in kN/m.

        Due to the mechanism of SciPy/NumPy/Matplotlib functions, you'll have to make your load functions adaptable to iterable
        objects(especially for boolean operations).

    Global variables:
        A_total
        H_g
        l1
        f
        m

        They are used everywhere in the program that you need to be careful about them.

    Drawbacks:
        I did not divide live loads into further classes like vehicular loads , ped. loads, wind loads, etc.
        because quasi-permanent coefficients of Chinese specifications for P, V loads are both 0.4. And in
        my project, these are the only two live loads considered.

        And global variables are used in the programme that it can be replaced by local variables, but I'm
        too lazy to do that.
    
"""
def y1_x(x, l1, f, m):
    k = np.arccosh(m)
    return f/(m-1) * (np.cosh(k*x/l1) - 1)
def dy1dx(x, l1, f, m):
    k = np.arccosh(m)
    return f*k/ (m-1)/l1 * np.sinh(k*x/l1)
def phi_x(x, l1, f, m):
    k = np.arccosh(m)
    return np.arctan(f*k/ (m-1)/l1 * np.sinh(k*x/l1))

def DL_arch(x, l1, f, m):
    """
        Calculate the dead load caused by the arch ring (continuous distributed)
        kN/m
        Characteristic Value
    """
    return A_total*25.5 / np.cos(phi_x(x, l1, f, m))

def DL(x, l1, f, m):
    return DL_arch(x, l1, f, m) + DL_spandrel(x, l1, f, m)

def DL_spandrel(x, l1, f, m):
    gd = 139.695 # kN/m
    if(isinstance(x, Iterable)):
        temp_func = lambda xx: gd + (y1_x(xx, l1, f, m) + d/2 - d/(2*np.cos(phi_x(xx, l1, f, m)) ) ) * 9 * 22.5 if(xx<l1-24) else 0
        return np.array(
            [temp_func(xx) for xx in x ]
        )
    else:
        if(x>l1-3*6):
            return 0
        return gd + (y1_x(x, l1, f, m) + d/2 - d/(2*np.cos(phi_x(x, l1, f, m)))) * 9 * 22.5


def LL(x, l1, f, m):
    # vehicular load :10.5kN/m
    # pedestrian load: 12kN/m
    if(isinstance(x, Iterable)):
        temp_func = lambda x: 10.5+12 if(x<l1-24) else 0
        return np.array([temp_func(xx) for xx in x])
    else:
        if(x>l1-18):
            return 0
        else:
            return 10.5+12


l0 = 70
f0 = 70/5.5
d = 1.5
m = 1.5
A_side = 836600e-6 # x2  m2
A_mid = 701200e-6  # x4  m2
A_wet = 377800e-6  # x5  m2

A_total = A_side*2+A_mid*4+A_wet*5  # m2
LL_COE = 0.0  # the coefficient for live loads, quasi-permanent combination uses 0.4
IntegralStep = 1e-2 # integration step for solvers
old_fy = 100
M_j = 10000
M_q = 1
l1 = l0/2
f = f0
oldratio = 100
oldl1 = 100
# %%
#while(np.abs (old_fy - M_j/M_q)>0.005):
for _ in range(10):

    #print("\tGetting the calculation span...")
    for _ in range(10):     # iterate for ten times to get the approximate solution
        phi_j = phi_x(l1, l1, f, m) # approximate phi_j

        f = f0+d/2 - d/2 * np.cos(phi_j)
        l1 = l0/2 + d/2 * np.sin(phi_j)
        #print(f"{f=:.8f}, {l1=:.8f}")

    #print("Calculation span is shown as the last one of the above. If you are not satisfied with the accuracy, you can change the time of iteration in the code.")
    print(f"{l1=:.3f}m")
    x_con_dl = np.array([l1-6, l1-12, l1-18, l1-24])
    F_dl = np.zeros(x_con_dl.shape)
    for i in range(len(x_con_dl)):
        # quasi-permanent combination:
        #   1.0 for dead load, 0.4 for vehicular load and 
        #   pedestrian load
        F_dl[i] = 0.3**2 * (y1_x(x_con_dl[i], l1, f, m)+d/2+0.19) * 25 + \
            (653.88 ) # 423.36+12+218.52 = 653.88
    F_dl[i] -=  653.88/2  # the pier on the most right only bear a quarter of the load
        # and there are two piers, so it has to be subtracted by 


    x_con_ll = np.array([l1-6, l1-12, l1-18, l1-24, 0])
    F_ll = np.zeros(x_con_ll.shape)
    for i in range(len(x_con_ll)):
        F_ll[i] = (670+72)  # characteristic values of vehicular loads and pedestrian loads
    # and notice that the last : F_ll[-1] is not loads transmitted by piers, but 
    #   the concentrated load defined in Code
    F_ll[i-1] -= (670+72)/2
    F_ll[i] = 360/2  # -> the concentrated load of the Class-I vehicular load for the whole bridge

    # Calculate the moment about the skewback caused by external loads, 
    #   using quasi-permanent combination, where 1.0 for dead loads and 0.4 for live loads (P, V)

    temp_func1 = lambda x: DL(x, l1, f, m)*(l1-x)
    temp_func2 = lambda x: LL(x, l1, f, m)*(l1-x) # these two functions is to calculate (static) moment
    ret_val1 = integrate.quad(temp_func1, 0, l1)
    ret_val2 = integrate.quad(temp_func2, 0, l1)
    #print(f"\tTwo errors: {ret_val1[1]:.3f}, {ret_val2[1]:.3f} (unit:kN*m). If they are too big, you'll have to check your code.")
    M_j = ((F_dl*(l1-x_con_dl)).sum()) + LL_COE*((F_ll*(l1-x_con_ll)).sum()) + \
        ret_val1[0] + LL_COE*ret_val2[0]
    H_g = M_j/f
    print(f"{M_j=:.3f} kN*m\n{H_g=:.3f}kN")

    I_dl = x_con_dl<l1/2
    I_ll = x_con_ll<l1/2
    temp_func1 = lambda x: DL(x, l1, f, m)*(l1/2-x)
    temp_func2 = lambda x: LL(x, l1, f, m)*(l1/2-x) # these two functions is to calculate (static) moment
    ret_val1 = integrate.quad(temp_func1, 0, l1/2)
    ret_val2 = integrate.quad(temp_func2, 0, l1/2)
    #print(f"\tTwo errors: {ret_val1[1]:.3f}, {ret_val2[1]:.3f} (unit:kN*m). If they are too big, you'll have to check your code.")
    M_q = ((F_dl[I_dl]*(l1/2-x_con_dl[I_dl])).sum()) + LL_COE*((F_ll[I_ll]*(l1/2-x_con_ll[I_ll])).sum()) + \
        ret_val1[0] + LL_COE*ret_val2[0]
    y_q = M_q/H_g
    print(f"{M_q=:.3f} kN*m\n{y_q=:.3f} m")

    m = 0.5 * (f/y_q-2)**2 - 1
    print(f"{m=:.5f}")
    print("\n")
    if(np.abs(oldratio-f/y_q)<0.0025 and np.abs(oldl1-l1)<0.001):
        print(f"The result error is less than a half level!")
        break
    oldl1 = l1
    oldratio = f/y_q



# %%
# The parameters of the arch bridge have been obtained:
# m, f, l1
# loads: DL, LL, F_con_ll, F_con_dl, x_con_ll, x_con_dl
# And here we are to calculate the internal forces (first, under dead load)

# First of all, we have to obtain the elastic center of the
# arch bridge.

# No deviation, no elastic compression:

temp_func1 = lambda x: DL(x, l1, f, m)*(l1-x)
temp_func2 = lambda x: LL(x, l1, f, m)*(l1-x) # these two functions is to calculate (static) moment
ret_val1 = integrate.quad(temp_func1, 0, l1)
ret_val2 = integrate.quad(temp_func2, 0, l1)
#print(f"\tTwo errors: {ret_val1[1]:.3f}, {ret_val2[1]:.3f} (unit:kN*m). If they are too big, you'll have to check your code.")
M_j = ((F_dl*(l1-x_con_dl)).sum()) + LL_COE*((F_ll*(l1-x_con_ll)).sum()) + \
    ret_val1[0] + LL_COE*ret_val2[0]
H_g = M_j/f # kN

x = np.r_[0:l1:0.1]
N_dl = H_g/np.cos(phi_x(x, l1, f, m))
plt.plot(x, N_dl)
plt.gca().invert_xaxis()
plt.xlabel("x(m)")
plt.ylabel("Axial Force (kN, No Compression)")
plt.grid()
plt.show()


# %%
# considering deviation (only dead loads):

V = F_dl.sum() + integrate.quad(lambda x: DL(x, l1, f, m), 0, l1)[0]
def filled_seg_func(x, y):
    ret = np.zeros((2,))
    ret[0] = y[1]
    ret[1] = 1/H_g*DL(x, l1, f, m)
    return ret


filled_ode = integrate.solve_ivp(filled_seg_func, (0, l1-24), [0, 0], t_eval=np.r_[0:l1-24:IntegralStep])
# plt.plot(filled_ode.t, filled_ode.y[0], label="Rational Arch Axis 2")
# plt.plot(x, y1_x(x, l1, f, m), label="Arch Axis from 5P Method")
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.legend()
# plt.show()
def open_seg_func(x, y):
    ret = np.zeros(y.shape)
    ret[0] = y[1]
    ret[1] = 1/H_g*DL(x, l1, f, m)
    return ret
"""# %%
# open segment bvp
# It turns out that this is not a bvp (boundary value problem)

def open_seg_func(x, y):
    ret = np.zeros(y.shape)
    ret[0] = y[1]
    ret[1] = 1/H_g*DL(x, l1, f, m)
    return ret

def open_seg_bc(ya, yb):
    return np.array([ya[0]-filled_ode.y[0, -1], yb[0]-f])

open_ode = integrate.solve_bvp(open_seg_func, open_seg_bc, 
    np.r_[l1-24:l1:IntegralStep], 
    [y1_x(np.r_[l1-24:l1:IntegralStep], l1, f, m), np.tan(phi_x(np.r_[l1-24:l1:IntegralStep], l1, f, m))])

RationalX = np.hstack((filled_ode.t, open_ode.x))
RationalY = np.hstack((filled_ode.y[0], open_ode.y[0]))

x = np.r_[0:l1:0.1]
plt.plot(RationalX, RationalY, label="Rational Arch Axis")
plt.plot(x, y1_x(x, l1, f, m), label="Arch Axis from 5P Method")
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.legend()
plt.show()"""
# %%
OpenBoundary = [l1-24, l1-18, l1-12, l1-6, l1]
# open segment ivp
lastY = filled_ode.y[0, -1]
open_ode_list = []
openX = np.array([])
openY = np.array([])
for i in range(len(OpenBoundary)-1):
    I = (x_con_dl)<=OpenBoundary[i] + 0.1
                # plus this 0.1 to avoid the problem that 
                # due to discontinuity of the data. And 0.1 is 
                # small enough, or you can use 1e-5, but I think
                # 0.1 is big enough to cover the discontinuity, and
                # too small to cause errors
    tanphi_0 = ( (F_dl[I].sum() 
    + integrate.quad(lambda x: DL(x, l1, f, m), 0, OpenBoundary[i])[0])/
    H_g) 
    open_ode_list.append(integrate.solve_ivp(open_seg_func, [OpenBoundary[i], OpenBoundary[i+1]], 
    [lastY, tanphi_0], t_eval = np.r_[OpenBoundary[i]: OpenBoundary[i+1]:IntegralStep]))
    lastY = open_ode_list[-1].y[0, -1]
    #plt.plot(open_ode_list[-1].t, open_ode_list[-1].y[0], c='b')
    openX = np.hstack((openX, open_ode_list[-1].t))
    openY = np.hstack((openY, open_ode_list[-1].y[0]))

RationalX = np.hstack((filled_ode.t, openX))
RationalY = np.hstack((filled_ode.y[0], openY))

x = np.r_[0:l1:0.1]
plt.plot(RationalX, RationalY, label="Rational Arch Axis")
plt.plot(x, y1_x(x, l1, f, m), label="Catenary Arch Axis from 5P Method")
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.legend()
plt.show()
# %%

# elastic center

temp_func = lambda x: y1_x(x, l1, f, m)*np.sqrt(1+dy1dx(x, l1, f, m)**2)

numerator = integrate.quad(temp_func, 0, l1)[0]
temp_func = lambda x: np.sqrt(1+dy1dx(x, l1, f, m)**2)
denominator = integrate.quad(temp_func, 0, l1)[0]

y_s = numerator/denominator

# %%
RationalY_x = interpolate.interp1d(RationalX, RationalY, bounds_error=False, 
    fill_value='extrapolate')

dy_x = lambda x:(y1_x(x, l1, f, m) - RationalY_x(x))
# secondary forces
temp_func = lambda x: (y1_x(x, l1, f, m) - RationalY_x(x))*np.sqrt(1+dy1dx(x, l1, f, m)**2)

numerator = integrate.quad(temp_func, 0, l1)[0]
temp_func = lambda x: np.sqrt(1+dy1dx(x, l1, f, m)**2)
denominator = integrate.quad(temp_func, 0, l1)[0]

dX1 = -H_g * numerator/denominator


temp_func = lambda x: dy_x(x)*\
    (y_s - y1_x(x, l1, f, m) ) * np.sqrt(1+dy1dx(x, l1, f, m)**2)
numerator = integrate.quad(temp_func, 0, l1)[0]

temp_func = lambda x: (y_s-y1_x(x, l1, f, m))**2*np.sqrt(1+dy1dx(x, l1, f, m)**2)
denominator = integrate.quad(temp_func, 0, l1)[0]
dX2 = H_g * numerator/denominator

# %%

# deviation bending moment
# y = ys - y1
dM_x = lambda x: dX1 - dX2*(y_s - y1_x(x, l1, f, m)) + H_g*dy_x(x)

x = np.r_[0:l1:0.1]
plt.plot(x, dM_x(x))
plt.axhline(c='k')
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x (m)")
plt.ylabel("ΔM (kN*m)")
plt.show()
# %%

# elastic compression


# 不会


# %%
# Coder Note:
# IntegrationWarning: The maximum number of subdivisions (50) has been achieved.
# the reason for the above warning I think, is the discontinuity 
# of the rational arch axis (I guess). I don't really know actually.
