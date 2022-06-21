# %%
import numpy as np
from scipy import optimize, integrate, interpolate
from collections.abc import Iterable
# ***from collections import Iterable*** will be removed after Python 3.10
# and this feature is induced since Python 3.3 .
import os
import matplotlib.pyplot as plt
import pandas as pd
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
    return (DL_arch(x, l1, f, m) + DL_spandrel(x, l1, f, m))

def DL_spandrel(x, l1, f, m):
    gd = 139.695 # kN/m
    if(isinstance(x, Iterable)):
        temp_func = lambda xx: gd + \
            (y1_x(xx, l1, f, m) + d/2 - d/(2*np.cos(phi_x(xx, l1, f, m)) ) )\
                 * 9 * 22.5 if(xx<l1-18) else 0
        return np.array(
            [temp_func(xx) for xx in x ]
        )
    else:
        if(x>l1-18):
            return 0
        return gd + (y1_x(x, l1, f, m) + d/2 - d/(2*np.cos(phi_x(x, l1, f, m)))) * 9 * 22.5


def LL(x, l1, f, m):
    # vehicular load :10.5kN/m
    # pedestrian load: 12kN/m
    if(isinstance(x, Iterable)):
        temp_func = lambda x: 2*10.5+12 if(x<l1-18) else 0
        return np.array([temp_func(xx) for xx in x])
    else:
        if(x>l1-18):
            return 0
        else:
            return 2*10.5+12


l0 = 70
f0 = 70/5.5
d = 1.5
m = 1.5
A_side = 836600e-6 # x2  m2
A_mid = 701200e-6  # x4  m2
A_wet = 377800e-6  # x5  m2
I_total = 1.786519# m4
I_new = I_total + np.pi/4 * 25e-3**2 * 0.65**2 * 15*2*6 * (2e5/3.25e4)
A_total = A_side*2+A_mid*4+A_wet*5  # m2
DL_COE = 1.0
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
with pd.ExcelWriter("Iteration_Steps.xlsx") as writer:
    for IterationTime in range(30):

        #print("\tGetting the calculation span...")
        for _ in range(20):     # iterate for ten times to get the approximate solution
            phi_j = phi_x(l1, l1, f, m) # approximate phi_j

            f = f0+d/2 - d/2 * np.cos(phi_j)
            l1 = l0/2 + d/2 * np.sin(phi_j)
            #print(f"{f=:.8f}, {l1=:.8f}")

        #print("Calculation span is shown as the last one of the above. If you are not satisfied with the accuracy, you can change the time of iteration in the code.")
        print(f"{l1=:.3f} m, {f=:.3f} m")
        x_con_dl = np.array([l1-6, l1-12, l1-18])
        F_dl = np.zeros(x_con_dl.shape)
        for i in range(len(x_con_dl)):
            # quasi-permanent combination:
            #   1.0 for dead load, 0.4 for vehicular load and 
            #   pedestrian load
            F_dl[i] = 0.3**2 * (y1_x(x_con_dl[i], l1, f, m)+d/2+0.19) * 25 * 2.000 + \
                (653.88 ) # 423.36+12+218.52 = 653.88
        F_dl[i] -=  653.88/2  # the pier on the most right only bear a quarter of the load
            # and there are two piers, so it has to be subtracted by 
        F_dl += 25*9*0.5*0.4 # plus the pier cap weight
        F_dl += 25*0.4*9*0.2
        x_con_diaphram = np.r_[0:32+1:4]
        F_diaphram = np.zeros(x_con_diaphram.shape)
        for i in range(len(x_con_diaphram)):
            F_diaphram[i] = 10.5474
        F_diaphram[0] /= 2
        
        x_con_dl = np.hstack((x_con_dl,x_con_diaphram))
        F_dl = np.hstack((F_dl, F_diaphram))

        x_con_ll = np.array([l1-6, l1-12, l1-18, 0])
        F_ll = np.zeros(x_con_ll.shape)
        for i in range(len(x_con_ll)):
            F_ll[i] = (10.5*6*2+12 * 6)  # characteristic values of vehicular loads and pedestrian loads
        # and notice that the last : F_ll[-1] is not loads transmitted by piers, but 
        #   the concentrated load defined in Code
        F_ll[i-1] -= (10.5*6*2+12*6)/2
        F_ll[i] = 360*2 /2  # -> the concentrated load of the Class-I vehicular load for the whole bridge


        #F_dl *= DL_COE
        # Calculate the moment about the skewback caused by external loads, 
        #   using quasi-permanent combination, where 1.0 for dead loads and 0.4 for live loads (P, V)

        temp_func1 = lambda x: DL(x, l1, f, m)*(l1-x)
        temp_func2 = lambda x: LL(x, l1, f, m)*(l1-x) # these two functions is to calculate (static) moment
        ret_val1 = integrate.quad(temp_func1, 0, l1)
        ret_val2 = integrate.quad(temp_func2, 0, l1)
        int1 = ret_val1[0]
        print("F_dl =", np.round(F_dl,3))
        print("x_dl =", np.round(x_con_dl,3))
        #print(f"\tTwo errors: {ret_val1[1]:.3f}, {ret_val2[1]:.3f} (unit:kN*m). If they are too big, you'll have to check your code.")
        M_j = DL_COE*((F_dl*(l1-x_con_dl)).sum()) + LL_COE*((F_ll*(l1-x_con_ll)).sum()) + \
            DL_COE*ret_val1[0] + LL_COE*ret_val2[0]
        H_g = M_j/f
        print(f"{M_j=:.3f} kN*m\n{H_g=:.3f}kN")
        I_dl = x_con_dl<l1/2
        I_ll = x_con_ll<l1/2
        temp_func1 = lambda x: DL(x, l1, f, m)*(l1/2-x)
        temp_func2 = lambda x: LL(x, l1, f, m)*(l1/2-x) # these two functions is to calculate (static) moment
        ret_val1 = integrate.quad(temp_func1, 0, l1/2)
        ret_val2 = integrate.quad(temp_func2, 0, l1/2)

        int2 = ret_val1[0]
        #print(f"\tTwo errors: {ret_val1[1]:.3f}, {ret_val2[1]:.3f} (unit:kN*m). If they are too big, you'll have to check your code.")
        M_q = DL_COE*((F_dl[I_dl]*(l1/2-x_con_dl[I_dl])).sum()) + LL_COE*((F_ll[I_ll]*(l1/2-x_con_ll[I_ll])).sum()) + \
            DL_COE*ret_val1[0] + LL_COE*ret_val2[0]
        y_q = M_q/H_g
        print(f"{M_q=:.3f} kN*m\n{y_q=:.3f} m")

        m = 0.5 * (f/y_q-2)**2 - 1
        DICT = {"F":F_dl, "x":x_con_dl,"F_q":F_dl[I_dl], "x_q":x_con_dl[I_dl],\
             "M_j":M_j, "M_q":M_q, "l1":l1, "f":f, "y_q":y_q,
             "m":m}
        #df = pd.DataFrame(DICT).transpose()
        table = np.ndarray((len(F_dl)+1, 10), object)
        table[0, 0] = "F(kN)"
        for i in range(0, len(F_dl)):
            table[i+1, 0] = f"{F_dl[i]:.3f}"

        table[0, 1] = "x(m)"
        for i in range(0, len(F_dl)):
            table[i+1, 1] = f"{x_con_dl[i]:.3f}"

        table[0, 2] = "F_q(kN)"
        for i in range(0, len(F_dl[I_dl])):
            table[i+1, 2] = f"{F_dl[I_dl][i]:.3f}"

        table[0, 3] = "x_q(kN)"
        for i in range(0, len(F_dl[I_dl])):
            table[i+1, 3] = f"{x_con_dl[I_dl][i]:.3f}"

        table[0, 4] = "M_j(kN*m)"
        table[1, 4] = f"{M_j:.3f}"
        table[0, 5] = "M_q(kN*m)"
        table[1, 5] = f"{M_q:.3f}"
        table[0, 6] = "f(m)"
        table[1, 6] = f"{f:.3f}"
        table[0, 7] = "y_q(m)"
        table[1, 7] = f"{y_q:.3f}"
        table[0, 8] = "l1(m)"
        table[1, 8] = f"{l1:.3f}"
        table[0, 9] = "m"
        table[1, 9] = f"{m:.3f}"

        table[3, 4] = "Int1(kN*m)"
        table[4, 4] = f"{int1:.3f}"
        table[3, 5] = "Int2(kN*m)"
        table[4, 5] = f"{int2:.3f}"
        df = pd.DataFrame(table)
        df.to_excel(writer, f"Iteration {IterationTime}", header=False, index=False)
        print(f"{m=:.5f}")
        if( m<=1 ):
            #print("m <= 1, you have to redesign the bridge layout!")
            raise Exception("m <= 1, you have to redesign the bridge layout!")
        print("\n")
        if(np.abs(oldratio-f/y_q)<0.005 and np.abs(oldl1-l1)<0.001):
            print(f"The result error is less than a half level!")
            break
        oldl1 = l1
        oldratio = f/y_q


#m -= 0.1
#m = 1.0001
#m = np.sqrt(m)
#m = 3
# m = 1.01
m = 1.35
x = np.linspace(0, l1, 100)
plt.plot(x, DL(x, l1, f, m))
plt.xlabel("x(m)")
plt.ylabel("kN/m")
plt.gca().invert_xaxis()
plt.grid()
plt.title("Dead Load (Continuously Distributed)")
plt.show()


x = np.linspace(0, l1, 20)
y = y1_x(x, l1, f, m )
pd.DataFrame({"x(m) ":x[0:10], "y(m) ":y[0:10], \
    "x(m)":x[10:], "y(m)":y[10:]}).to_excel("Coordination.xlsx")

x = np.linspace(0, l1, 20)
y = y1_x(x, l1, f, m )+d/(2*np.cos(phi_x(x, l1, f, m )))
pd.DataFrame({"x(m) ":x[0:10], "y(m) ":y[0:10], \
    "x(m)":x[10:], "y(m)":y[10:]}).to_excel("Coordination+.xlsx")

x = np.linspace(0, l1, 20)
y = y1_x(x, l1, f, m )-d/(2*np.cos(phi_x(x, l1, f, m )))
pd.DataFrame({"x(m) ":x[0:10], "y(m) ":y[0:10], \
    "x(m)":x[10:], "y(m)":y[10:]}).to_excel("Coordination-.xlsx")
# %%
# Change Load Combination:
LL_COE = 1.4
DL_COE = 1.2
#F_ll[F_ll<200] = 0
# Remember to modify LL()!
integral_accuracy = 1e-3
# The parameters of the arch bridge have been obtained:
# m, f, l1
# loads: DL, LL, F_con_ll, F_con_dl, x_con_ll, x_con_dl
# And here we are to calculate the internal forces (first, under dead load)

# First of all, we have to obtain the elastic center of the
# arch bridge.

# No deviation, no elastic compression:
# I_dl_sort = np.argsort(x_con_dl)
# x_con_dl = x_con_dl[I_dl_sort]
# F_dl = F_dl[I_dl_sort]

# I_ll_sort = np.argsort(x_con_ll)
# x_con_ll = x_con_ll[I_ll_sort]
# F_ll = F_ll[I_ll_sort]

x_con_all = np.hstack((x_con_ll, x_con_dl))
F_all = np.hstack((F_ll*LL_COE, DL_COE*F_dl))

I_all_sort = np.argsort(x_con_all)
x_con_all = x_con_all[I_all_sort]
F_all = F_all[I_all_sort]



temp_func1 = lambda x: DL(x, l1, f, m)*(l1-x)
temp_func2 = lambda x: LL(x, l1, f, m)*(l1-x) # these two functions is to calculate (static) moment
ret_val1 = integrate.quad(temp_func1, 0, l1)
ret_val2 = integrate.quad(temp_func2, 0, l1)
#print(f"\tTwo errors: {ret_val1[1]:.3f}, {ret_val2[1]:.3f} (unit:kN*m). If they are too big, you'll have to check your code.")
M_j = DL_COE*((F_dl*(l1-x_con_dl)).sum()) + LL_COE*((F_ll*(l1-x_con_ll)).sum()) + \
    DL_COE*ret_val1[0] + LL_COE*ret_val2[0]
H_g = M_j/f # kN

x = np.r_[0:l1:0.1]
N_dl = H_g/np.cos(phi_x(x, l1, f, m))
plt.plot(x, N_dl)
plt.gca().invert_xaxis()
plt.xlabel("x(m)")
plt.ylabel("Axial Force (kN, No Compression)")
plt.grid()
plt.title("No Deviation, No Elastic Compression")
plt.show()


# %%
# considering deviation (only dead loads):

def rational_axis_ode(x, y):
    ret = np.zeros(y.shape)
    ret[0] = y[1]
    ret[1] = 1/H_g*(DL_COE*DL(x, l1, f, m)+LL_COE*LL(x, l1, f, m))
    return ret

integral_bound = np.hstack(([0], np.unique(np.round(x_con_all, 3)), [l1]))
lastY = 0
RationalX = np.zeros((0,))
RationalY = np.zeros((0,))

for i in range(len(integral_bound)-1):
    I = x_con_all<=integral_bound[i]+1e-3
    tan0 = (F_all[I].sum() + \
        integrate.quad(lambda x: DL_COE*DL(x, l1, f, m)+LL_COE*LL(x, l1, f, m), \
            0, integral_bound[i])[0]) / H_g

    ode_ret = integrate.solve_ivp(rational_axis_ode, [integral_bound[i], integral_bound[i+1]+1e-6], \
        [lastY, tan0], 'RK45', \
            np.r_[integral_bound[i]:integral_bound[i+1]:integral_accuracy])
    if(len(ode_ret.t)):
        RationalX = np.hstack((RationalX, ode_ret.t))
        RationalY = np.hstack((RationalY, ode_ret.y[0]))
        lastY = RationalY[-1]

plt.plot(RationalX, RationalY, label="Rational Arch Axis")
plt.plot(np.r_[0:l1:1e-2], y1_x(np.r_[0:l1:1e-2], l1, f, m), label = "Catenary Arch Axis from 5P Method")
plt.plot(np.r_[0:l1:1e-2], y1_x(np.r_[0:l1:1e-2], l1, f, m)-0.75/np.cos(phi_x(np.r_[0:l1:1e-2], l1, f, m)), label = "Intrados Arch Web", c = 'k')
plt.plot(np.r_[0:l1:1e-2], y1_x(np.r_[0:l1:1e-2], l1, f, m)+0.75/np.cos(phi_x(np.r_[0:l1:1e-2], l1, f, m)), label = "Extrados Arch Web", c = 'k')

plt.title("Rational Arch Axis")
plt.grid()
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.legend()
plt.show()

RationalY_x = interpolate.interp1d(RationalX, RationalY, 
bounds_error=False, fill_value='extrapolate')
dy_x = lambda x: y1_x(x, l1, f, m )-RationalY_x(x)

# %%
# elastic center:

y_s = integrate.quad(lambda x: y1_x(x, l1, f, m)* np.sqrt(1+dy1dx(x, l1, f, m)**2), 0, l1)[0]\
 / integrate.quad(lambda x:np.sqrt(1+dy1dx(x, l1, f, m)**2), 0, l1)[0]
y_x = lambda x: y_s - y1_x(x, l1, f, m)

# %%
# deviation secondary forces

dX1 = -H_g *\
    integrate.quad(
        lambda x:dy_x(x)*np.sqrt(1+dy1dx(x, l1, f, m)**2), 0, l1)[0]\
    / integrate.quad(lambda x:np.sqrt(1+dy1dx(x, l1, f, m)**2), 0, l1)[0]

dX2 = H_g * \
    integrate.quad(\
        lambda x:y_x(x)*dy_x(x)*np.sqrt(1+dy1dx(x, l1, f, m)**2) , 0, l1)[0]\
    / integrate.quad(lambda x: y_x(x)**2, 0, l1)[0]

dM_x = lambda x: dX1 - dX2*y_x(x) + H_g*dy_x(x)
x = np.r_[0:l1:1e-2]
plt.plot(x, dM_x(x))
plt.gca().invert_xaxis()
plt.axhline(c='k')
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Bending Moment (kN*m)")
plt.title("Deviation Effect")
plt.show()

x = np.linspace(0, l1)
dN_x = lambda x:dX2*np.cos(phi_x(x, l1, f, m ))

plt.plot(x, dN_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Axial Force (kN)")
plt.title("Deviation Effect")
plt.show()

dQ_x = lambda x:dX2*np.sin(phi_x(x, l1, f, m))
plt.plot(x, dQ_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Shear Force (kN)")
plt.title("Deviation Effect")
plt.show()

# %%
# Elastic compression
#I_total/=2
mu = integrate.quad(lambda x: np.cos(phi_x(x, l1, f, m))**2 *np.sqrt(1+dy1dx(x, l1, f, m)**2)/(A_total + (2e5/3.25e4-1)*30*6*np.pi*25e-3**2/4), 0, l1)[0] / \
integrate.quad(lambda x: y_x(x)**2 * np.sqrt(1+dy1dx(x, l1, f, m)**2)/(I_new*0.8), 0, l1)[0]


mu1 = integrate.quad(lambda x:1/np.cos(phi_x(x, l1, f, m))/(A_total + (2e5/3.25e4-1)*6*30*np.pi*25e-3**2/4), 0, l1)[0] * 1/ \
integrate.quad(lambda x: y_x(x)**2 * np.sqrt(1+dy1dx(x, l1, f, m)**2)/(I_new*0.8), 0, l1)[0]

Nc_x = lambda x:-mu1/(1+mu) * H_g * np.cos(phi_x(x, l1, f, m))
Mc_x = lambda x: mu1/(1+mu) * H_g * (y_s - y1_x(x, l1, f, m))
Qc_x = lambda x: mu1/(1+mu) * H_g * np.sin(phi_x(x, l1, f, m))

plt.plot(x, Qc_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Shear Force (kN)")
plt.title("Elastic Compression Effect")
plt.show()

plt.plot(x, Mc_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Bending Moment (kN*m)")
plt.title("Elastic Compression Effect")
plt.show()

plt.plot(x, Nc_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Axial Force (kN)")
plt.title("Elastic Compression Effect")
plt.show()

# %%
# total internal forces
plt.plot(x, Nc_x(x)+dN_x(x)+H_g/np.cos(phi_x(x, l1, f, m)))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Axial Force (kN)")
plt.title("Total Internal Force")
plt.show()

plt.plot(x, Mc_x(x)+dM_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Bending Moment (kN*m)")
plt.title("Total Internal Force")
plt.show()

plt.plot(x, Qc_x(x)+dQ_x(x))
plt.gca().invert_xaxis()
plt.grid()
plt.xlabel("x(m)")
plt.ylabel("Shear Force (kN)")
plt.title("Total Internal Force")
plt.show()

#%%

La = integrate.quad(lambda x: 2*np.sqrt(1+dy1dx(x, l1, f, m)**2), 0, l1)[0]
print(f"Length of arch axis: {La:.3f} m")

x = np.r_[0:l1+0.1:0.1]
y = y1_x(x, l1, f, m)
dat = ""
for i in range(len(x)):
    dat += f"{x[i]:.3f},{y[i]:.3f}\n"
with open("arch_axis.txt", "wt", encoding='utf-8') as file:
    file.write(dat)
# %%
# Coder Note:
# IntegrationWarning: The maximum number of subdivisions (50) has been achieved.
# the reason for the above warning I think, is the discontinuity 
# of the rational arch axis (I guess). I don't really know actually.
