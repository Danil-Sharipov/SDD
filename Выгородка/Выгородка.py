D=3574.4
D/=1000
l=1546.6
l/=1000
a=138.3331245
a/=1000
d=70
d/=1000
delta=1.14*d/2
x_0=27
y_0=48
x_0/=1000
y_0/=1000
d_holes=130
d_holes/=1000
from functools import lru_cache
from matplotlib import pyplot as plt
import numpy as np
import warnings
import scipy.special as sc
from mpmath import*
import math
from tqdm import tqdm
from ansys.mapdl.core import launch_mapdl
import scipy.optimize as optimize
import scipy.optimize as opt
import time
def gr(d,delta,a,l):
    x=l-math.cos(math.pi/6)*a
    y=a/2+3*a*math.sin(math.pi/6)+3*a
    k=math.tan(math.pi/6)
    c=y-k*x
    x_1=x-math.cos(math.pi/6)*a
    y_1=y+a*math.sin(math.pi/6)
    x_1+=(delta+d/2)*math.cos(math.pi/6)
    y_1+=(delta+d/2)*math.sin(math.pi/6)
    @lru_cache()
    def check_delta2(delta2,eps):
        x_2=x_1+(delta2+d)*math.cos(math.pi/6)
        y_2=y_1+(delta2+d)*math.sin(math.pi/6)
        x_3=x_2
        y_3=y_2-(delta2+d)
        y_30=k*x_3+c
        x_30=(y_3-c)/k
        u=(x_30-x_3)
        i=(y_3-y_30)
        delta2_new=2*(u*i/math.sqrt(u**2+i**2)-d/2)
        if abs(delta2_new-delta2)>eps:
            return check_delta2((delta2_new+delta2)/2,eps)
        else:
            return [x_1,y_1,x_2,y_2,x_3,y_3,delta2_new]
    arrr=check_delta2(0,1/1e8)
    x_1,y_1,x_2,y_2,x_3,y_3,delta2=arrr
    x_4=x_3+math.cos(math.pi/3)*(delta2+d)
    y_4=y_3-math.sin(math.pi/3)*(delta2+d)
    x_5=x_4
    y_5=y_4-delta2-d
    x_6=x_5+math.cos(math.pi/6)*(delta2+d)
    y_6=y_5-math.sin(math.pi/6)*(delta2+d)
    x_7=x_6+math.cos(math.pi/3)*(delta2+d)
    y_7=y_6-math.sin(math.pi/3)*(delta2+d)
    x_8=x_7
    y_8=y_7-delta2-d
    return [x_1,y_1,x_2,y_2,x_3,y_3,x_4,y_4,x_5,y_5,x_6,y_6,x_7,y_7,x_8,y_8]
def main(d,delta,arr,show,delta3):
    mapdl = launch_mapdl(loglevel="ERROR")
    mapdl.prep7()
    mapdl.units("SI")
    def Co(D,t):
        return math.sqrt(D**2/4-t**2)
    x=Co(D,y_0)-x_0
    y=0
    k1=mapdl.k(x=x,y=y)
    x=x
    y+=y_0
    k2=mapdl.k(x=x,y=y)
    l1=mapdl.l(k1,k2)
    
    x+=x_0
    y=y
    k3=mapdl.k(x=x,y=y)
    l2=mapdl.l(k2,k3)
    
    x_1=x
    x=D/2*math.cos(math.pi/6)
    y=D/2*math.sin(math.pi/6)
    k4=mapdl.k(x=x,y=y)
    x=(x+x_1)/2
    y=Co(D,x)
    k5=mapdl.k(x=x,y=y)
    l3=mapdl.larc(k3,k4,k5)
    
    x=l-2*math.cos(math.pi/6)*a
    y=math.tan(math.pi/6)*x
    k6=mapdl.k(x=x,y=y)
    l4=mapdl.l(k4,k6)
    
    x+=math.cos(math.pi/6)*a
    y-=math.sin(math.pi/6)*a
    k7=mapdl.k(x=x,y=y)
    l5=mapdl.l(k6,k7)
    
    x=x
    y-=a
    k8=mapdl.k(x=x,y=y)
    l6=mapdl.l(k7,k8)
    
    x+=math.cos(math.pi/6)*a
    y-=math.sin(math.pi/6)*a
    k9=mapdl.k(x=x,y=y)
    l7=mapdl.l(k8,k9)
    
    x=x
    y-=a
    k10=mapdl.k(x=x,y=y)
    l8=mapdl.l(k9,k10)
    
    x+=math.cos(math.pi/6)*a
    y-=math.sin(math.pi/6)*a
    k11=mapdl.k(x=x,y=y)
    l9=mapdl.l(k10,k11)
    
    x=x
    y-=a
    k12=mapdl.k(x=x,y=y)
    l10=mapdl.l(k11,k12)
    
    x-=math.cos(math.pi/6)*a
    y-=math.sin(math.pi/6)*a
    k13=mapdl.k(x=x,y=y)
    l11=mapdl.l(k12,k13)
    
    x=x
    y-=a/2
    k14=mapdl.k(x=x,y=y)
    l12=mapdl.l(k13,k14)
    
    l13=mapdl.l(k14,k1)
    
    
    
    area=mapdl.al('ALL')
    circ_anum = mapdl.cyl4(l+delta3, 0, d_holes/2)
    plate_with_hole_anum = mapdl.asba(area, circ_anum)
    for i in range(int(len(arr)/2)):
        xhole=arr[2*i]
        yhole=arr[2*i+1]
        circ_anum = mapdl.cyl4(xhole, yhole, d/2)
        plate_with_hole_anum = mapdl.asba(plate_with_hole_anum, circ_anum)
    if show==True:
        mapdl.lplot(vtk=True,cpos="xy",notebook=False)
    
    
    mapdl.et(1, 75)
    mapdl.mp("DENS",' ', 7800) 
    mapdl.mp("kxx",' ', 19)
    mapdl.lsel("ALL")
    mapdl.esize(5/1000)
    mapdl.amesh("ALL")
    mapdl.mopt("EXPND", 0.7)
    if show==True:
        mapdl.eplot(vtk=True, cpos="xy", show_edges=True, show_axes=False, line_width=2, background="w",notebook=False)
    num=mapdl.mesh.nodes
    mapdl.lsel("S","LINE",' ',5)
    for i in range(6,13):
        mapdl.lsel("A","LINE",' ',i)
    if show==True:
        mapdl.lplot(vtk=True,cpos="xy",notebook=False)
    mapdl.sfl("ALL","conv",vali=30000,val2i=320)
    mapdl.allsel("ALL","LINE")
    for i in range(5,13):
        mapdl.lsel("U","LINE",' ',i)
    mapdl.lsel("U","LINE",' ',19)
    mapdl.lsel("U","LINE",' ',18)
    mapdl.lsel("U","LINE",' ',20)
    mapdl.lsel("U","LINE",' ',24)
    mapdl.lsel("U","LINE",' ',27)
    mapdl.lsel("U","LINE",' ',29)
    if show==True:
        mapdl.lplot(vtk=True,cpos="xy",notebook=False)
    mapdl.sfl("ALL","conv",vali=3000,val2i=280)
    mapdl.allsel()
    Reff=1700.096/1000
    q_v_max=247769.2736129015*1000
    q_0=q_v_max
    alpha=sc.jn_zeros(0, 1)[0]
    for q in tqdm(range(0,len(num))):
        x, y, z=num[q]
        q+=1
        r=math.sqrt(x**2+y**2)
        q_n=max([q_0*besselj(0,alpha/Reff*r),0])
        mapdl.bf(q,"HGEN",q_n) 
    mapdl.finish()
    mapdl.slashsolu()
    mapdl.solve()
    mapdl.post1()
    mapdl.set(1, 1)
    if show==True:
        mapdl.post_processing.plot_nodal_temperature(notebook=False,cpos="xy")
    result = mapdl.result
    nnum, temp = result.nodal_temperature(0)
    mapdl.exit()
    print(max(temp))
    print(delta*1000)
    print(delta3*1000)
    return max(temp)
def find_delta(x):
    delta=x[0]
    delta3=x[1]
    return main(d,delta,gr(d,delta,a,l),False,delta3)
x0 =[26.279707105888303/1000,0.1]
bnds = ((1/1000, 100/1000),(d_holes/2+1/1000, 40/1000))
print(main(d,x0[0],gr(d,x0[0],a,l),True,x0[1]))
res = optimize.minimize(find_delta,x0, bounds=bnds)
print(res.x)
print(main(d,res.x[0],gr(d,res.x[0],a,l),True,res.x[1]))