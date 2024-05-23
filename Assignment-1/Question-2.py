import numpy as np
from functions import *

import numpy as np


E1 = 145.3e3
E2 = 8.5e3
G12 = 4.58e3
v12 = 0.31
theta = [0,90,+45,-45,-45,+45,90,0,0,90,+45,-45,-45,+45,90,0]
thickness = [0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]
properties_o = [[E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12]]
Nx = 0
Ny = 0
Ns = 0
Mx = 0
My = 0
Ms = 0
Force_per_length = np.array([[Nx],[Ny],[Ns],[Mx],[My],[Ms]])

strain_fpf_x_MS = []
strain_fpf_y_MS = []
strain_fpf_s_MS = []
strain_lpf_x_MS = []
strain_lpf_y_MS = []
strain_lpf_s_MS = []

strain_fpf_x_PK = []
strain_fpf_y_PK = []
strain_fpf_s_PK = []
strain_lpf_x_PK = []
strain_lpf_y_PK = []
strain_lpf_s_PK = []



fpf_Ns_plot_MaxStress = []
fpf_Ny_plot_MaxStress = []
lpf_Ns_plot_MaxStress = []
lpf_Ny_plot_MaxStress = []
fpf_Ns_plot_Puck = []
fpf_Ny_plot_Puck = []
lpf_Ns_plot_Puck = []
lpf_Ny_plot_Puck = []



def failure_envelope_puck(thickness, theta, properties, Force_per_length, Xt, Xc, Yt, Yc, S, r, the):
    failed_laminates = []
    shear_failure = []
    fpf = False
    lpf = False
    while fpf == False or lpf == False:
        A, B, D, z = ABD_matrix_differentMaterial(thickness, theta, properties)
        stress_laminas = stress_per_lamina_differentMaterials(Force_per_length, theta, A, B, D, properties, z)
        strain_laminas = strain_per_lamina_differentMaterials(Force_per_length, theta, A, B, D, properties, z)
        failure_index = Puck_envelope(stress_laminas, Xt, Xc, Yt, Yc, S)
        for i, sub_failure_index in enumerate(failure_index):
            if i in failed_laminates:
                continue
            elif sub_failure_index[2] >= 1 or sub_failure_index[1] >=1  or sub_failure_index[3] >=1:
                if i not in shear_failure:
                    properties[i] = [E1, 0.1 * E2, 0.1 * v12, 0.1 * G12]
                    shear_failure.append(i)
                    print("Degradation in laminate " + str(i))
                else:
                    properties[i] = [1e-10] * 4
                    failed_laminates.append(i)
                    print("Inter laminar failure " + str(i))
            elif sub_failure_index[0] >= 1:
                properties[i] = [1e-10] * 4
                failed_laminates.append(i)
                print("Failure in laminate " + str(i))
        #if len(failed_laminates) > 0 and fpf == False:
        if len(shear_failure) > 0 and fpf == False:   
            fpf = True
            fpf_Ns = (r)*np.sin(the)
            fpf_Ny = (r)*np.cos(the)
            strain_fpf_x_PK.append(strain_laminas[0][0])
            strain_fpf_y_PK.append(strain_laminas[0][1])
            strain_fpf_s_PK.append(strain_laminas[0][2])
            print("First ply failure in " + str(i))
            
        failed_laminates_f = set(failed_laminates)
        if len(failed_laminates_f) == len(theta):
            lpf = True
            lpf_Ns = r*np.sin(the)
            lpf_Ny = r*np.cos(the)
            strain_lpf_x_PK.append(strain_laminas[0][0])
            strain_lpf_y_PK.append(strain_laminas[0][1])
            strain_lpf_s_PK.append(strain_laminas[0][2])
            break
        r = r+1
        Force_per_length[2] = r*np.sin(the)
        Force_per_length[1] = r*np.cos(the)
    return fpf_Ns, fpf_Ny, lpf_Ns, lpf_Ny, strain_fpf_x_PK, strain_fpf_y_PK, strain_fpf_s_PK, strain_lpf_x_PK, strain_lpf_y_PK, strain_lpf_s_PK

def failure_envelope_maxStress(thickness, theta, properties, Force_per_length, Xt, Xc, Yt, Yc, S, r, the):
    failed_laminates = []
    shear_failure = []
    fpf = False
    lpf = False
    while fpf == False or lpf == False:
        A, B, D, z = ABD_matrix_differentMaterial(thickness, theta, properties)
        stress_laminas = stress_per_lamina_differentMaterials(Force_per_length, theta, A, B, D, properties, z)
        strain_laminas = strain_per_lamina_differentMaterials(Force_per_length, theta, A, B, D, properties, z)
        failure_index = max_stress_criteria_envelope(stress_laminas, Xt, Xc, Yt, Yc, S)
        for i, sub_failure_index in enumerate(failure_index):
            if i in failed_laminates:
                continue
            if any(value >= 1 for value in sub_failure_index):
                if sub_failure_index[2] >= 1 or sub_failure_index[1] >=1:
                    if i not in shear_failure:
                        properties[i] = [E1, 0.1 * E2, 0.1 * v12, 0.1 * G12]
                        shear_failure.append(i)
                        print("Shear degradation in laminate " + str(i))
                    else:
                        properties[i] = [1e-10] * 4
                        failed_laminates.append(i)
                if sub_failure_index[0] >= 1:
                    properties[i] = [1e-10] * 4
                    failed_laminates.append(i)
                    print("Failure in laminate " + str(i))
        #if len(failed_laminates) > 0 and fpf == False:
        if len(shear_failure) > 0 and fpf == False:    
            fpf = True
            fpf_Ns = (r)*np.sin(the)
            fpf_Ny = (r)*np.cos(the)
            strain_fpf_x_MS.append(strain_laminas[0][0])
            strain_fpf_y_MS.append(strain_laminas[0][1])
            strain_fpf_s_MS.append(strain_laminas[0][2])
            print("First ply failure in " + str(i))
            
        failed_laminates_f = set(failed_laminates)
        if len(failed_laminates_f) == len(theta):
            lpf = True
            lpf_Ns = r*np.sin(the)
            lpf_Ny = r*np.cos(the)
            strain_lpf_x_MS.append(strain_laminas[0][0])
            strain_lpf_y_MS.append(strain_laminas[0][1])
            strain_lpf_s_MS.append(strain_laminas[0][2])
            break
        r = r+1
        Force_per_length[2] = r*np.sin(the)
        Force_per_length[1] = r*np.cos(the)
    return fpf_Ns, fpf_Ny, lpf_Ns, lpf_Ny, strain_fpf_x_MS, strain_fpf_y_MS, strain_fpf_s_MS, strain_lpf_x_MS, strain_lpf_y_MS, strain_lpf_s_MS
#the1 = np.linspace(np.pi/8.5,np.pi/8,3)
the1 = np.arange(0,2*np.pi+np.pi/180,np.pi/180)
for the in the1:
    properties_o = [[E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12]]
    Nx = 0
    Ny = 0
    Ns = 0
    Mx = 0
    My = 0
    Ms = 0
    Force_per_length = np.array([[Nx],[Ny],[Ns],[Mx],[My],[Ms]])
    fpf_Ns, fpf_Ny, lpf_Ns, lpf_Ny, strain_fpf_x_PK, strain_fpf_y_PK, strain_fpf_s_PK, strain_lpf_x_PK, strain_lpf_y_PK, strain_lpf_s_PK = failure_envelope_puck(thickness, theta, properties_o, Force_per_length, 1932, 1480, 108, 220, 132.8, 0, the)
    print("Iteration done for " + str(the))
    fpf_Ns_plot_Puck.append(fpf_Ns/sum(thickness))
    fpf_Ny_plot_Puck.append(fpf_Ny/sum(thickness))
    lpf_Ns_plot_Puck.append(lpf_Ns/sum(thickness))
    lpf_Ny_plot_Puck.append(lpf_Ny/sum(thickness))

for the in the1:
    properties_o = [[E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12],
              [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12], [E1, E2, v12, G12]]
    Nx = 0
    Ny = 0
    Ns = 0
    Mx = 0
    My = 0
    Ms = 0
    Force_per_length = np.array([[Nx],[Ny],[Ns],[Mx],[My],[Ms]])
    fpf_Ns, fpf_Ny, lpf_Ns, lpf_Ny, strain_fpf_x_MS, strain_fpf_y_MS, strain_fpf_s_MS, strain_lpf_x_MS, strain_lpf_y_MS, strain_lpf_s_MS = failure_envelope_maxStress(thickness, theta, properties_o, Force_per_length, 1932, 1480, 108, 220, 132.8, 0, the)
    print("Iteration done for " + str(the))
    fpf_Ns_plot_MaxStress.append(fpf_Ns/sum(thickness))
    fpf_Ny_plot_MaxStress.append(fpf_Ny/sum(thickness))
    lpf_Ns_plot_MaxStress.append(lpf_Ns/sum(thickness))
    lpf_Ny_plot_MaxStress.append(lpf_Ny/sum(thickness))


plt.figure(figsize=(8, 6))  
plt.plot(fpf_Ny_plot_Puck, fpf_Ns_plot_Puck, marker='x',label = "First ply failure Puck")
plt.plot(lpf_Ny_plot_Puck, lpf_Ns_plot_Puck, marker='x',label = "Last ply failure Puck")
plt.plot(fpf_Ny_plot_MaxStress, fpf_Ns_plot_MaxStress, marker='.',label = "First ply failure MaxStres")
plt.plot(lpf_Ny_plot_MaxStress, lpf_Ns_plot_MaxStress, marker='.',label = "Last ply failure MaxStress")
plt.xlabel("Sigmay in MPa")
plt.ylabel("Sigmas in MPa")
plt.title("Max Stress vs Puck criteria")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))  
plt.plot(180/np.pi*the1,strain_fpf_x_MS,marker = '.',label = "strain_x for fpf")
plt.plot(180/np.pi*the1,strain_fpf_y_MS,marker = '.',label = "strain_y for fpf")
plt.plot(180/np.pi*the1,strain_fpf_s_MS,marker = '.',label = "strain_s for fpf")
#plt.plot(strain_fpf_y,strain_fpf_s, marker = 'x')
plt.xlabel("Iterations in Theta")
plt.ylabel("Strain")
plt.title("Max stress criteria")
plt.legend()
plt.grid(True)
plt.show()


plt.figure(figsize=(8, 6))
plt.plot( 180/np.pi*the1,strain_lpf_x_PK,marker = '.', label = "strain_x for lpf")
plt.plot( 180/np.pi*the1,strain_lpf_y_PK,marker = '.', label = "strain_y for lpf")
plt.plot( 180/np.pi*the1,strain_lpf_s_PK,marker = '.', label = "strain_s for lpf")
#plt.plot(strain_lpf_y,strain_lpf_s, marker = 'x')
plt.xlabel("Iterations in Theta")
plt.ylabel("Strain")
plt.title("Max stress criteria")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(8, 6))  
plt.plot(180/np.pi*the1,strain_fpf_x_MS,marker = 'x', label = "strain_x for fpf")
plt.plot(180/np.pi*the1,strain_fpf_y_MS,marker = 'x', label = "strain_y for fpf")
plt.plot(180/np.pi*the1,strain_fpf_s_MS,marker = 'x', label = "strain_s for fpf")
#plt.plot(strain_fpf_y,strain_fpf_s, marker = 'x')
plt.xlabel("Theta")
plt.ylabel("strain")
plt.title("Puck criteria")
plt.legend()
plt.grid(True)
plt.show()


plt.figure(figsize=(8, 6))
plt.plot( 180/np.pi*the1,strain_lpf_x_PK,marker = 'x',  label = "strain_x for lpf")
plt.plot( 180/np.pi*the1,strain_lpf_y_PK,marker = 'x',  label = "strain_y for lpf")
plt.plot( 180/np.pi*the1,strain_lpf_s_PK,marker = 'x',  label = "strain_s for lpf")
#plt.plot(strain_lpf_y,strain_lpf_s, marker = 'x')
plt.xlabel("Theta")
plt.ylabel("strain")
plt.title("Puck criteria")
plt.legend()
plt.grid(True)
plt.show()




