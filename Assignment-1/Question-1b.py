from functions import *
#Inputs
E1 = 145.3e9
E2 = 8.5e9
v12 = 0.31
G12 = 4.58e9
theta = [0,0,90,30,90]
thickness = [0.125e-3,0.125e-3,0.125e-3,0.125e-3,0.125e-3]
properties = [[E1,E2,v12,G12],[E1,E2,v12,G12],[E1,E2,v12,G12],[E1,E2,v12,G12],[E1,E2,v12,G12]]
Nx = 0.2e2
Ny = 1.8e4
Ns = 0
Mx = 18e3
My = 0
Ms = 0
Force_per_length = np.array([[Nx],[Ny],[Ns],[Mx],[My],[Ms]])

#computations
A,B,D,z = ABD_matrix_differentMaterial(thickness,theta,properties)
print(z)

strains_in_laminas = strain_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z)
stress_in_laminas = stress_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z)

#post-processing 

print("The stress in lamina 0 degrees is "+ str(stress_in_laminas[0]/1e9)+" GPa")
print("The stress in lamina 0 degrees is "+ str(stress_in_laminas[1]/1e9)+" GPa")
print("The stress in lamina 90 degrees is "+ str(stress_in_laminas[2]/1e9)+" GPa")
print("The stress in lamina 30 degrees is "+ str(stress_in_laminas[3]/1e9)+" GPa")
print("The stress in lamina 90 degrees is "+ str(stress_in_laminas[4]/1e9)+" GPa")

print("The strains in lamina 0 degrees is "+ str(strains_in_laminas[0]))
print("The strains in lamina 0 degrees is "+ str(strains_in_laminas[1]))
print("The strains in lamina 90 degrees is "+ str(strains_in_laminas[2]))
print("The strains in lamina 30 degrees is "+ str(strains_in_laminas[3]))
print("The strains in lamina 90 degrees is "+ str(strains_in_laminas[4]))
# print(strains_in_laminas)


z1 = np.arange(-3.125e-4,3.125e-4,6.25e-6)
#strains_in_laminas, stress_in_laminas = stress_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z1)
stress_1 = []
stress_2 = []
stress_12 = []
strain_1 = []
strain_2 = []
strain_12 = []
for i in range(len(strains_in_laminas)):
    stress_1.append(stress_in_laminas[i][0]/1e9)
    stress_2.append(stress_in_laminas[i][1]/1e9)
    stress_12.append(stress_in_laminas[i][2]/1e9)
    strain_1.append(strains_in_laminas[i][0])
    strain_2.append(strains_in_laminas[i][1])
    strain_12.append(strains_in_laminas[i][2])
    
stress_1 = np.tile(np.array(stress_1),20)
stress_1 = stress_1.reshape(-1, len(stress_1)).flatten()
stress_2 = np.tile(np.array(stress_2), 20)
stress_2 = stress_2.reshape(-1, len(stress_2)).flatten()
stress_12 = np.tile(np.array(stress_12), 20)
stress_12 = stress_12.reshape(-1, len(stress_12)).flatten()
strain_1 = np.tile(np.array(strain_1), 20)
strain_1 = strain_1.reshape(-1, len(strain_1)).flatten()
strain_2 = np.tile(np.array(strain_2), 20)
strain_2 = strain_2.reshape(-1, len(strain_2)).flatten()
strain_12 = np.tile(np.array(strain_12), 20)
strain_12 = strain_12.reshape(-1, len(strain_12)).flatten()

y_ticks_positions = [-2.5e-4, -1.25e-4, 0, 1.25e-4, 2.5e-4]
y_tick_labels = ["Ply 1: 0°", "Ply 2: 0°", "Ply 3: 90°", "Ply 4: 30°", "Ply 5: 90°"]
grid_lines = [-3.125e-4, -1.875e-4, -6.25e-5, 6.25e-5, 1.875e-4, 3.125e-4]
print(z1)
custom_plot_plies(stress_1, z1, y_ticks_positions, y_tick_labels, grid_lines,"Stress (GPa)", "Stress in principal direction 1 along different laminates")
custom_plot_plies(stress_2, z1, y_ticks_positions, y_tick_labels, grid_lines,"Stress (GPa)", "Stress in principal direction 2 along different laminates")
custom_plot_plies(stress_12, z1, y_ticks_positions, y_tick_labels,grid_lines, "Stress (GPa)", "Stress in principal direction 12 along different laminates")
custom_plot_plies(strain_1, z1, y_ticks_positions, y_tick_labels, grid_lines,"Strain (no units)",  "Strain in principal direction 1 along different laminates")
custom_plot_plies(strain_2, z1, y_ticks_positions, y_tick_labels, grid_lines,"Strain (no units)",  "Strain in principal direction 2 along different laminates")
custom_plot_plies(strain_12, z1, y_ticks_positions, y_tick_labels,grid_lines, "Strain (no units)", "Strain in principal direction 12 along different laminates")