import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import random
from scipy.stats import norm


def stress(S1, S2, S3, S4, S5, S6):
    stress = np.array([
        [S1],
        [S2],
        [S3],
        [S4],
        [S5],
        [S6]
    ])
    return stress

def stress_inplane(S1,S2,S6):
    stress = np.array([
        [S1],
        [S2],
        [S6]
    ])
    return stress

def compliance_orthotropic(E1,E2,E3, v12, v23, v13, G23, G13, G12):
    compliance_matrix = np.array([
        [1/E1, -v12/E1, -v13/E1, 0, 0, 0],
        [-v12/E1, 1/E2, -v23/E2, 0, 0, 0],
        [-v13/E1, -v23/E2, 1/E3, 0, 0, 0],
        [0, 0, 0, 1/G23, 0, 0],
        [0, 0, 0, 0, 1/G13, 0],
        [0, 0, 0, 0, 0, 1/G12]])
    return compliance_matrix

def compliance_orthotropic_inplane(E1,E2,v12, G12):
    compliance_matrix = np.array([
        [1/E1, -v12/E1, 0],
        [-v12/E1, 1/E2, 0],
        [0, 0, 1/G12]])
    return compliance_matrix

def compliance_transverselyIsotropic(E1, E2, v12, v23, G12):
    compliance_matrix = np.array([
        [1/E1, -v12/E1, -v12/E1, 0, 0, 0],
        [-v12/E1, 1/E2, -v23/E2, 0, 0, 0],
        [-v12/E1, -v23/E2, 1/E2, 0, 0, 0],
        [0, 0, 0, 2*(1+v23)/E2, 0, 0],
        [0, 0, 0, 0, 1/G12, 0],
        [0, 0, 0, 0, 0, 1/G12]
    ])
    return compliance_matrix

def compliance_transverselyIsotropic_inplane(E1, E2, v12, G12):
    compliance_matrix = np.array([
        [1/E1, -v12/E1, 0],
        [-v12/E1, 1/E2, 0],
        [0, 0, 1/G12]
    ])
    return compliance_matrix

def compliance_isotropic(E,v):
    compliance_matrix = np.array([
        [1/E, -v/E, -v/E, 0, 0, 0],
        [-v/E, 1/E, -v/E, 0, 0, 0],
        [-v/E, -v/E, 1/E, 0, 0, 0],
        [0, 0, 0, 2*(1+v)/E, 0, 0],
        [0, 0, 0, 0, 2*(1+v)/E, 0],
        [0, 0, 0, 0, 0, 2*(1+v)/E]
    ])
    return compliance_matrix

def compliance_isotropic_inplane(E,v):
    compliance_matrix = np.array([
        [1/E, -v/E, 0],
        [-v/E, 1/E, 0],
        [0, 0, 2*(1+v)/E]
    ])
    return compliance_matrix

def compliance_transform(theta, compliance):
    theta = theta * np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    
    T_stress = np.array([
        [m**2,n**2,0,0,0,2*m*n],
        [n**2,m**2,0,0,0,-2*m*n],
        [0,0,1,0,0,0],
        [0,0,0,m,-n,0],
        [0,0,0,n,m,0],
        [-m*n,m*n,0,0,0,m**2-n**2]])

    T_strain = np.array([
        [m**2,n**2,0,0,0,m*n],
        [n**2,m**2,0,0,0,-m*n],
        [0,0,1,0,0,0],
        [0,0,0,m,-n,0],
        [0,0,0,n,m,0],
        [-2*m*n,2*m*n,0,0,0,m**2-n**2]])

    compliance_transformed = np.matmul(np.matmul(np.linalg.inv(T_strain), compliance), T_stress)

    return compliance_transformed

def compliance_transform_inplane(theta, compliance):
    theta = theta * np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    
    T_stress = np.array([
        [m**2,n**2,2*m*n],
        [n**2,m**2,-2*m*n],
        [-m*n,m*n,m**2-n**2]])

    T_strain = np.array([
        [m**2,n**2, m*n],
        [n**2,m**2,-m*n],
        [-2*m*n,2*m*n,m**2-n**2]])

    compliance_transformed = np.matmul(np.matmul(np.linalg.inv(T_strain), compliance), T_stress)

    return compliance_transformed

def stiffness_transform(theta, stiffness):
    theta = theta * np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    
    T_stress = np.array([
        [m**2,n**2,0,0,0,2*m*n],
        [n**2,m**2,0,0,0,-2*m*n],
        [0,0,1,0,0,0],
        [0,0,0,m,-n,0],
        [0,0,0,n,m,0],
        [-m*n,m*n,0,0,0,m**2-n**2]])

    T_strain = np.array([
        [m**2,n**2,0,0,0,m*n],
        [n**2,m**2,0,0,0,-m*n],
        [0,0,1,0,0,0],
        [0,0,0,m,-n,0],
        [0,0,0,n,m,0],
        [-2*m*n,2*m*n,0,0,0,m**2-n**2]])

    stiffness_transformed = np.matmul(np.matmul(np.linalg.inv(T_stress), stiffness), T_strain)

    return stiffness_transformed

def stiffness_transform_inplane(theta, stiffness):
    theta = theta * np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    
    T_stress = np.array([
        [m**2,n**2,2*m*n],
        [n**2,m**2,-2*m*n],
        [-m*n,m*n,m**2-n**2]])

    T_strain = np.array([
        [m**2,n**2, m*n],
        [n**2,m**2,-m*n],
        [-2*m*n,2*m*n,m**2-n**2]])

    stiffness_transformed = np.matmul(np.matmul(np.linalg.inv(T_stress), stiffness), T_strain)

    return stiffness_transformed


def stress_transform(theta, stress):
    theta = theta*np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    T_stress = np.array([
        [m**2,n**2,0,0,0,2*m*n],
        [n**2,m**2,0,0,0,-2*m*n],
        [0,0,1,0,0,0],
        [0,0,0,m,-n,0],
        [0,0,0,n,m,0],
        [-m*n,m*n,0,0,0,m**2-n**2]])
    stress_transformed = np.matmul(T_stress, stress)
    return stress_transformed

def stress_transform_inplane(theta, stress):
    theta = theta*np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    T_stress = np.array([
        [m**2,n**2,2*m*n],
        [n**2,m**2,-2*m*n],
        [-m*n,m*n,m**2-n**2]])
    stress_transformed = np.matmul(T_stress, stress)
    return stress_transformed

def strain_transform(theta, strain):
    theta = theta*np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    T_strain = np.array([
        [m**2,n**2,0,0,0,m*n],
        [n**2,m**2,0,0,0,-m*n],
        [0,0,1,0,0,0],
        [0,0,0,m,-n,0],
        [0,0,0,n,m,0],
        [-2*m*n,2*m*n,0,0,0,m**2-n**2]])
    strain_transformed = np.matmul(T_strain, strain)
    return strain_transformed

def strain_transform_inplane(theta, strain):
    theta = theta*np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    T_strain = np.array([
        [m**2,n**2, m*n],
        [n**2,m**2,-m*n],
        [-2*m*n,2*m*n,m**2-n**2]])
    strain_transformed = np.matmul(T_strain, strain)
    return strain_transformed

def error_percent(actual, foundout):
    error_percent = abs((foundout-actual)/actual) * 100
    return error_percent

def custom_plot(x1, x2, label_x1 ="X-axis", label_x2 = "Y-axis", title = "Title"):
    random_color = (random.random(), random.random(), random.random())
    plt.figure(figsize=(8, 6))  
    plt.plot(x1, x2, color=random_color,label=label_x2)
    plt.xlabel(label_x1)
    plt.ylabel(label_x2)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()
    
def ABD_matrix(thickness, theta, E1, E2, v12, G12):
    number_of_laminas = len(thickness)
    compliance_initial = compliance_orthotropic_inplane(E1, E2, v12, G12)
    z = np.zeros((len(thickness)+1))
    middle_index = len(thickness) // 2
    new_thickness = []
    if len(thickness) % 2 ==0:
        for i in range(number_of_laminas+1):
                if i<len(thickness)/2:
                    z[i] = - sum(thickness[i:middle_index])
                else:
                    z[i] = sum(thickness[middle_index:i])
    else:
        if len(thickness) % 2 !=0:
            for i in range(len(thickness)):
                if i != middle_index:
                    new_thickness.append(thickness[i])
                else:
                    new_thickness.append(thickness[i]/2)
                    new_thickness.append(0)
                    new_thickness.append(thickness[i]/2)
        middle_index = len(new_thickness) // 2
        for i in range(number_of_laminas+1):
                if i<len(thickness)/2:
                    z[i] = - sum(new_thickness[i:middle_index])
                elif i == len(thickness)/2:
                    continue
                else:
                    z[i] = sum(new_thickness[middle_index:i+2])
    z = np.array(z)
    Q = []

    for i in theta:
        compliance_temp = compliance_transform_inplane(i, compliance_initial)
        Q_temp = np.linalg.inv(compliance_temp)
        Q.append(Q_temp)

    Q = np.array(Q)

    A = np.zeros((3,3))
    B = np.zeros((3,3))
    D = np.zeros((3,3))

    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,number_of_laminas):
                # if k<number_of_laminas//2:
                #     A[i][j] += Q[k][i][j] * (- z[k] + z[k+1])
                #     B[i][j] += Q[k][i][j] * (- z[k]**2 + z[k+1]**2)/2
                #     D[i][j] += Q[k][i][j] * (- z[k]**3 + z[k+1]**3)/3
                # else:
                A[i][j] += Q[k][i][j] * (z[k+1] - z[k])
                B[i][j] += Q[k][i][j] * (z[k+1]**2 - z[k]**2)/2
                D[i][j] += Q[k][i][j] * (z[k+1]**3 - z[k]**3)/3
    return A,B,D,z

def ABD_matrix_differentMaterial(thickness, theta, properties):
    number_of_laminas = len(thickness)
    compliance_initial = []
    for i in properties:
        compliance_initial.append(compliance_orthotropic_inplane(i[0], i[1], i[2], i[3])) 
    z = np.zeros((len(thickness)+1))
    middle_index = len(thickness) // 2
    new_thickness = []
    if len(thickness) % 2 ==0:
        for i in range(number_of_laminas+1):
                if i<len(thickness)/2:
                    z[i] = - sum(thickness[i:middle_index])
                else:
                    z[i] = sum(thickness[middle_index:i])
    else:
        if len(thickness) % 2 !=0:
            for i in range(len(thickness)):
                if i != middle_index:
                    new_thickness.append(thickness[i])
                else:
                    new_thickness.append(thickness[i]/2)
                    new_thickness.append(0)
                    new_thickness.append(thickness[i]/2)
        middle_index = len(new_thickness) // 2
        for i in range(number_of_laminas+1):
                if i<len(thickness)/2:
                    z[i] = - sum(new_thickness[i:middle_index])
                elif i == len(thickness)/2:
                    continue
                else:
                    z[i] = sum(new_thickness[middle_index:i+2])
    z = np.array(z)
    
    Q = []

    for i in range(len(theta)):
        compliance_temp = compliance_transform_inplane(theta[i], compliance_initial[i])
        Q_temp = np.linalg.inv(compliance_temp)
        Q.append(Q_temp)

    Q = np.array(Q)

    A = np.zeros((3,3))
    B = np.zeros((3,3))
    D = np.zeros((3,3))

    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,number_of_laminas):
                # if k<number_of_laminas//2:
                #     A[i][j] += Q[k][i][j] * (- z[k] + z[k+1])
                #     B[i][j] += Q[k][i][j] * (- z[k]**2 + z[k+1]**2)/2
                #     D[i][j] += Q[k][i][j] * (- z[k]**3 + z[k+1]**3)/3
                # else:
                A[i][j] += Q[k][i][j] * (z[k+1] - z[k])
                B[i][j] += Q[k][i][j] * (z[k+1]**2 - z[k]**2)/2
                D[i][j] += Q[k][i][j] * (z[k+1]**3 - z[k]**3)/3
    return A,B,D,z

def laminate_properties(total_thickness, A):
    Q = A/total_thickness
    compliance = np.linalg.inv(Q)
    Ex = round(1/compliance[0][0],2)
    Ey = round(1/compliance[1][1],2)
    Gxy = round(1/compliance[2][2],2)
    vxy = round(-Ex * compliance[1][0],2)
    vyx = round(-Ey * compliance[0][1],2)
    # A = np.linalg.inv(A)
    # Ex = 1/(total_thickness*A[0][0])
    # Ey = 1/(total_thickness*A[1][1])
    # Gxy = 1/(total_thickness* A[2][2])
    # vxy = -A[0][1]/A[0][0]
    # vyx = -A[0][1]/A[1][1]
    return Ex,Ey,vxy,vyx,Gxy

def laminate_properties_bending(total_thickness, D):
    D = np.linalg.inv(D)
    Exb = 12/(total_thickness**3 * D[0][0])
    Eyb = 12/(total_thickness**3 * D[1][1])
    Gxyb = 12/(total_thickness**3 *D[2][2])
    vxyb = - D[0][1]/D[0][0]
    vyxb = - D[1][0]/D[1][1]
    return Exb,Eyb,vxyb,vyxb,Gxyb

def stress_per_lamina(Force_per_length, theta, A,B,D, E1, E2, v12, G12, z):
    compliance = compliance_orthotropic_inplane(E1, E2, v12, G12)
    Q = np.linalg.inv(compliance_transform_inplane(0,compliance))
    ABD = np.zeros((6,6))
    ABD[0:3,0:3] = A
    ABD[0:3,3:6] = B
    ABD[3:6,0:3] = B
    ABD[3:6,3:6] = D
    
    strains_curvature = np.linalg.inv(ABD) @ Force_per_length
    stress_laminas = []
    for i in range(len(theta)):
        stress_laminas.append(Q @ strain_transform_inplane(theta[i], (strains_curvature[0:3] + 0.5*(z[i+1] + z[i]) * strains_curvature[3:6])))
    return stress_laminas

def stress_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z):
    Q = []
    for i in properties:
        compliance = compliance_orthotropic_inplane(i[0], i[1], i[2], i[3])
        Q_temp = np.linalg.inv(compliance_transform_inplane(0,compliance))
        Q.append(Q_temp)
    ABD = np.zeros((6,6))
    ABD[0:3,0:3] = A
    ABD[0:3,3:6] = B
    ABD[3:6,0:3] = B
    ABD[3:6,3:6] = D
    strains_curvature = np.linalg.inv(ABD) @ Force_per_length
    stress_laminas = []
    for i in range(len(theta)):
        stress_laminas.append(Q[i] @ strain_transform_inplane(theta[i], (strains_curvature[0:3] + 0.5*(z[i+1] + z[i]) * strains_curvature[3:6])))
    return stress_laminas

def strain_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z):
    Q = []
    for i in properties:
        compliance = compliance_orthotropic_inplane(i[0], i[1], i[2], i[3])
        Q_temp = np.linalg.inv(compliance_transform_inplane(0,compliance))
        Q.append(Q_temp)
    ABD = np.zeros((6,6))
    ABD[0:3,0:3] = A
    ABD[0:3,3:6] = B
    ABD[3:6,0:3] = B
    ABD[3:6,3:6] = D
    strains_curvature = np.linalg.inv(ABD) @ Force_per_length
    strains_laminas = []
    for i in range(len(theta)):
        # result = Q[i] @ strain_transform_inplane(theta[i], (strains_curvature[0:3] + 0.5*(z[i+1] + z[i]) * strains_curvature[3:6]))
        # rounded_result = np.round(result)
        strains_laminas.append(strain_transform_inplane(theta[i], (strains_curvature[0:3] + 0.5*(z[i+1] + z[i]) * strains_curvature[3:6])))
        # stress_laminas.append(rounded_result)
    return strains_laminas

def max_stress_criteria(stress_laminas, Xt, Xc, Yt, Yc, S, theta, properties, thickness, Force_per_length):
    
    failed_laminates = []
    number_of_laminas = len(stress_laminas)
    while number_of_laminas > 0:
        failure_index = []
        # print(stress_laminas)
        for i in range(len(stress_laminas)):
            if stress_laminas[i][0]>=0:
                fi_1 = stress_laminas[i][0]/Xt
            else:
                fi_1 = stress_laminas[i][0]/-Xc
            if stress_laminas[i][1]>=0:
                fi_2 = stress_laminas[i][1]/Yt
            else:
                fi_2 = stress_laminas[i][1]/-Yc
            if stress_laminas[i][2]>=0:
                fi_3 = abs(stress_laminas[i][2]/S)
            else:
                fi_3 = abs(stress_laminas[i][2]/S)
            failure_index.append([fi_1, fi_2, fi_3])
        
        
        for i, sublist in enumerate(failure_index):
            if any(value >= 1 for value in sublist):
                failed_laminates.append(i)
                print(str(theta[i]) + " laminate has failed at " +str(Force_per_length[0]) )

        for i in failed_laminates:
            properties[i] = [1e-10,1e-10,1e-10,1e-10]

        number_of_laminas =  len(stress_laminas) - len(failed_laminates)
        Force_per_length[0] = Force_per_length[0] + 1
        

        A,B,D,z = ABD_matrix_differentMaterial(thickness, theta, properties)
        stress_laminas = stress_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z)

    return failed_laminates

def custom_plot_plies(x1, x2, y_ticks_positions, y_tick_labels, grid_lines, label_x1 ="X-axis", title = "Title"):
    plt.figure(figsize=(8, 6))  
    plt.plot(x1, x2,label=label_x1)
    plt.xlabel(label_x1)
    plt.title(title)
    plt.legend()
    
    plt.yticks(y_ticks_positions, y_tick_labels)
    plt.ylim(-3.125e-4, 3.125e-4)
    for line in grid_lines:
        plt.axhline(y=line, color='gray', linestyle='--', linewidth=0.5)
    plt.show()



def normal_distribution(mean, std, values):
    x = np.linspace(mean - 4 * std, mean + 4 * std, values)
    #x = np.random.normal(mean, std, values)
    cdf = norm.cdf(x, mean, std)
    return x, cdf

def generate_random_from_cdf(x_values, cdf_values, num_samples):
    random_numbers = np.random.rand(num_samples)
    inverse_cdf = np.interp(random_numbers, cdf_values, x_values)    
    return inverse_cdf

def Puck(stress_laminas, Xt, Xc, Yt, Yc, S):
    p12t = 0.3 #pt | ||
    p12c = 0.25 #pc | ||
    p23c = 0.25 #pc | |
    FI = np.zeros([len(stress_laminas), 4])
    v21 = 0.31 * 8.5/145.3
    v21f = 0.2
    m_sigma_f = 1.1
    E1 = 145.3
    E1f = 230
    Ra = Yc/(2*(1+p23c))
    S12c = S * np.sqrt(1+2*p23c) * (np.sqrt(1+2*p12c*Yc/S) -1)
    failed = 0
    for i, stress in enumerate(stress_laminas):
        puck_stress = stress[0]-(v21-v21f* m_sigma_f * E1/E1f)*stress[1]
        if puck_stress[0]>=0:
            FI[i][0] = puck_stress[0]/Xt
        elif puck_stress[0] <0:
            FI[i][0] = abs(puck_stress[0])/Xc
        #Mode A
        if stress[1] > 0:
            FI[i][1] = (np.sqrt((stress[2]/S)**2 + (1 - p12t*Yt/S)**2 * (stress[1]/Yt)**2) + p12t*stress[1]/S)
        #IFF Mode B
        elif stress[1]<0 and abs(stress[1]/stress[2]) <= Ra/abs(S12c) and 0 <= abs(stress[1]/stress[2]):
            FI[i][2] = 1/S*(np.sqrt(stress[2]**2 + (p12c*stress[1])**2) + p12c*stress[1])
            #IFF Mode C
            # if stress[1] >= -1:
            #     FI[i][3] = 0
        elif stress[1]<0 and  abs(stress[2]/stress[1]) <= abs(S12c)/Ra and 0 <= abs(stress[2]/stress[1]):
            FI[i][3] = ((stress[2]/(2*(1 + p23c)*S))**2 + (stress[1]/Yc)**2)*(Yc)/(abs(stress[1]))
    for i, sublist in enumerate(FI):
        if any(value >= 1 for value in sublist):
            failed = 1
    return failed

def max_stress_criteria(stress_laminas, Xt, Xc, Yt, Yc, S):
    failure_index = np.zeros([len(stress_laminas), 3])
    failed = 0
    for i, stress in enumerate(stress_laminas):
        failure_index[i][0] = stress[0] / Xt if stress[0] >= 0 else abs(stress[0]) / Xc
        failure_index[i][1] = stress[1] / Yt if stress[1] >= 0 else abs(stress[1]) / Yc
        failure_index[i][2] = abs(stress[2]) / S
    for i, sublist in enumerate(failure_index):
        if any(value >= 1 for value in sublist):
            failed = 1
    return failed

def selecting_samples(value,order):
    condition = True
    E1_x, E1_cdf = normal_distribution(145.3e3, 3.28e3, int(1e6))
    E2_x, E2_cdf = normal_distribution(8.5e3, 1.28e3, int(1e6))
    v12_x, v12_cdf = normal_distribution(0.31, 0.018, int(1e6))
    G12_x, G12_cdf = normal_distribution(4.58e3, 0.83e3, int(1e6))
    Xt_x, Xt_cdf = normal_distribution(1932, 128.3, int(1e6))
    Xc_x, Xc_cdf = normal_distribution(1480, 128.3, int(1e6))
    Yt_x, Yt_cdf = normal_distribution(108, 8.2, int(1e6))
    Yc_x, Yc_cdf = normal_distribution(220, 8.2, int(1e6))
    S_x, S_cdf = normal_distribution(132.8, 6.21, int(1e6))

    while condition:
        E1_sampled = generate_random_from_cdf(E1_x, E1_cdf, value)
        E2_sampled = generate_random_from_cdf( E2_x, E2_cdf, value)
        v12_sampled = generate_random_from_cdf( v12_x, v12_cdf, value)
        G12_sampled = generate_random_from_cdf( G12_x, G12_cdf, value)
        Xt_sampled = generate_random_from_cdf( Xt_x, Xt_cdf, value)
        Xc_sampled = generate_random_from_cdf( Xc_x, Xc_cdf, value)
        Yt_sampled = generate_random_from_cdf(Yt_x, Yt_cdf, value)
        Yc_sampled = generate_random_from_cdf(Yc_x, Yc_cdf, value)
        S_sampled = generate_random_from_cdf(S_x, S_cdf, value)
        E1_obtained_mean = np.mean(E1_sampled)
        E2_obtained_mean = np.mean(E2_sampled)
        v12_obtained_mean = np.mean(v12_sampled)
        G12_obtained_mean = np.mean(G12_sampled)
        Xt_obtained_mean = np.mean(Xt_sampled)
        Xc_obtained_mean = np.mean(Xc_sampled)
        Yt_obtained_mean = np.mean(Yt_sampled)
        Yc_obtained_mean = np.mean(Yc_sampled)
        S_obtained_mean = np.mean(S_sampled)
        if abs(Xc_obtained_mean - 1480) >= 128.3/order and abs(Yc_obtained_mean - 220) >= 8.2/order and abs(E1_obtained_mean - 145.3) >= 3.28/order and abs(E2_obtained_mean - 8.5) >= 1.28/order and abs(v12_obtained_mean - 0.31) >= 0.018/order and abs(G12_obtained_mean - 4.58) >= 0.83/order and abs(Xt_obtained_mean - 1932) >= 128.3/order and abs(Yt_obtained_mean - 145.3) >= 8.2/order and abs(S_obtained_mean - 132.8) >= 6.21/order:
            # value += 1000
            continue
        else:
            condition = False  
            print(value)

    Sampled_properties = np.array([E1_sampled, E2_sampled, v12_sampled, G12_sampled, Xt_sampled, Xc_sampled, Yt_sampled, Yc_sampled, S_sampled])
    return Sampled_properties


def max_stress_criteria_envelope(stress_laminas, Xt, Xc, Yt, Yc, S):
    failure_index = np.zeros([len(stress_laminas), 3])
    for i, stress in enumerate(stress_laminas):
        failure_index[i][0] = stress[0] / Xt if stress[0] >= 0 else abs(stress[0]) / Xc
        failure_index[i][1] = stress[1] / Yt if stress[1] >= 0 else abs(stress[1]) / Yc
        failure_index[i][2] = abs(stress[2]) / S
    return failure_index

def Puck_envelope(stress_laminas, Xt, Xc, Yt, Yc, S):
    p12t = 0.3 #pt | ||
    p12c = 0.25 #pc | ||
    p23c = 0.25 #pc | |
    FI = np.zeros([len(stress_laminas), 4])
    v21 = 0.31 * 8.5/145.3
    v21f = 0.2
    m_sigma_f = 1.1
    E1 = 145.3e3
    E1f = 230e3
    Ra = Yc/(2*(1+p23c))
    S12c = S * np.sqrt(1+2*p23c)
    for i, stress in enumerate(stress_laminas):
        puck_stress = stress[0]-(v21-v21f* m_sigma_f * E1/E1f)*stress[1]
        if puck_stress[0]>=0:
            FI[i][0] = puck_stress[0]/Xt
        elif puck_stress[0] <0:
            FI[i][0] = abs(puck_stress[0])/Xc
        #Mode A
        if stress[1] > 0:
            FI[i][1] = (np.sqrt((stress[2]/S)**2 + (1 - p12t*Yt/S)**2 * (stress[1]/Yt)**2) + p12t*stress[1]/S)
        #IFF Mode B
        elif stress[1]<0 and abs(stress[1]/stress[2]) <= Ra/abs(S12c) and 0 <= abs(stress[1]/stress[2]):
            FI[i][2] = 1/S*(np.sqrt(stress[2]**2 + (p12c*stress[1])**2) + p12c*stress[1])
        #IFF Mode C
        elif stress[1]<0 and  abs(stress[2]/stress[1]) <= abs(S12c)/Ra and 0 <= abs(stress[2]/stress[1]):
            FI[i][3] = ((stress[2]/(2*(1 + p23c)*S))**2 + (stress[1]/Yc)**2)*(Yc)/(abs(stress[1]))
    return FI



