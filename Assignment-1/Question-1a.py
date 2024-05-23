from functions import *

E1 = 145.3e3
E2 = 8.5e3
v12 = 0.31
G12 = 4.58e3

theta_varied = np.linspace(0,180,90)

n_varied = np.linspace(1,20,20)



#For number of symmetries
Ex = []
Ey = []
vxy = []
vyx = []
Gxy = []
Exb = []
Eyb = []
vxyb = []
vyxb = []
Gxyb = []
for i in n_varied:
    theta = [15,0,-0,75,75,75,75,-0,0,15]
    theta = theta * int(i)
    thickness = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]
    thickness = thickness * int(i)
    A,B,D,z = ABD_matrix(thickness, theta, E1, E2, v12, G12)
    Ex1,Ey1,vxy1,vyx1,Gxy1 = laminate_properties(sum(thickness),A)
    Ex1b,Ey1b,vxy1b,vyx1b,Gxy1b =laminate_properties_bending(sum(thickness),D)
    Ex.append(Ex1/1e3)
    Ey.append(Ey1/1e3)
    vxy.append(vxy1)
    vyx.append(vyx1)
    Gxy.append(Gxy1/1e3)
    Exb.append(Ex1b/1e3)
    Eyb.append(Ey1b/1e3)
    vxyb.append(vxy1b)
    vyxb.append(vyx1b)
    Gxyb.append(Gxy1b/1e3)

# Create a figure and set of subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 8))

# Plot each graph on its respective subplot
axs[0, 0].plot(n_varied, Ex)
axs[0, 0].set_title("In-plane Youngs modulus in x-direction")
axs[0, 0].set_xlabel("Number of Repetitions")
axs[0, 0].set_ylabel("E in x-direction (GPa)")

axs[0, 1].plot(n_varied, Ey)
axs[0, 1].set_title("In-plane Youngs modulus in y-direction")
axs[0, 1].set_xlabel("Number of Repetitions")
axs[0, 1].set_ylabel("E in y-direction (GPa)")

axs[0, 2].plot(n_varied, Gxy)
axs[0, 2].set_title("In-plane Rigidity modulus")
axs[0, 2].set_xlabel("Number of Repetitions")
axs[0, 2].set_ylabel( "Rigidity Modulus (GPa)")

axs[1, 0].plot(n_varied, vxy)
axs[1, 0].set_title("In-plane Possion's ratio in xy-direction")
axs[1, 0].set_xlabel("Number of Repetitions")
axs[1, 0].set_ylabel("Poisson Ratio vxy")

axs[1, 1].plot(n_varied, vyx)
axs[1, 1].set_title("In-plane Poisson's ratio in yx-direction")
axs[1, 1].set_xlabel("Number of Repetitions")
axs[1, 1].set_ylabel("Poisson Ratio vyx")

axs[1, 2].axis('off')

# Adjust layout
fig.suptitle("Variation of in-plane properties with number of symmetry repetitions")
plt.tight_layout()

# Show the plot
plt.show()


# Create a figure and set of subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 8))

# Plot each graph on its respective subplot
axs[0, 0].plot(n_varied, Exb)
axs[0, 0].set_title("Flexural Youngs modulus in x-direction")
axs[0, 0].set_xlabel("Number of Repetitions")
axs[0, 0].set_ylabel("E in x-direction (GPa)")

axs[0, 1].plot(n_varied, Eyb)
axs[0, 1].set_title("Flexural Youngs modulus in y-direction")
axs[0, 1].set_xlabel("Number of Repetitions")
axs[0, 1].set_ylabel("E in y-direction (GPa)")

axs[0, 2].plot(n_varied, Gxyb)
axs[0, 2].set_title("Flexural Rigidity modulus")
axs[0, 2].set_xlabel("Number of Repetitions")
axs[0, 2].set_ylabel( "Rigidity Modulus (GPa)")

axs[1, 0].plot(n_varied, vxyb)
axs[1, 0].set_title("Flexural Possion's ratio in xy-direction")
axs[1, 0].set_xlabel("Number of Repetitions")
axs[1, 0].set_ylabel("Poisson Ratio vxy")

axs[1, 1].plot(n_varied, vyxb)
axs[1, 1].set_title("Flexural Poisson's ratio in yx-direction")
axs[1, 1].set_xlabel("Number of Repetitions")
axs[1, 1].set_ylabel("Poisson Ratio vyx")

axs[1, 2].axis('off')

# Adjust layout
fig.suptitle("Variation of flexural properties with number of symmetry repetitions")
plt.tight_layout()

# Show the plot
plt.show()


#For thetas

Ex = []
Ey = []
vxy = []
vyx = []
Gxy = []
Exb = []
Eyb = []
vxyb = []
vyxb = []
Gxyb = []
for j in theta_varied:
    theta = np.array([15,j,-j,75,75,75,75,-j,j,15])
    thickness = np.array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125])
    A,B,D,z = ABD_matrix(thickness, theta, E1, E2, v12, G12)
    Ex1,Ey1,vxy1,vyx1,Gxy1 = laminate_properties(sum(thickness),A)
    Ex1b,Ey1b,vxy1b,vyx1b,Gxy1b =laminate_properties_bending(sum(thickness),D)
    Ex.append(Ex1/1e3)
    Ey.append(Ey1/1e3)
    vxy.append(vxy1)
    vyx.append(vyx1)
    Gxy.append(Gxy1/1e3)
    Exb.append(Ex1b/1e3)
    Eyb.append(Ey1b/1e3)
    vxyb.append(vxy1b)
    vyxb.append(vyx1b)
    Gxyb.append(Gxy1b/1e3)


# Create a figure and set of subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 8))

# Plot each graph on its respective subplot
axs[0, 0].plot(theta_varied, Ex)
axs[0, 0].set_title("In-plane Youngs modulus in x-direction")
axs[0, 0].set_xlabel("Theta in degrees")
axs[0, 0].set_ylabel("E in x-direction (GPa)")

axs[0, 1].plot(theta_varied, Ey)
axs[0, 1].set_title("In-plane Youngs modulus in y-direction")
axs[0, 1].set_xlabel("Theta in degrees")
axs[0, 1].set_ylabel("E in y-direction (GPa)")

axs[0, 2].plot(theta_varied, Gxy)
axs[0, 2].set_title("In-plane Rigidity modulus")
axs[0, 2].set_xlabel("Theta in degrees")
axs[0, 2].set_ylabel( "Rigidity Modulus (GPa)")

axs[1, 0].plot(theta_varied, vxy)
axs[1, 0].set_title("In-plane Possion's ratio in xy-direction")
axs[1, 0].set_xlabel("Theta in degrees")
axs[1, 0].set_ylabel("Poisson Ratio vxy")

axs[1, 1].plot(theta_varied, vyx)
axs[1, 1].set_title("In-plane Poisson's ratio in yx-direction")
axs[1, 1].set_xlabel("Theta in degrees")
axs[1, 1].set_ylabel("Poisson Ratio vyx")

axs[1, 2].axis('off')

# Adjust layout
fig.suptitle("Variation of in-plane properties with Theta")
plt.tight_layout()

# Show the plot
plt.show()


# Create a figure and set of subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 8))

# Plot each graph on its respective subplot
axs[0, 0].plot(theta_varied, Exb)
axs[0, 0].set_title("Flexural Youngs modulus in x-direction")
axs[0, 0].set_xlabel("Theta in degrees")
axs[0, 0].set_ylabel("E in x-direction (GPa)")

axs[0, 1].plot(theta_varied, Eyb)
axs[0, 1].set_title("Flexural Youngs modulus in y-direction")
axs[0, 1].set_xlabel("Theta in degrees")
axs[0, 1].set_ylabel("E in y-direction (GPa)")

axs[0, 2].plot(theta_varied, Gxyb)
axs[0, 2].set_title("Flexural Rigidity modulus")
axs[0, 2].set_xlabel("Theta in degrees")
axs[0, 2].set_ylabel( "Rigidity Modulus (GPa)")

axs[1, 0].plot(theta_varied, vxyb)
axs[1, 0].set_title("Flexural Possion's ratio in xy-direction")
axs[1, 0].set_xlabel("Theta in degrees")
axs[1, 0].set_ylabel("Poisson Ratio vxy")

axs[1, 1].plot(theta_varied, vyxb)
axs[1, 1].set_title("Flexural Poisson's ratio in yx-direction")
axs[1, 1].set_xlabel("Theta in degrees")
axs[1, 1].set_ylabel("Poisson Ratio vyx")

axs[1, 2].axis('off')

# Adjust layout
fig.suptitle("Variation of flexural properties with Theta")
plt.tight_layout()

# Show the plot
plt.show()

#Variation of both theta and n-varied

for i in n_varied:
    Ex = []
    Ey = []
    vxy = []
    vyx = []
    Gxy = []
    Exb = []
    Eyb = []
    vxyb = []
    vyxb = []
    Gxyb = []
    for j in theta_varied:
        theta = [15,j,-j,75,75,75,75,-j,j,15]
        theta = theta * int(i)
        thickness = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]
        thickness = thickness * int(i)
        A,B,D,z = ABD_matrix(thickness, theta, E1, E2, v12, G12)
        Ex1,Ey1,vxy1,vyx1,Gxy1 = laminate_properties(sum(thickness),A)
        Ex1b,Ey1b,vxy1b,vyx1b,Gxy1b =laminate_properties_bending(sum(thickness),D)
        Ex.append(Ex1/1e3)
        Ey.append(Ey1/1e3)
        vxy.append(vxy1)
        vyx.append(vyx1)
        Gxy.append(Gxy1/1e3)
        Exb.append(Ex1b/1e3)
        Eyb.append(Ey1b/1e3)
        vxyb.append(vxy1b)
        vyxb.append(vyx1b)
        Gxyb.append(Gxy1b/1e3)

    # Create a figure and set of subplots
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))

    # Plot each graph on its respective subplot
    axs[0, 0].plot(theta_varied, Ex)
    axs[0, 0].set_title("In-plane Youngs modulus in x-direction")
    axs[0, 0].set_xlabel("Theta (°)")
    axs[0, 0].set_ylabel("E in x-direction (GPa)")

    axs[0, 1].plot(theta_varied, Ey)
    axs[0, 1].set_title("In-plane Youngs modulus in y-direction")
    axs[0, 1].set_xlabel("Theta (°)")
    axs[0, 1].set_ylabel("E in y-direction (GPa)")

    axs[0, 2].plot(theta_varied, Gxy)
    axs[0, 2].set_title("In-plane Rigidity modulus")
    axs[0, 2].set_xlabel("Theta (°)")
    axs[0, 2].set_ylabel( "Rigidity Modulus (GPa)")

    axs[1, 0].plot(theta_varied, vxy)
    axs[1, 0].set_title("In-plane Possion's ratio in xy-direction")
    axs[1, 0].set_xlabel("Theta (°)")
    axs[1, 0].set_ylabel("Poisson Ratio vxy")

    axs[1, 1].plot(theta_varied, vyx)
    axs[1, 1].set_title("In-plane Poisson's ratio in yx-direction")
    axs[1, 1].set_xlabel("Theta (°)")
    axs[1, 1].set_ylabel("Poisson Ratio vyx")

    axs[1, 2].axis('off')

    # Adjust layout
    fig.suptitle("Variation of in-plane properties with Theta " + str(i) +" repetitons")
    plt.tight_layout()

    # Show the plot
    plt.show()


    # Create a figure and set of subplots
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))

    # Plot each graph on its respective subplot
    axs[0, 0].plot(theta_varied, Exb)
    axs[0, 0].set_title("Flexural Youngs modulus in x-direction")
    axs[0, 0].set_xlabel("Theta (°)")
    axs[0, 0].set_ylabel("E in x-direction (GPa)")

    axs[0, 1].plot(theta_varied, Eyb)
    axs[0, 1].set_title("Flexural Youngs modulus in y-direction")
    axs[0, 1].set_xlabel("Theta (°)")
    axs[0, 1].set_ylabel("E in y-direction (GPa)")

    axs[0, 2].plot(theta_varied, Gxyb)
    axs[0, 2].set_title("Flexural Rigidity modulus")
    axs[0, 2].set_xlabel("Theta (°)")
    axs[0, 2].set_ylabel( "Rigidity Modulus (GPa)")

    axs[1, 0].plot(theta_varied, vxyb)
    axs[1, 0].set_title("Flexural Possion's ratio in xy-direction")
    axs[1, 0].set_xlabel("Theta (°)")
    axs[1, 0].set_ylabel("Poisson Ratio vxy")

    axs[1, 1].plot(theta_varied, vyxb)
    axs[1, 1].set_title("Flexural Poisson's ratio in yx-direction")
    axs[1, 1].set_xlabel("Theta (°)")
    axs[1, 1].set_ylabel("Poisson Ratio vyx")

    axs[1, 2].axis('off')

    # Adjust layout
    fig.suptitle("Variation of flexural properties with Theta " + str(i) +" repetitons")
    plt.tight_layout()

    # Show the plot
    plt.show()

#     custom_plot(theta_varied, Ex, "Theta (°)", "E in x-direction (GPa)", "Theta vs Youngs modulus in x-direction " + str(i) + " symmetry")
#     custom_plot(theta_varied, Ey, "Theta (°)", "E in y-direction (GPa)", "Theta vs Youngs modulus in y-direction")
#     custom_plot(theta_varied, Gxy, "Theta (°)", "Rigidity Modulus (GPa)", "Theta vs Rigidity modulus")
#     custom_plot(theta_varied, vxy, "Theta (°)", "Poisson Ratio vxy", "Theta vs Poisson ratio")
#     custom_plot(theta_varied, vyx, "Theta (°)", "Poisson Ratio vyx", "Theta vs Poisson ratio")
#     custom_plot(theta_varied, Exb, "Theta (°)", "E in x-direction (GPa)", "Theta vs Youngs modulus in x-direction " + str(i) + " symmetry")
#     custom_plot(theta_varied, Eyb, "Theta (°)", "E in y-direction (GPa)", "Theta vs Youngs modulus in y-direction")
#     custom_plot(theta_varied, Gxyb, "Theta (°)", "Rigidity Modulus (GPa)", "Theta vs Rigidity modulus")
#     custom_plot(theta_varied, vxyb, "Theta (°)", "Poisson Ratio vxy", "Theta vs Poisson ratio")
#     custom_plot(theta_varied, vyxb, "Theta (°)", "Poisson Ratio vyx", "Theta vs Poisson ratio")