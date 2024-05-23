import numpy as np
from functions import *
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy import stats
import random

theta = [0,90,+45,-45,-45,+45,90,0,0,90,+45,-45,-45,+45,90,0]
thickness = [0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]
#Force = 850
Force = 1200
Force_per_length = np.array([[Force*np.cos(30*np.pi/180)],[Force*np.sin(30*np.pi/180)],[0],[0],[0],[0]]) 

def calculate_ci(data, confidence=0.95):
    data = np.array(data)
    mean = np.mean(data)
    n = len(data)
    stderr = np.std(data, ddof=1) / np.sqrt(n)
    margin = stderr * np.abs(stats.t.ppf((1 - confidence) / 2, n - 1))
    lower_bound = mean - margin
    upper_bound = mean + margin
    return mean, lower_bound, upper_bound

def convergence_reliability(values, order):
    Sampled_properties = selecting_samples(values,order)
    number_of_fails = []
    i = 0
    while i < Sampled_properties.shape[1]:
        properties = [[Sampled_properties[0][i],Sampled_properties[1][i], Sampled_properties[2][i], Sampled_properties[3][i]]]*16
        A,B,D,z = ABD_matrix_differentMaterial(thickness, theta, properties)
        stress_laminas = stress_per_lamina_differentMaterials(Force_per_length, theta, A,B,D, properties,z)
        #number_of_fails.append(max_stress_criteria_reliability(stress_laminas, Xt, 1480, Yt, 220, S))
        number_of_fails.append(Puck(stress_laminas, Sampled_properties[4][i], Sampled_properties[5][i], Sampled_properties[6][i], Sampled_properties[7][i], Sampled_properties[8][i]))
        i = i +1
    print("Number of failures " + str(number_of_fails.count(1)))
    print("Total number of runs " + str(len(number_of_fails)))
    prob_failure = number_of_fails.count(1)/len(number_of_fails) * 100
    print(prob_failure)
    
    return prob_failure
x_iterations = []
y_values = []
#while upper_bound - lower_bound >=0:
values = 1000
while values <= 1e5:
    upper = 100
    lower = 0
    mean = 100
    count  =  0
    prob_failure = []
    iterations = 0
    while (upper-lower)/mean * 100 > 5:
        count = 0
        while count <=5:
            prob_failure.append(convergence_reliability(values,10))
            count += 1
            iterations += 1
        confidence_level = 0.95
        mean, lower, upper = calculate_ci(prob_failure, confidence_level)
        print(f"Mean ({mean}) [{lower}, {upper}] : Confidence interval ({confidence_level * 100}%)")
    y_values.append(values) #take from initial thing
    if values >= 10000:
        values = values + 10000
    else:
        values += 1000
    x_iterations.append(iterations)


custom_plot(y_values, x_iterations , "No. of Samples taken", "No. of iterations required", "Iterations vs Samples taken for error of 15%")





# Example usage:



#Simple. Write a code for checking multiple times at same value. Count how many times done to get confidence level withing +-1 %