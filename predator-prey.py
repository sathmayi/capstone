
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

y0 = [10, 1] #cannot be 0 as populations have to grow
#y0 = [gamma/delta, alpha/beta] to check steady state values 

#prey and predator in hundreds of units

t = np.linspace(0, 50, num=1000)
# we are not specifying units of time

#define constant parameters

alpha = 1.1
beta = 0.4
delta = 0.1
gamma = 0.4

params = [alpha, beta, delta, gamma] #put into an array

def model(variables, t, params):

    #define all things locally
    x = variables[0] #prey pop level
    y = variables[1] #predator pop level

    alpha = params[0]
    beta = params[1]
    delta = params[2]
    gamma = params[3]

    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y

    return([dxdt, dydt])

#get output

y = odeint(model, y0, t, args = (params,))
#creating y matrix, where a column for each species, and a row for each timepoint

f, (plt1, plt2) = plt.subplots(2) #creating two diff subplots

line1 = plt1.plot(t,y[:,0], color = "b")
#plot time points, and prey population at each timepoint
#saying we want every row in the first column to plot each timepoint
line2 = plt2.plot(t,y[:,1], color = "r")
#to get every row in second column

plt1.set_ylabel('Prey (hundereds)')
plt2.set_ylabel('Predators (hundereds)')
plt2.set_xlabel("Time")
plt.show()
