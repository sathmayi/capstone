import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#initial states, values of k
A0 = 10
B0 = 0
C0 = 0
k1 = 1
k2 = 0.5
y0 = [A0, B0, C0] #list of initial states

#creating the ODE system
def model(y,t):
  A, B, C = y
  dAdt = -k1 * A
  dBdt = k1 * A - k2 * B
  dCdt = k2 * B
  return [dAdt, dBdt, dCdt]

#creating a time array
t = np.linspace(0,10)

eq = odeint(model, y0, t) #solving ODE
A = eq[:, 0]
B = eq[:, 1]
C = eq[:, 2]

plt.plot(t, A, 'r-', linewidth=2, label='A')
plt.plot(t, B, 'b--', linewidth=2, label='B')
plt.plot(t, C, 'g:', linewidth=2, label='C')
plt.xlabel('time')
plt.ylabel('Concentration of A, B, C')
plt.legend()
plt.show()
