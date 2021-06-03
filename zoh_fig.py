import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

lp_filter = sig.butter(1, 5, "low", analog=True)
t_cont = np.linspace(0, 1, 10000)

f_t = sig.step(lp_filter, T=t_cont)

x_n = f_t[0]

plt.figure()
plt.plot(f_t[0], f_t[1])

plt.show()