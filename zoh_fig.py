import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

lp_filter = sig.butter(1, 5, "low", analog=True)
t_cont = np.linspace(0, 1, 10000)

f_t = sig.step(lp_filter, T=t_cont)

x_n = np.zeros(11)
t = np.zeros(11)

for x in np.arange(x_n.size-1):
    x_n[x] = f_t[1][x*1000]
    t[x] = f_t[0][x*1000]
t[-1] = 1
x_n[-1] = f_t[1][-1]

x_2 = x_n[0:-1]
t_2 = t[1:]

t_2 = t_2*10
f_t = f_t
t = t *10
#x_n = (f_t[1][0], f_t[1][1000], f_t[1][2000], f_t[1][3000], f_t[1][4000], f_t[1][5000], f_t[1][6000, f_t[1][7000], f_t[1][8000], f_t[1][9000], f_t[1][10000])

print(x_n[:-1])

plt.figure()
plt.title("Zero-order hold")
plt.plot(f_t[0]*10, f_t[1], label="$E_k(t)$")
plt.stem(t, x_n, "r--", markerfmt="kx", basefmt="k", label="sample point")
#plt.vlines(t_2,  x_2, x_n[:-1])
#plt.stem(t_2, x_2, "r")
plt.hlines(x_n[:-1], t[:-1], t_2, "k", linestyles="dotted", label=r"$E_k[k]$ (zoh)")
plt.xlabel("Time index [k]")
plt.ylabel("Amplitude")
plt.legend()
plt.savefig("zoh_fig")
plt.show()
