
import numpy as np
import matplotlib.pyplot as plt
a = 1  
x = np.linspace(0, a, 1000) 

v2 = np.array([124.033, 1])
v2_1 = np.array([24.8453,1])


v2_norm = v2 / np.linalg.norm(v2)
v2_1norm = v2_1 / np.linalg.norm(v2_1)
c1,c2 = v2_norm

print (c1,c2)

f_x = c1 * np.sqrt(2 / a) * np.sin(np.pi * x / a) + c2 * np.sqrt(2 / a) * np.sin(3 * np.pi * x / a)

psi_ground = np.sqrt(2 / a) * np.sin(np.pi * x / a)

plt.figure(figsize=(10, 6))
plt.plot(x, f_x, label="Trial Wavefunction $f(x)$", linewidth=2)
plt.plot(x, psi_ground, label="Unperturbed Ground State $\psi_1(x)$", linestyle='--', linewidth=2)
plt.axhline(0, color='black', linewidth=0.8, linestyle='dashed')

plt.title("Comparison of Trial Wavefunction $f(x)$ and Unperturbed Ground State", fontsize=14)
plt.xlabel("$x$ (position)", fontsize=12)
plt.ylabel("Wavefunction Amplitude", fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)

plt.show()
