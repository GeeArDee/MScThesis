import numpy as np
import matplotlib.pyplot as plt

exp_125mm_failure = np.array([[0.14, 0.14, 0.14, 0.11, 0.11, 0.1, 0.11, 0.11, 0.11, 0.12, 0.12, 0.11, 0.11], [128, 128.1, 128.2, 127.9, 127.8, 127.7, 127.6, 127.5, 127.4, 127.3, 127.2, 127.75, 127.65]])
exp_125mm_success= np.array([[0.16, 0.15, 0.15, 0.15, 0.13, 0.12, 0.13, 0.12, 0.12, 0.15, 0.13, 0.12, 0.11, 0.12, 0.12, 0.12, 0.13, 0.13], [128, 128, 	128.1, 	127.9, 	127.9, 	127.9, 	127.8, 	127.8, 	127.8, 	127.7, 	127.7, 	127.7, 	127.7, 	127.6, 	127.5, 	127.4, 	127.3, 	127.2]])

exp_duallens_success = np.array([[0.12 ,0.12 ,0.11 ,0.11 ,0.11 ,0.12 ,0.12 ,0.11 ,0.11 ,0.11 ,0.11 ,0.11 ,0.11, 0.12],[117.1, 117.2, 117.2, 117.2, 117.2, 117.3, 116.9, 117.2, 117.2, 117.2, 117.2, 117.2, 117.2, 117.0]])
exp_duallens_failure = np.array([[0.11, 0.10, 0.11, 0.11, 0.11, 0.11],[117.1, 117.2, 117.3, 116.9, 117.2, 117.0]])


plt.figure()
plt.scatter(exp_125mm_failure[1,:], exp_125mm_failure[0,:]*3080, marker='x', color='red', label='Failures')
plt.scatter(exp_125mm_success[1,:], exp_125mm_success[0,:]*3080, marker='o', color='blue', label='Success')

plt.xlabel('Distance $d$ [mm]')
plt.ylabel('Laser power $P$ [W]')
plt.ylim(300,500)
plt.xlim(127.1,128.3)
plt.legend()
plt.grid(which='both')
plt.savefig('Thesis/assets/5 results/125mm_focus_threshold.pdf')
plt.show()


plt.figure()
plt.scatter(exp_duallens_failure[1,:], exp_duallens_failure[0,:]*3080, marker='x', color='red', label='Failures')
plt.scatter(exp_duallens_success[1,:], exp_duallens_success[0,:]*3080, marker='o', color='blue', label='Success')

plt.xlabel('Distance $d$ [mm]')
plt.ylabel('Laser power $P$ [W]')
plt.ylim(300,400)
plt.xlim(116.8,117.4)
plt.legend()
plt.grid(which='both')
plt.savefig('Thesis/assets/5 results/duallens_focus_threshold.pdf')
plt.show()