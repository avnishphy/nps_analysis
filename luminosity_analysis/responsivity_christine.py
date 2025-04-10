import matplotlib.pyplot as plt
import numpy as np

unser_current = np.arange(2, 40, 1)

bcm2_gain = 5707
bcm2_offset = 249300

bcm4a_gain = 9597
bcm4a_offset = -1839

bcm2_current = bcm2_gain * unser_current + bcm2_offset
bcm4a_current = bcm4a_gain * unser_current + bcm4a_offset

bcm2_response = bcm2_current/unser_current
bcm4a_response = bcm4a_current/unser_current

# plt.plot(unser_current, bcm2_response/1000, label='BCM2')
plt.plot(unser_current, bcm4a_response, label='BCM4A')
plt.ylim(7300,11400)
plt.xlabel('Unser Current')
plt.ylabel('Response')
plt.title('Response vs Unser Current')
plt.legend()
plt.grid(True)
plt.show()
