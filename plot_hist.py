import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
avec_data = np.loadtxt("avec_output.txt")
print(avec_data)

positions = np.arange(len(avec_data))
binary_positions = [format(i, '0' + str(int(np.log2(len(avec_data)))) + 'b') for i in range(len(avec_data))]

# Create a histogram
plt.bar(binary_positions, avec_data, color='blue')
plt.xlabel('Final states', fontsize = 20)
plt.xticks(rotation=45, ha='right', fontsize = 10)
plt.yticks(fontsize = 20)
plt.ylabel('Probabilities', fontsize = 20)
plt.title('Probability distribution',fontsize = 40)
plt.grid(True)
plt.show()
