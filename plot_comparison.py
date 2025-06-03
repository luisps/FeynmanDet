import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
def plot_comparison(file1_path, file2_path):
    data1 = np.loadtxt(file1_path)
    data2 = np.loadtxt(file2_path)

    values1 = data1
    values2 = data2
    print(len(values1))
    print(len(values2))
    bar_width = 0.35
    index = np.arange(len(values1))

    plt.bar(index, (values1), bar_width, label='Red/Green', color='orange')
    plt.bar([i + bar_width for i in index], (values2), bar_width, label='Brute Force', color='blue')

    CircuitNames = ["3x3 trivial circuit", "4x4 highly branching circuit", "3x3 Haddamard circuit", "6 qubit Random circuit with depth 5", "7 qubit random circuit with depth 5"]
    plt.xlabel('Circuits', fontsize=20)
    plt.ylabel('$\log_{10}$(Timing)', fontsize=20)
    plt.title('Comparison of Red/Green and Brute Force Timings in seconds', fontsize=20)
    plt.xticks([i + bar_width/2 for i in index], CircuitNames, rotation=10)
    plt.legend()

    # Add labels for each pair of bars
    for i, (v1, v2) in enumerate(zip(values1, values2)):
        if v1 >= 0:
            plt.text(i, v1, f'{v1:.2f}', ha='center', va='bottom', color='black', fontsize=8)
        else:
            plt.text(i, v1, f'{v1:.2f}', ha='center', va='top', color='black', fontsize=8)

        if v2 >= 0:
            plt.text(i + bar_width, v2, f'{v2:.2f}', ha='center', va='bottom', color='black', fontsize=8)
        else:
            plt.text(i + bar_width, v2, f'{v2:.2f}', ha='center', va='top', color='black', fontsize=8)

    plt.yscale('log')  # Set the y-axis to log scale
    plt.tight_layout()
    plt.show()

# Substitua 'timing_output.txt' e 'timing_output_unoptimized.txt' pelos caminhos reais dos seus arquivos
plot_comparison('timing_output.txt', 'timing_output_unoptimized.txt')
