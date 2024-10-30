import matplotlib.pyplot as plt
import numpy as np
import sys

# Function to read data from the output file
def read_data(filename):
    data = []

    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            row = [float(part) for part in parts]
            data.append(row)

    data = np.array(data)
    sample_sizes = data[:, 0]
    datasets = [data[:, i:i+2] for i in range(1, data.shape[1], 2)]
    
    return sample_sizes, datasets

# Function to plot the data with error bars (standard deviation)
def plot_data(sample_sizes, datasets):
    plt.figure(figsize=(10, 6))

    labels = ['All', 'Core', 'Softcore', 'Dispensable', 'Private']
    colors = ['b', 'g', 'r', 'c', 'm']
    
    for i, dataset in enumerate(datasets):
        avg_values = dataset[:, 0]
        std_devs = dataset[:, 1]
        plt.errorbar(sample_sizes, avg_values, yerr=std_devs, fmt='o-', capsize=5, label=labels[i], ecolor='gray', elinewidth=2, markeredgewidth=2, color=colors[i])

    # Add labels and title
    plt.xlabel('Sample size', fontsize=22)
    plt.ylabel('Gene count', fontsize=22)
    plt.grid(True)

    # Increase font size for tick labels
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)    

    # Show legend
    plt.legend(fontsize=22, loc='upper right', bbox_to_anchor=(1, 0.80))

    # Show the plot
    plt.tight_layout()
    plt.show()

# Main function to run the script
def main():
    if len(sys.argv) != 2:
        print("Usage: python curve_plot2.py <filename>")
        return

    output_file = sys.argv[1]

    # Read data from file
    sample_sizes, datasets = read_data(output_file)

    # Plot the data
    plot_data(sample_sizes, datasets)

if __name__ == '__main__':
    main()