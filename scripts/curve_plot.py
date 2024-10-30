import matplotlib.pyplot as plt
import numpy as np

# Function to read data from the output file
def read_data(filename):
    sample_sizes = []
    avg_gene_counts = []
    std_devs = []

    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                sample_sizes.append(int(parts[0]))        # Sample Size
                avg_gene_counts.append(float(parts[1]))   # Average Gene Count
                std_devs.append(float(parts[2]))          # Standard Deviation

    return np.array(sample_sizes), np.array(avg_gene_counts), np.array(std_devs)

# Function to plot the data with error bars (standard deviation)
def plot_data(sample_sizes, avg_gene_counts, std_devs):
    plt.figure(figsize=(10, 6))

    # Plot the curve with error bars
    plt.errorbar(sample_sizes, avg_gene_counts, yerr=std_devs, fmt='o-', capsize=5, label='All', ecolor='gray', elinewidth=2, markeredgewidth=2, color='b')

    # Add labels and title
    plt.xlabel('Sample size', fontsize=22)
    plt.ylabel('Gene count', fontsize=22)
#    plt.title('Gene Count vs. Sample Size with Standard Deviation')
    plt.grid(True)

    # Increase font size for tick labels
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)   

    # Show legend
    plt.legend(fontsize=22)

    # Show the plot
    plt.tight_layout()
    plt.show()

# Main function to run the script
def main():
    # Path to the output file from the Perl script
    output_file = 'output.txt'  # Change this to the actual path

    # Read data from file
    sample_sizes, avg_gene_counts, std_devs = read_data(output_file)

    # Plot the data
    plot_data(sample_sizes, avg_gene_counts, std_devs)

if __name__ == '__main__':
    main()
