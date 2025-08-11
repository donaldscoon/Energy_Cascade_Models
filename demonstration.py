from Dynamic_ECs import CAV, BOS, AMI, EC
import os

# Get the absolute path to the directory containing the running script
base_dir = os.path.dirname(os.path.abspath(__file__))

# Input file and results directory located in the same directory as the script
input_file = os.path.join(base_dir, 'input_file.csv')
results_dir = os.path.join(base_dir, 'EC_results')

# Runs a demonstration of the Cavazzoni 2004 et. al Modified Energy Cascade
CAV(input_file=f'{dir}/input_file.csv', TCB=10, t_A=28, 
    output=results_dir, save_output=True)

# Runs a demonstration of the Boscheri 2012 et. al Modified Energy Cascade
BOS(input_file=f'{dir}/input_file.csv', TCB=10, t_A=28, 
    output=results_dir, save_output=True)

# Runs a demonstration of the Amitrano 2020 et. al Modified Energy Cascade
AMI(input_file=f'{dir}/input_file.csv', TCB=10, t_A=28, 
    output=results_dir, save_output=True)

# Runs a demonstration of the Volk 1995 et. al Energy Cascade
EC(input_file=f'{dir}/input_file.csv', TCB=10, t_A=28, 
   output=results_dir, save_output=True)