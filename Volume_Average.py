from odbAccess import *
import numpy as np
import os

Target_Strain = 0.035
Epsilon = np.log(1+Target_Strain)

# Get any .odb file in the current directory
odb_files = [f for f in os.listdir('.') if f.endswith('.odb')]
if len(odb_files) == 0:
    raise FileNotFoundError("No .odb files found in the current directory.")
odb_name = odb_files[0]  # Use the first .odb file found
output_file_name = "SS_Curve.txt"

# Open the output database
odb = openOdb(path=odb_name)

# Access the first step
step_name = list(odb.steps.keys())[0]
step = odb.steps[step_name]

# Prepare for output
output = []

# Loop through each frame
for frame_idx, frame in enumerate(step.frames):
    print(f"Processing frame {frame_idx + 1}/{len(step.frames)}")
    
    total_stress_volume = 0.0
    total_volume = 0.0

    # Get the stress field and element volume field
    stress_field = frame.fieldOutputs['S']
    volume_field = frame.fieldOutputs['EVOL']

    # Loop through all elements with stress and volume data
    stress_values = stress_field.values
    volume_values = volume_field.values

    # Create a dictionary for quick volume lookup by element label
    volume_dict = {v.elementLabel: v.data for v in volume_values}

    for stress in stress_values:
        element_label = stress.elementLabel
        if element_label in volume_dict:
            element_volume = volume_dict[element_label]
            s_mises = stress.mises  # Von Mises stress
            total_stress_volume += s_mises * element_volume
            total_volume += element_volume

    # Compute volume-averaged stress for this frame
    volume_avg_stress = total_stress_volume / total_volume if total_volume > 0 else 0.0
    # Multiply time by 0.035
    strain = frame.frameValue * Epsilon
    output.append((strain, volume_avg_stress))

# Write results to file
with open(output_file_name, "w") as file:
    file.write("strain\tVolume_Averaged_Stress\n")
    for strain, avg_stress in output:
        file.write(f"{strain:.4f}\t{avg_stress:.4f}\n")

print(f"Results written to {output_file_name}")

# Close the ODB
odb.close()
