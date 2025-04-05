# -*- coding: mbcs -*-
import os
import numpy as np

# Abaqus-specific imports
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import visualization
import xyPlot
from odbAccess import openOdb

showStopButtonInGui()
def generate_ss_curve(EngStrain, Odb_Path, Output_Path, S_Comp, Plot_Flag):
    """
    A function to compute volume-averaged stress from an ODB file,
    write results to a text file, and optionally plot them in Abaqus/CAE.
    
    Parameters:
    -----------
    EngStrain : float
        Engineering strain to be converted to true strain.
    Odb_Path : str
        Full path to the ODB file.
    Output_Path : str
        Full path to the output text file.
    S_Comp : str
        Stress component to extract (e.g., "S11", "S22", etc.).
    Plot_Flag : bool
        If True, create an XY plot in Abaqus/CAE using the generated data.
    """
    # Safety check: this script must run in Abaqus' Python environment
    showStopButtonInGui()

    # Convert engineering strain to true strain
    Epsilon = np.log(1.0 + EngStrain)

    # Open the ODB
    if not os.path.isfile(Odb_Path):
        raise FileNotFoundError("The specified ODB file does not exist: " + Odb_Path)

    odb = openOdb(path=Odb_Path)

    # Get the first step in the ODB
    step_name = list(odb.steps.keys())[0]
    step = odb.steps[step_name]

    # Collect results (strain, volume-averaged stress)
    output_data = []

    # Loop over frames
    for frame_idx, frame in enumerate(step.frames):
        print("Processing frame {}/{}".format(frame_idx + 1, len(step.frames)))

        total_stress_volume = 0.0
        total_volume = 0.0

        # Get stress and element volume fields
        stress_field = frame.fieldOutputs['S']
        volume_field = frame.fieldOutputs['EVOL']

        stress_values = stress_field.values
        volume_values = volume_field.values

        # Dictionary for volume lookups by element label
        volume_dict = {v.elementLabel: v.data for v in volume_values}

        # Accumulate stress * volume
        for stress_value in stress_values:
            elem_label = stress_value.elementLabel
            if elem_label in volume_dict:
                elem_volume = volume_dict[elem_label]

                # Pick correct stress component
                s_comp_data = stress_value.data  # tuple (S11, S22, S33, S12, S13, S23)
                if S_Comp.upper() == 'S11':
                    s = s_comp_data[0]
                elif S_Comp.upper() == 'S22':
                    s = s_comp_data[1]
                elif S_Comp.upper() == 'S33':
                    s = s_comp_data[2]
                elif S_Comp.upper() == 'S12':
                    s = s_comp_data[3]
                elif S_Comp.upper() == 'S13':
                    s = s_comp_data[4]
                elif S_Comp.upper() == 'S23':
                    s = s_comp_data[5]
                else:
                    raise ValueError("Unsupported S_Comp: {}".format(S_Comp))

                total_stress_volume += s * elem_volume
                total_volume += elem_volume

        # Compute volume-averaged stress for this frame
        if total_volume > 0.0:
            volume_avg_stress = total_stress_volume / total_volume
        else:
            volume_avg_stress = 0.0

        # Compute the true strain for this frame
        # frameValue is typically the step time, so we scale it by Epsilon
        strain_val = frame.frameValue * Epsilon

        output_data.append((strain_val, volume_avg_stress))

    # Close the ODB
    odb.close()

    # Write results to the output file
    with open(Output_Path, "w") as f:
        f.write("strain\tVolume_Averaged_Stress\n")
        for strain_val, avg_stress in output_data:
            f.write(f"{strain_val:.6f}\t{avg_stress:.6f}\n")

    print("Results written to {}".format(Output_Path))

    # ----------------------------------------------------
    # If requested, read the file back and create XY Plot
    # ----------------------------------------------------
    if Plot_Flag:
        # Make/read XY data from the output file
        xy_data = []
        with open(Output_Path, "r") as f:
            # Skip header
            next(f)
            for line in f:
                parts = line.split()
                if len(parts) < 2:
                    continue
                strain_val = float(parts[0])
                stress_val = float(parts[1])
                xy_data.append((strain_val, stress_val))

        # Create (or reuse) a viewport
        if 'Viewport: 1' not in session.viewports:
            session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200, height=150)

        session.viewports['Viewport: 1'].makeCurrent()
        session.viewports['Viewport: 1'].maximize()

        # Create a unique name for the XYPlot if one already exists
        base_plot_name = 'XYPlot-VolumeAveraged'
        plot_name = base_plot_name
        i = 1
        while plot_name in session.xyPlots.keys():
            plot_name = base_plot_name + '_' + str(i)
            i += 1

        # Create the XYPlot object with the new unique name
        xyp = session.XYPlot(name=plot_name)
        chartName = xyp.charts.keys()[0]
        chart = xyp.charts[chartName]

        # Define x- and y-axis quantity types (for a typical stress-strain plot)
        xQuantity = visualization.QuantityType(type=STRAIN)
        yQuantity = visualization.QuantityType(type=STRESS)

        # Convert our Python list to an XYData object
        myXYData = xyPlot.XYData(
            data=xy_data,
            sourceDescription='Volume-averaged stress data',
            axis1QuantityType=xQuantity,
            axis2QuantityType=yQuantity
        )

        # Create a Curve
        curve = session.Curve(xyData=myXYData)
        chart.setValues(curvesToPlot=(curve,))
        session.charts[chartName].autoColor(lines=True, symbols=True)

        # Display the XYPlot
        session.viewports['Viewport: 1'].setValues(displayedObject=xyp)

        # Optionally adjust symbol style, legend, etc.
        curve.symbolStyle.setValues(show=True, size=2)
        curve.setValues(useDefault=False, legendLabel='Vol_Avg_Stress')

        print("Plot created in Abaqus/CAE window: {}".format(plot_name))
