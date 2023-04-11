# Acoustic Plate Reader

# Introduction

# Construction

# Running a scan
## Types of scans
The most common types of scans we run are:
* Pre-/post-collapse imaging: this is the fastest scan type, and is good for checking if samples are producing GVs or not. It requires that you know the appropriate voltage to image your GV type of interest without collapsing them, which can be determined from an imaging voltage ramp.
* Imaging voltage ramps: this scan type is useful for characterizing the way a GV type responds to varying acoustic pressure, either to find a pressure that will image the GVs without collapsing them or for other applications like multiplexed imaging.
* Collapse voltage ramps: this scan type will acquire acoustic collapse data for your samples. Imaging is performed using whatever "beam type" is specified (xAM recommended), and collapse is performed using a parabolic beam.

## Starting a scan via the GUI
1. Enter the name you wish the data folder to have in the "save name" field
2. Click "Select Vantage directory" and navigate to the Vantage directory you wish to run the data acquisition script from (this will likely have a name of the form "Vantage-4.6.2-XXX" and should contain a script called "activate.p"
3. Select the transducer, beam type, plate parameters, and number of equilibration frames you want. You will rarely need to change the starting well/row/column fields, and you should make sure that the numbers of rows/columns matches the number of plates you wish to scan. If you set these numbers to be too high, the motor stage will not automatically recognize this error, and may damage the transducer. The number of equilibration frames is the number of frames that occur between a change in voltage and the imaging being saved for that voltage, and is intended to make sure that if any GVs are going to collapse at that imaging voltage, they do so before the final image for that voltage is aquired.
4. Select the scan type, and set the voltages for the scan as desired.


# Processing data
The Data Processing GUI was designed to process data for any of the three scan types that was acquired using the Data Acquisition GUI. 

## Processing data via the GUI

# Common issues
## Unable to connect to motor stage
Open Device Manager, and under Ports, look for "Prolific USB-to-Serial Comm Port (COM2)". If the device doesn't show up, restart and replug the motor stage. If it shows up in the list but isn't on COM2, go to Properties->Port Settings->Advanced and switch it. If this doesn't work, try re-installing the driver for the motor stage from the Velmex website.
## Scan stops before finishing
This is usually caused by a problem connecting with the motor stage, or rarely by a Windows Update that forces the computer to restart.
