# Acoustic Plate Reader

# Introduction
![APR Workflow](https://github.com/shapiro-lab/acoustic-plate-reader/assets/14302923/3c5e6829-79da-44d7-9f6f-d468f0007d1f)


# Construction
## Motorized stage
The motorized stage used to scan the transducer over the samples consists of three components: the stage itself, the acrylic base that supports it, and the Matlab scripts used to control it. The scripts can be found [here](https://github.com/drmittelstein/ultrasound_hardware_control?tab=readme-ov-file#scripts-in-this-repository), but are provided only for troubleshooting purposes; the control of the stage during a scan is built into the Data Acquisition GUI. The parts required to build the stage and base are provided in the tables below; the assembly of the stage is self-explanatory, and the CAD files for machining the base can be found [here](https://github.com/drmittelstein/ultrasound_rigs/tree/master/Mech%20Translating%20Stage%20Setup).

**Motorized stage: BiSlide Motorized Translating System: 24”, 12”, 8”, (5 um resolution) from Velmex**

<img width="112" alt="Screenshot 2024-04-23 at 2 57 05 PM" src="https://github.com/shapiro-lab/acoustic-plate-reader/assets/14302923/5f860fbe-5f89-4a33-a814-68213c507a29">

| Velmex catalog number | Part name                                                   | Part role      |
|-----------------------|-------------------------------------------------------------|----------------|
| MN10-0240-M02-21      | BiSlide, T=24", 2mm/rev, limits, NEMA 23 mount             | X-axis BiSlide |
| PK266-03A-P1         | Vexta Type 23T2, Single Shaft Stepper Motor                | X-axis BiSlide |
| MC-3                  | Cleat, 2 hole BiSlide for English or Metric optical tables | X-axis BiSlide |
| MB-1                  | BiSlide Bolt 1/4-20x3/4" Socket cap screw (10 pack)        | X-axis BiSlide |
| MN10-0120-M02-21      | BiSlide, travel=12 inch, 2 mm/rev, limits, NEMA 23        | Y-axis BiSlide |
| PK266-03A-P1         | Vexta Type 23T2, Single Shaft Stepper Motor                | Y-axis BiSlide |
| MC-2                  | Cleat, 2 hole BiSlide                                      | Y-axis BiSlide |
| MN10-0080-M02-21      | BiSlide, travel=8 inch, 2 mm/rev, limits, NEMA 23         | Z-axis BiSlide |
| PK266-03B-P2         | Vexta Type 23T2, Double Shaft Stepper Motor               | Z-axis BiSlide |
| D6CL-6.3F             | D6CL Damper for Type 23 Double Shaft Stepper Motor         | Z-axis BiSlide |
| VXM-3              | Control, 3-axis programmable stepping motor control         | Motor controller |
| VXM-USB-RS232         | USB to RS232 Serial Communication Cable 10ft              | Cable |

**Acrylic base**

| McMaster-Carr catalog number | Material name                                                    | Material role        |
|------------------------------|------------------------------------------------------------------|----------------------|
| 8560K325                     | Clear Scratch- and UV-Resistant Cast Acrylic Sheet (24" x 36" x 1") | Base                 |
| 8560K914                     | Clear Scratch- and UV-Resistant Cast Acrylic Sheet (12" x 12" x 2") | Material to make pillars |
| 1227T459                    | Clear Scratch- and UV-Resistant Cast Acrylic Bar (3" Wide x 1/2" Thick) | Material to make adapters |

Dimensions of mounted stage: 18” x 28”

<img width="392" alt="Screenshot 2024-04-23 at 2 59 03 PM" src="https://github.com/shapiro-lab/acoustic-plate-reader/assets/14302923/ac41ef90-97af-40c7-ad17-acd20938071e">


## Imaging water tank
The water tank used in this study is designed to hold 12 of the 96-well phantoms. The parts were made by laser cutting, assembled with xxx(glue) and water-proofed with xxx(glue2?). For the walls (Long_v4.DXF and Short_v4.DXF) and the base insert (Base_horiz.DXF), 0.25 inch acrylic sheets were used, while for the bottom piece (Inner.DXF), 0.5 in ones were used for improved rigidity. Alternatively, two pieces of the bottom cut from 0.25 in acrylic sheets can be stacked and glued together. The simulated full assembly of the tank is shown in Full_Assemble.SLDASM for the reference.
## 96-well phantoms
The phantom molds were 3D-printed (96-well_PhantomMold_v1.STL and 96-well_PhantomMold_body_shallow.STL). Specifically, the molding parts (96-well_PhantomMold_v1.STL) were printed with Nylon 12 using Selective Laser Sintering through Shapeways for the finer features of the wells. The walls (96-well_PhantomMold_body_shallow.STL) were printed with PLA in-house. 
## Transducer holders
The holder for the transducer used in this study (Verasonics L22-14vX) was 3D-printed with PLA in-house (L22_holder_15deg.SLDPRT and L22_holder_top.SLDPRT).

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
5. Confirm that everything is correct, and click "Start Scan"

<img width="774" alt="Screenshot 2023-04-12 at 12 07 20 PM" src="https://user-images.githubusercontent.com/14302923/231559930-73bba2a2-abcb-42e9-86f7-af7056e918a1.png">

6. Looking at the gray Bmode iamge, posisiton the first well as it appears in the screenshot. The depth will depend on the transducer used: the sample should be at 5 mm for the L22, and 20 mm for the L10.
7. In the gray VSX GUI, click "Start Ramp" (you may have to click it many times for the scan to start).

# Processing data
The Data Processing GUI was designed to process data for any of the three scan types that was acquired using the Data Acquisition GUI. 

## Processing data via the GUI

# Common issues
## Unable to connect to motor stage
Open Device Manager, and under Ports, look for "Prolific USB-to-Serial Comm Port (COM2)". If the device doesn't show up, restart and replug the motor stage. If it shows up in the list but isn't on COM2, go to Properties->Port Settings->Advanced and switch it. If this doesn't work, try re-installing the driver for the motor stage from the Velmex website.
## Scan stops before finishing
This is usually caused by a problem connecting with the motor stage, or rarely by a Windows Update that forces the computer to restart.
