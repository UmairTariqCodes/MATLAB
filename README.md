# MATLAB Portfolio Problems
These Portfolio problems demonstrate comprehensive knowledge and understanding of mathematical and computational models relevant to the engineering discipline, and an appreciation of their limitations. They underline the ability monitor, interpret and apply the results of analysis and modelling in order to bring about continuous improvement.

Candidate Number: 249383

## Portfolio Problem: [01 - Network_Problem (2).pptx](https://github.com/user-attachments/files/16193376/01.-.Network_Problem.2.pptx) 
Water Supply Network Design and Analysis

Objective:
To design and analyse a water supply network servicing a small industrial estate, focusing on determining the unknown volumetric flow rates in each pipe and the pressure loss from the source to each unit.

Network Specifications: 
6 nodes, each linked to an industrial unit. 7 internal pipes and 2 supply pipes from an off-site network with known source pressure. 6 industrial units (Unit 1 to Unit 6).

Tasks:
Determine Flow Rates: Calculate the flow rates (Q) for all pipes in the network.
Calculate Pressure Loss: Determine the pressure loss from the source to each unit.

## Portfolio Problem: [02 - Supply_Variation.pptx](https://github.com/user-attachments/files/16193378/02.-.Supply_Variation.pptx) 
Water Demand Variation Analysis

Objective:
Investigate the supply pressure variations due to changes in water demand across six industrial units, using flow rate data recorded every minute over a 24-hour period.

Task 1: Smooth Curve Fitting and Pressure Analysis
Fit a smooth curve to the noisy flow rate data for each unit.
For each minute with demand, re-run the calculations from Problem 1 to determine pipe flow rates and supply pressure loss.
Present a graph of daily variation in water pressure for each unit.

Task 2: Demand Exceedance Analysis
Identify times when each unit's water usage exceeds its declared demand using numerical root finding techniques.
Determine intersection points of the smooth curve with the declared usage value.
Graph unit supply pressure versus time of day.

## Portfolio Problem: [03 - Cooling_Fin (1).pptx](https://github.com/user-attachments/files/16193379/03.-.Cooling_Fin.1.pptx)
Cooling Fin Thermal Simulation

Objective:
To set up and solve a 2D numerical simulation for the internal transient temperature distribution of a cooling fin used in an industrial process within one of the buildings.

Network Specifications:
Geometry:
Length (L): 4 cm
Thickness (th): 1 cm
Material Properties:
Thermal Conductivity (k): 1 W/m-K
Density (ρ): 2100 kg/m³
Specific Heat Capacity (c): 850 J/kg-K
Parameters Based on Candidate Number:

Heat Transfer Coefficient (h): Derived from the 3rd, 4th, and 5th digits of the candidate number.
Heat Flux (q''): Calculated by dividing the first 5 digits of the candidate number by the last digit.

Tasks:
Define equations for nine types of nodes (internal, boundaries, and corners) using appropriate boundary conditions.
Apply finite differencing (1st order forward for time, 2nd order central for spatial).

Output Results:
Generate plots of temperature distribution at five equally spaced points along the mid-plane and top boundary.
Create contour plots of the 2D temperature distribution.
Include detailed solution steps and MATLAB code in the portfolio write-up.

## Portfolio Problem: [04 - Mixing_Tanks (1).pptx](https://github.com/user-attachments/files/16193382/04.-.Mixing_Tanks.1.pptx)
Mixing Tanks – Initial Value Problems

Objective:
To derive and solve a system of first-order differential equations for a network of mixing tanks involving time-dependent injections of three chemical components.

Network Specifications:

Tanks: 4 tanks with equal volume.
Chemical Concentrations: Each tank has three chemicals (α, β, γ) with time-dependent injections.
Flow Rates:
Various flow rates between tanks and inlets as depicted in the network diagram.
Parameters Based on Candidate Number:

Injection Rates:
𝐴1=3rd digit of candidate number×100
𝐵1=4th digit of candidate number×20
C1=5th digit of candidate number×100
B2=Tuning parameter (initially set to 40, adjust for balance).

Tasks:
Derive System of Equations:
Set notation for a system of first-order differential equations for each chemical in each tank.
Derive 12 equations (3 per tank) based on the network layout and time-varying concentration inputs.
Represent equations using α, β, γ, and then in a simplified notation for coding.

Output Results:
Produce plots showing the concentration of chemicals over time in each tank.
Validate results and ensure steady state is reached before finalizing the timespan.

