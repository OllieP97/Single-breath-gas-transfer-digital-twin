This repository contains MATLAB code used to implement and reproduce the mechanistic digital-twin simulations described in the
manuscript 'Evaluation of a digital twin of blood-gas barrier function in respiratory diseases'.

The model simulates a 10-s single-breath carbon monoxide gas-transfer manoeuvre using a four-compartment representation of the lungs and returns simulated DLCO and KCO.

**'Simulate_KCO.m'** 
is the main function for core four-compartment CO gas-transfer simulation.

**'Figure4_Sweep.m'** 
reproduces the mechanistic parameter sweeps used for Figure 4, calling Simulate_KCO and generating the corresponding DLCO
and KCO response surfaces.

**'MRI.m'**
runs participant-specific simulations using hyperpolarised Xe-129 MRI-derived inputs. Subject-specific healthy reference simulations are first calibrated to GLI-predicted DLCO. Patient simulations then incorporate available physiological and MRI-derived measurements, including RV and TLC, functional residual capacity, ventilation defect percentage, regional lung-volume distribution and Xe-129 RBC:barrier ratio. The RBC:barrier ratio is used as a relative scaling term for gas-transfer efficiency.

**'Simulate_KCO_qCT_withSA.m'**
is the main function with the addition of a dimensionless surface-area scaling term ('SA_Scale') representing the relative functional gas-exchange surface area retained within the model. This can be used to run ILD patient simulations with qCT derived inputs.

**'QCT.m'**
runs participant-specific qCT-informed simulations for the ILD cohort incorporating subject specific inputs including GLI healthy reference values, RV and TLC alveolar-volume fraction, regional lung-volume distribution, regional pulmonary vascular distribution, alveolar-capillary membrane thickness scaling and functional gas-exchange surface-area scaling.

**Running the Simulations:**

No specialist MATLAB toolboxes are required for the core simulations.
Some plotting or analysis functions may require standard MATLAB statistical functionality depending on the installed MATLAB version. All .m files can be placed in the same MATLAB working directory.
Do not run 'Simulate_KCO.m' or 'Simulate_KCO_qCT_withSA.m' directly. These files define functions that are called by the driver scripts.
