# NonLinDRHAMatlab
Matlab code for non-linear dynamic response history analysis of one-dimensional multi-degree of freedom systems

Matlab code for the application of the non-linear dynamic response history analysis (Non Lin-ear DRHA) of multi-degree of freedom (MDOF) structures is presented. Direct integration of the equations of motion is used for the calculation of the non-linear dynamic response of a MDOF system subject to dynamic loading. The dynamic response history analysis procedure with direct integration proceeds incrementally, by solving the MDOF system equations for each time step. 

For the direct integration of equations of motion, the function NLDRHA_MDOF_1D.m is used. See the following examples:
•	example_SDOF_Clough_Penzien.m
•	example_SDOF_Rajasekaran.m
•	example_Shear_Building_2_NPTEL.m
•	example_Energy_SDOF_Chopra.m
•	example_SDOF_Chopra.m
for more details. 

The function NLDRHA_MDOF_1D.m needs the function Bilinear_Kin_1D.m in order to run properly. In the last function the bilinear elastoplastic hysteretic model with elastic viscous damping is implemented for SDOF or MDOF systems. The MDOF structure modeled with this function consists of lumped masses connected with stiffness and damping elements in series, along a single dimension (i.e. one-dimensional MDOF system, hence the function is named as '1D'). Each lumped mass has one degree of freedom. The first degree of freedom is at the top of the structure and the last at its fixed base. However, the last degree of freedom
is not included in the input arguments of the function, i.e. not contained in ndof, as it is al-ways fixed. The nonlinear stiffness is virtually of the bilinear type, where an initial stiffness and a post-yield stiffness are defined. The unloading or reloading curve of this model are par-allel to the initial loading curve, and a hysteresis loop is created by continuously loading and unloading the structure above its yield limit. This behavior can be viewed as hardening of the kinematic type.

All the above functions can be used for acceleration time histories of a constant time step size. If this is not the case, then the acceleration time history needs to be resampled by using the MATLAB program file function RESAMPLE.m. The user is encouraged to see the example
•	example_Resampling_Nonuniform_Time_History.m
in this last case. 


The present code is accompanied by 6 examples in which its application is presented. These examples are taken from various standard textbooks or other material. The results of the ex-amples are verified by the results of the application of the present code.
The author is open to any suggestions or recommendations that the users may have.
