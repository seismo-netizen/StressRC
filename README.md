Readme

This package of codes read in the simulated dynamic wavefield of two major earthquakes of the 2019 Ridgecrest Earthquake sequence, compute and plot stress on Parkfield section of the San Andrews fault.

Simulation code: SPECFEM3D
3D velocity model: CVMS4.26
Finite fault model: Yue et al. 2020, in preparation

Parkfield fault SE: lat = 35.0859, lon =  -119.6080;
Parkfield fault NW: lat = 36.4681, lon = -121.0391;

Reference frame:
x: east; y: north; z: up
Tensile positive, compression negative

Contents:
Matlab code: S1_stress.m, read in data, compute stress, plot figures
Matlab code: S2_img2vdo.m, merge figures to avi video

Matlab data: pk_str_RC71.mat, dynamic strain of the Mw 7.1 mainshock
Matlab data: pk_str_RC64.mat, dynamic strain of the Mw 6.4 foreshock

Covert video from avi to mov via online tool: https://www.media.io/convert/avi-to-mov.html

Video: RC71_PK_STRESS.mov, animation of dynamic stress of the Mw 7.1 mainshock on fault
Video: RC64_PK_STRESS.mov, animation of dynamic stress of the Mw 6.4 foreshock on fault

Figure: RC64_PK_STRESS_095.png, dynamic stress of the Mw 7.1 mainshock on fault at t = 86.8 sec
Figure: RC71_PK_STRESS_095.png, dynamic stress of the Mw 6.4 foreshock on fault at t = 86.8 sec
Figure: Wvfld_fault.png, Z component displacement at t = 78.8 sec, bold red line shows the Parkfield section 

Folder: figures, example figures generated by S1_stress.m
