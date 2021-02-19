There are 3 Simulink Model, 4 MATLAB code, and 1 MATLAB data file is attached. 

Keep all the files in same folder.

1. Project.m - to run the roll control channel response and its robustness analysis.

2. To run the autopilot
	- run the Project.m file first.
	- run alt_turn.slx for the simulation without noise or disturbance.
	- run alt_turn_noise.slx for the simulation with noise or disturbance.

3. To run the S-function
	- Keep CRJ330Config.m, CRJ330.m, Data.mat, get_x_dot2.m and alt_turn_sfun.slx file in same folder.
	- run Data.mat file first in MATLAB command window.
	- run the Simulink file alt_turn_sfun.slx
