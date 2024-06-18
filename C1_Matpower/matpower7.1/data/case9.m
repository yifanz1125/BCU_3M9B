function mpc = case9
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from p. 70 of:
%
%   Chow, J. H., editor. Time-Scale Modeling of Dynamic Networks with
%   Applications to Power Systems. Springer-Verlag, 1982.
%   Part of the Lecture Notes in Control and Information Sciences book
%   series (LNCIS, volume 46)
%
%   which in turn appears to come from:
%
%   R.P. Schulz, A.E. Turner and D.N. Ewart, "Long Term Power System
%   Dynamics," EPRI Report 90-7-0, Palo Alto, California, 1974.

%   MATPOWER
%% MATPOWER Case Format : Version 2
mpc.version = '2';


%% Interaction settings
[Datatitle,Datapath]=Fun_getGlobal;
% Datapath='E:\OneDrive - zju.edu.cn\A_FilesCloud\A1_System\A_Simulations\B2_3M9B\Data\';
addpath(genpath(Datapath));
% Datatitle='Temp';%'BCUTb121_EP_Prefault';
prefile=['Preset_' Datatitle,'.mat'];
load(prefile);
if (preset.flagdata==0)

    %%-----  Power Flow Data  -----%%
    %% system MVA base
    baseMVA = 100;
    
    %% bus data
    %	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
    bus = [
	    1	3	0	0	0	0	1	1	0	16.5	1	1.1	0.9;
	    2	2	0	0	0	0	1	1	0	18	1	1.1	0.9;
	    3	2	0	0	0	0	1	1	0	13.8	1	1.1	0.9;
	    4	1	0	0	0	0	1	1	0	230	1	1.1	0.9;
	    5	1	125	50	0	0	1	1	0	230	1	1.1	0.9;
	    6	1	90	30	0	0	1	1	0	230	1	1.1	0.9;
	    7	1	0	0	0	0	1	1	0	230	1	1.1	0.9;
	    8	1	100	35	0	0	1	1	0	230	1	1.1	0.9;
	    9	1	0	0	0	0	1	1	0	230	1	1.1	0.9;
    ];
    
    %% generator data
    %	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    % mpc.gen = [
    % 	1	72.3	0	300	-300	1.04	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
    % 	2	163	0	300	-300	1.025	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
    % 	3	85	0	300	-300	1.025	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
    % ];
    % mpc.gen = [
    % 	1	89.8	0	300	-300	1.1	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
    % 	2	134.32	0	300	-300	1.0974	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
    % 	3	94.19	0	300	-300	1.0866	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
    % ];
    gen = [
	    1	89.8	0	300	-300	1.1083	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
	    2	134.32	0	300	-300	1.1071	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	    3	94.19	0	300	-300	1.0606	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
    ];
    %% branch data
    %	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
    % mpc.branch = [
    % 	1	4	0	0.0576	0	250	250	250	0	0	1	-360	360;
    % 	4	6	0.017	0.092	0.158	250	250	250	0	0	1	-360	360;
    % 	6	9	0.039	0.17	0.358	150	150	150	0	0	1	-360	360;
    % 	3	9	0	0.0586	0	300	300	300	0	0	1	-360	360;
    % 	8	9	0.0119	0.1008	0.209	150	150	150	0	0	1	-360	360;
    % 	8	7	0.0085	0.072	0.149	250	250	250	0	0	1	-360	360;
    % 	7	2	0	0.0625	0	250	250	250	0	0	1	-360	360;
    % 	7	5	0.032	0.161	0.306	250	250	250	0	0	1	-360	360;
    % 	5	4	0.01	0.085	0.176	250	250	250	0	0	1	-360	360;
    % ];
    branch = [
	    1	4	0	0.0576+0.0608	0	250	250	250	0	0	1	-360	360;
	    4	6	0.017	0.092	0.158	250	250	250	0	0	1	-360	360;
	    6	9	0.039	0.17	0.358	150	150	150	0	0	1	-360	360;
	    3	9	0	0.0586+0.1813	0	300	300	300	0	0	1	-360	360;
	    8	9	0.0119	0.1008	0.209	150	150	150	0	0	1	-360	360;
	    8	7	0.0085	0.072	0.149	250	250	250	0	0	1	-360	360;
	    7	2	0	0.0625+0.1198	0	250	250	250	0	0	1	-360	360;
	    7	5	0.032	0.161	0.306	250	250	250	0	0	1	-360	360;
	    5	4	0.01	0.085	0.176	250	250	250	0	0	1	-360	360;
    ];
    
else
    bus=preset.bus;
    gen=preset.gen;
    branch=preset.branch;
    baseMVA=preset.Sbase;
end


%%-----  OPF Data  -----%%
mpc.bus=bus;
mpc.gen=gen;
mpc.branch=branch;
mpc.baseMVA=baseMVA;

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
];
