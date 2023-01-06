clear all; close all; clc;

%% Load Model
load cc_model
cf=feplot; cf.model=model; cf.def=def;fecom('ColorDataEvalY');
model
def
%% Q:2 Damping
load cc_model
cf=feplot; cf.model=model; cf.def=def;fecom('ColorDataEvalY');
damp=.01;
f=linspace(0,5000,2048)'; % frequencies
sys=nor2ss(def,damp, 426.02 , 1976.02 ,'Hz acc');
%qbode(sys,f*2*pi,'iiplot "simul1" -po');
ci=iiplot;iicom('sub1 1'); % pointer to responses
ci.Stack{'simul1'} % Data structure containing response

% Different damping ratios change the values of the peak and also sometimes
% displaces them, the bandwidth is also changed


sys1=nor2ss(def, 0.1 , 426.02 , 1976.02 ,'Hz acc');
qbode(sys1,f*2*pi,'iiplot "simul 0.1" -po');
sys2=nor2ss(def, 0.5 , 426.02 , 1976.02 ,'Hz acc');
qbode(sys2,f*2*pi,'iiplot "simul 0.5" -po');

%% Raleigh

damp_rat = 0.01;
freq = def.data(7);
beta = 2*damp_rat/freq

ral_freq = beta*def.data/2;

sys3=nor2ss(def, 0.01 , 911.02, 1976.02 ,'Hz acc');
qbode(sys3,f*2*pi,'iiplot "simul Ral 911" -po');


sys4=nor2ss(def, 0.01 , 426.02, 1976.02 ,'Hz acc');
qbode(sys4,f*2*pi,'iiplot "simul Ral 426.02" -po');

% The first part of the equation 2 is the [c][phi], which is the 
% observability and [phi]'[b] is the controllability, the force thingy

%% Comandibility Computation
% commandability computation 
[(1:size(def.def,2))' def.def'*fe_c(def.DOF,[426;911]+.02)']
% Show node numbers
fecom('ColorDataEvalY -alpha.3');%fecom('showline');
fecom('textnode 426 911','FontSize',16,'Color','r')

%% Time Response
load cc_model
cf=feplot; cf.model=model; cf.def=def;
uf.t=linspace(0,.5,2048*3)';  % time vector
uf.u0=uf.t*0;uf.u0(find(uf.t<=1e-4))=1;  % ch1 (input =impact 0.1 ms);
uf.window='None';  uf.noise=0; uf.filt=[];
damp=.001; 
uf.sys=nor2ss(def,damp,426.02,1976.02,'Hz acc'); % State space model
uf=cc_simul('simul',uf);cc_simul('plot',uf);

%% Time Response Buttered (Anti-Aliased)
load cc_model; cf=feplot; cf.model=model; cf.def=def;
uf.t=linspace(0,.5,2048)';  % time vector
fmax=1500;[num,den]=butter(16,fmax*2*pi,'s');uf.filt=tf(num,den);
uf.u0=uf.t*0;uf.u0(find(uf.t<=1e-4))=1;  % ch1 (input impact 0.1 ms);
uf.window='None';  uf.noise=0.01;damp=.001; 
uf.sys=nor2ss(def,damp,426.02,1976.02,'Hz acc'); % state-space model
uf=cc_simul('simul',uf);cc_simul('plot',uf)