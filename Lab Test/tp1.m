M = 1;
K = 4e2;
gamma = 0.01;


Fs = 100;            % Sampling frequency            
fsin = 5;
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
tspan = (0:L-1)*T;        % Time vector
b=1


c = gamma*2*((K*M)^0.5);
A = [0 1
    -K/M -c/M]
v = [0 1*b]';
B = [0;1]

% Initial Conditions
u= [0.00;0];

% Derivative Function
fmin=15;fmax=25;
PA.u=@(t) cos((fmin+(fmax-fmin)/2*t/tspan(end))*t); % Q2.2 2. No Leakage
%PA.u= @(t) sin(198*pi*2*t); % Q2.2 2. Leakage
der = @(t,x) A*x + v*(PA.u(t)); % Q2
%der = @(t,x) A*x + v'*u; % Q1


% Response

[ti,yl] = ode45(der, tspan, u)
%ind=(ti>6); ti=ti(ind); yl=yl(ind,:); % Truncation

%plot(ti, yl(:,1));
%hold on 
%% vectorized
close all; clear all; clc;
M = 1;
K = 4e2;


Fs = 100;            % Sampling frequency            
fsin = 5;
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
tspan = (0:L-1)*T;        % Time vector
b=1;
gamma = [0.01 0.03 0.05 0.1];

yl = zeros(length(tspan), length(gamma)*2);
ti = zeros(length(tspan), length(gamma));
ft = zeros(length(tspan), length(gamma));
vari = [];
u= [0.00;0];
v = [0;1];
B = [0;1];
fmin=2*pi*2;fmax=2*pi*10;
PA.u=@(t) cos((fmin+(fmax-fmin)/2*t/tspan(end))*t); % Q2.2 2. No Leakage

for i=1:length(gamma)
    c = gamma(i)*2*((K*M)^0.5);
    A = [0 1
        -K/M -c/M];
 
 

% Initial Conditions


% Derivative Function
%PA.u= @(t) sin(198*pi*2*t); % Q2.2 2. Leakage
    der = @(t,x) A*x + v*(PA.u(t)); % Q2
%der = @(t,x) A*x + v'*u; % Q1


% Response

    [ti(:,i),yl(:, 2*i-1:2*i)] = ode45(der, tspan, u);
%ind=(ti>6); ti=ti(ind); yl=yl(ind,:); % Truncation

    plot(ti(:,i), yl(:,2*i-1));
%hold on 
% FFT

    ft(:,i) = fft(yl(:,2*i-1));


% % Plot Fourier Transform
% P2 = abs(ft/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% plot(f,P1)
% hold on;


    N=size(ft,1);
    freq=0 : 1 /(N*T) : (N-1)/N/T;
    
    %plot(freq,(abs(ft(:,i))));
    hold on;
    vari(end+1) = var(yl(:,2*i-1));
end
legend('\zeta = 0.01', '\zeta = 0.03', '\zeta = 0.05', '\zeta = 0.1');
%xlabel('Frequency')
xlabel('Time')
ylabel('Amplitude');
title('Effects of Damping on Variance')
zeta = gamma
vari
%%
M = 1;
K = 4e2;
gamma = 0.05;

Fs = 100;            % Sampling frequency            
fsin = 5;
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
tspan = (0:L-1)*T;



c = gamma*2*((K*M)^0.5);
A = [0 1
    -K/M -c/M]

B = [0;1]

Yo=[0 ; 0];

PA=struct('A', A ,'B',[0;1],'u',[]);
fmin=15;fmax=25;PA.u=@(t) cos((fmin+(fmax-fmin)/2*t/tspan(end))*t);
der = @(t,x) PA.A*x + B*(PA.u(t)); % Q2
%der = @(t,x) A*x + v'*u; % Q1


% Response

[ti,yl] = ode45(der, tspan, PA.u)
