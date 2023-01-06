%-- mr.ms2sc.16-17@listes.centralesupelec.fr
%-- elie-abiraad@ensta-paristech.fr
%-- marco.merlotti@student.ecp.fr
%-- thanh-tung.nguyen@ensta-paristech.fr

m=1;
k=400;
c=0.4;

%-----------------------------------%
%--                               --%
%-- Reponse libre - CI non nulles --%
%--                               --%
%-----------------------------------%

PA.A=[-c/m -k/m ; 1 0]; %-- Y=[ q' ; q];
PA.B=[1/m;0];
PA.alpha=0.; %-- facteur pour la non linéarité
PA.f0=1.; %-- fréquence de l'excitation
PA.u=@(t,PA)0; %-- terme d'excitation nul
Yo=[0 ; 1]; %-- condition initiale en déplacement non nul
dt=0.01;
Tspan=[0 : dt : 10]; %-- echantillonnage à pas constant
dt=diff(Tspan(1:2));
[t,Y]=ode45(@deriv,Tspan,Yo,[],PA);

figure;
plot(t,Y(:,2));
title('1DDL')
xlabel('Temps [s]')
ylabel('Amplitude [m]');

figure;
plot(t,Y(:,1));
title('1DDL')
xlabel('Temps [s]')
ylabel('Vitesse [m.s^-^1]');

Yf=fft(Y);
N=size(Y,1);
freq=0 : 1 /(N*dt) : (N-1)/N/dt;

figure;
plot(freq,(abs(Yf(:,2))),'-o');

%%
%----------------------------------%
%--                              --%
%-- Influence de l'amortissement --%
%--                              --%
%----------------------------------%

% Damping study, store multiple results as a structure, then finalize plot
damp=[ 0.04 0.4 4 ];
C1=struct('X',[],'Y',[]); % prepare structure for result
N=1000;
for j1=1:length(damp)
    PA.A=[-damp(j1)/m -k/m ; 1 0]; %-- Y=[ q' ; q];
    fp=1/(2*pi)*sqrt(k/m);
    N=1000;
    Tspan=linspace(0,20/fp,N+1);
    Tspan=Tspan(1:N);
    dt=diff(Tspan(1:2));
    [t,Y]=ode45(@deriv,Tspan,Yo,[],PA);
    freq=0 : 1 /(N*dt) : (N-1)/N/dt;
    
    C1.X=t;
    C1.Y(:,j1)=Y(:,2); % store first column for each result in loop
end

% Finalize plot
figure;
plot(C1.X,C1.Y);
xlabel('x');
ylabel('y');
set(gca,'xlim',[0 max(Tspan)])

figure;
semilogx(freq,log10(abs(fft(C1.Y))));
xlabel('x');
ylabel('y');
set(gca,'xlim',[0 max(freq)/2])

%%
%--------------------%
%--                --%
%-- Reponse forcée --%
%--                --%
%--------------------%

PA.A=[-c/m -k/m ; 1 0]; %-- Y=[ q' ; q];
PA.B=[1/m;0];
PA.f0=1.;
PA.u=@(t)cos((fmin+(fmax-fmin)/2*t/10)*t); %-- terme d'excitation sinusoidal
Yo=[0 ; 0]; %-- condition initiale nulle
dt=0.01;
Tspan=[0 : dt : 10]; %-- echantillonnage à pas constant
[t,Y]=ode45(@deriv,Tspan,Yo,[],PA);

figure;
plot(t,Y(:,2));
title('1DDL')
xlabel('Temps [s]')
ylabel('Amplitude [m]');

figure;
plot(t,Y(:,1));
title('1DDL')
xlabel('Temps [s]')
ylabel('Vitesse [m.s^-^1]');

Yf=fft(Y);
N=size(Y,1);
freq=0 : 1 /(N*dt) : (N-1)/N/dt;

figure;
plot(freq,(abs(Yf(:,2))),'-o');
%%
%--------------------------%
%--                      --%
%-- reponse non linéaire --%
%--                      --%
%--------------------------%

PA.A=[-c/m -k/m ; 1 0]; %-- Y=[ q' ; q];
PA.B=[1/m;0];
PA.alpha=100.; %-- facteur pour la non linéarité
PA.f0=1.; %-- fréquence de l'excitation
PA.amp=1000;
PA.u=@(t,PA)PA.amp*sin(2*pi*PA.f0*t); %-- terme d'excitation nul
Yo=[0 ; 0]; %-- condition initiale en déplacement non nul
dt=0.01;
Tspan=[0 : dt : 100]; %-- echantillonnage à pas constant
dt=diff(Tspan(1:2));
[t,Y]=ode45(@deriv,Tspan,Yo,[],PA);

figure;
plot(t,Y(:,2));
title('1DDL')
xlabel('Temps [s]')
ylabel('Amplitude [m]');

figure;
plot(t,Y(:,1));
title('1DDL')
xlabel('Temps [s]')
ylabel('Vitesse [m.s^-^1]');

Yf=fft(Y);
N=size(Y,1);
freq=0 : 1 /(N*dt) : (N-1)/N/dt;

figure;
plot(freq,(abs(Yf(:,2))),'-o');

