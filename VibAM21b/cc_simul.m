function [out,out1] = cc_simul(varargin)

% SUPPORT FUNCTION FOR VibAM project. For other data see
% http://savoir.ensam.eu/moodle/course/view.php?id=1874 
%
% Supported commands are
%
%  <a href="matlab:cc_simul('tuto')">List available tutorials (class support)</a>
%  <a href="matlab:cc_simul('webtp1DOF')">Open TP1 (one degree of freedom)</a>
%  <a href="matlab:cc_simul('webtpModal')">Open TP2 (transfers in time/freq)</a>
%  <a href="matlab:cc_simul('webtpCorrelDB')">Open TP3 (param/id/correlation)</a>
%
%  uf = cc_simul('simul',uf)  % time simulation
%  cc_simul('plot',uf)        % display simulation results

%       Etienne Balmes
%       Copyright (c) 1990-2019 by SDTools, All Rights Reserved.
%       See cc_simul('cvs') for revision information

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

persistent UseGui
if isempty(UseGui); UseGui=1;end

CAM=varargin{1}; 
if ~ischar(varargin{1});
  obj=varargin{1};evt=varargin{2};CAM=varargin{3};carg=4;
end
[CAM,Cam]=comstr(CAM,1); 
if nargin>1;uf=varargin{2}; carg=3;else;uf=[];carg=2;end

% --------------------------------------------------------------------------
%% #Simul  Simulation de la reponse temporelle (partie 2 TP2) ----------------
if comstr(Cam,'simul')  

% Frequency response (continuous time model)
uf.f=(0:length(uf.t)-1)/length(uf.t)/diff(uf.t(1:2));
uf.H=qbode(uf.sys,uf.f*2*pi);  

% When noise is present do averaging
if ~isfield(uf,'noise'); uf.noise=0;end
if uf.noise~=0;   uall=uf.u0(:)*[1 1 1 1 1];  % if noise 5 averages
else uall=uf.u0; 
end
filt=[];if isfield(uf,'filt');filt=uf.filt;end
sys=uf.sys;
try; warning('off','Control:analysis:LsimUndersampled');end

%% Loop on input signals (column of uf.u0) - - - - - - - - - - - - - - - - -
uf.frame={};
for j1=1:size(uall,2);  u0=uall(:,j1);

if isempty(filt)&&exist('lsim','file') 
%% Without_filter - - - - - - - - - - - - - - - - - - - - - - -
  disp('Simulation sans filtre anti-repliement');
  [y,t]=lsim(sys,u0(:),uf.t(:));  
  
  y=y+randn(size(y))*uf.noise;
  u=u0+randn(size(u0))*uf.noise;
  
  % store result with a check on sampling
  i1=find(t==uf.t(2))-1;
  if norm(t(1:i1:end)-uf.t(:))/t(2)<sqrt(eps)
   frame=struct('X',t,'Y',[u y(1:i1:end,:)]); 
  else error('problem with subsampling'); 
  end

elseif exist('lsim','file') 
%% With anti-aliasing filter - - - - - - - - - - - - - - - - - - - - - - -
  disp('Simulation avec filtre anti-repliement');
  warning('off','Control:analysis:LsimUndersampled');
  
   [y,t]=lsim(series(sys,filt),u0,uf.t); 
   y=y+randn(size(y))*uf.noise;u=u0+randn(size(u0))*uf.noise;
   
   [u,t]=lsim(filt,u0,uf.t);
   %[y,t]=lsim(filt,y,uf.t);
   
   i1=find(t==uf.t(2))-1; % gestion du sous-echantillonage
   if norm(t(1:i1:end)-uf.t(:))/t(2)>sqrt(eps); 
       error('problem with subsampling'); 
   else; u=u(1:i1:end);y=y(1:i1:end);
   end
  frame=struct('X',t,'Y',[u y]);  
  warning('on','Control:analysis:LsimUndersampled');

 elseif isempty(filt) % No filter ode45 (much slower than ode45)

  sys.u=@(t)interp1(uf.t,uf.u0,t);
  [t,x]=ode45(@(t,y)ssEval(t,y,sys),uf.t,zeros(size(sys.a,1),1),odeset);   
  y=x*sys.c';
  frame=struct('X',t,'Y',[uf.u0 y]);  

 end % with or without anti-aliasing

 if ~isfield(uf,'frame'); uf.frame={};end
 uf.frame{j1}=frame;

end % look on input signals
out=uf;
sdtweb('_link','uf=cc_simul(''simul'',uf);cc_simul(''plot'',uf);','Simu')

% --------------------------------------------------------------------------
%% #Plot : display commands ------------------------------------
elseif comstr(Cam,'plot'); [CAM,Cam]=comstr(CAM,5); 

if isempty(CAM)
%%  Ceci est un prototype de generation de trace necessaires pour le TP.
% Vous etes censes le modifier pour generer les traces que vous estimerez
% necessaires

% Entree et sorties sur les differentes mesures

y=uf.frame{1}.Y;y(1,length(uf.frame)*2)=0;
for j1=2:length(uf.frame);y(:,j1*2+(-1:0))=uf.frame{1}.Y;end

u=y(:,1:2:end); y=y(:,2:2:end); 

t=uf.t; f=(0:length(uf.t)-1)/length(uf.t)/diff(uf.t(1:2));
N=length(t);

gf=figure(1);set(gf,'name','Time response'); % traces temporels

subplot(221);	plot(t,u);	xlabel('t [s]');ylabel('u');axis tight
subplot(223);	plot(t,y);	xlabel('t [s]');ylabel('y');axis tight

% fenetrage 
if ~isfield(uf,'window'); uf.window='None';disp('Using natural window');end
win=fe_curve(['window' uf.window],N);
u=u.*(win(:,ones(size(u,2),1)));  y=y.*(win(:,ones(size(u,2),1)));

subplot(222);	plot(t,u);	xlabel('t [s]');ylabel('u*w');axis tight
subplot(224);	plot(t,y);	xlabel('t [s]');ylabel('y*w');axis tight


Ruu=fft(u).*conj(fft(u));  Ruy=fft(y).*conj(fft(u)); Ryy=fft(y).*conj(fft(y));
H1=mean(Ruy,2)./mean(Ruu,2);  H2=conj(mean(Ryy,2)./mean(Ruy,2));

gf=figure(4);set(gf,'name','Frequency response');
if isfield(uf,'H')   
 %% cas : simulation
 subplot(411);	semilogy(f,abs(fft(y)));	
 xlabel('Frequency [Hz]');ylabel('FFT(y)');axis('tight');
 subplot(412);	semilogy(f,abs(fft(u)));	ylabel('FFT(u)');axis('tight');
 subplot(212); semilogy(f,abs([H1 H2 uf.H]));axis('tight');
 ii_plp(ii_pof(eig(uf.sys.a)/2/pi,3))
 legend('H1','H2','Exact'); 
 setlines;
 
 try % Clean link axes
  if sdtdef('verm')>900
   set(gca,'userdata',linkprop(findobj(gf,'type','axes'),'xlim'));
  end
 end
 
else % cas : mesure reelle
 subplot(221);	plot(f,db(fft(y)));	ylabel('FFT(y)');
 subplot(222); plot(f,db([H1 H2]));  legend('H1','H2')
 subplot(223); plot(f,real(H2./conj(H1))); ylabel('Coherence');
 subplot(224); plot(f,db(Ryy./Ruy)); ylabel('H_1^n');
end
iimouse
out=struct('H1',H1,'H2',H2,'COH',real(H1./H2));

if sdtkey('cvsnum','iimouse')>1.252&&UseGui&&sdtdef('verm')>900
cingui('objset',[1 4],{'@dock',{'name','Simu','arrangement',[1 2], ...
    'tileWidth',[.3 .7]}});drawnow;
end

elseif comstr(Cam,'hist')
%% #PlotHist : display freq/damp history ------------------------------------

Res=varargin{2}; % Get Results from call cc_simul('plotHist',Res);

% Generate curves to based on pole history
figure(1);
subplot(121);semilogx(Res.Range.val(:,1),Res.hpo',':+');
 axis('tight');
 set(gca,'ylim',[1000 3200]);
 xlabel(Res.Range.lab{1});
 ylabel('Frequency [Hz]')
subplot(122);plot(Res.hpo',Res.hda'*100,'+');
 set(gca,'xlim',[1000 3000],'ylim',[0 30]); 
 xlabel('Frequency [Hz]'); ylabel('Damping ratio [%]');
 
% Display modes at each design point (interactive with datatip)
 def=Res.def;
 def.LabFcn='sprintf(''%.0f Hz %.2f %%, alpha=%.2f, beta%.2f'',def.data(ch,1:4).*[1 100 1 1])';
 cf=feplot;cf.def=def; fecom('colordataevala');fecom('coloralpha');
 cc_simul('rangedatacursor'); % open data tip in figure(1)

%%
else; error('Plot%s unknown',CAM);
end

elseif comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);
%% #Script : basic simulations to be edited

if comstr(Cam,'dfrf')
%% #TutoDfrf : tutorial steps for VibAM course -2
  
%% Step1 Initialization

if ~exist('cc_simul','file');vishno19('path');end
model=cc_simul('LoadCCasSDT');format('short','e');
% reorient for easier view
model.Node(:,5:7)=model.Node(:,5:7)*basis([0 0 1],[1 0 0]);

cf=sdth.urn('feplot(2)');cf.model=model;fecom(';showFiPro');
fecom('curtab Model');
[SE,CE]=fe_case(model,'assemble -matdes 2 1 -SE');

% Talk gives 
% - short description of .Node, .Elt, .pl, .il, .Stack
sdtweb('_link','sdtweb node');sdtweb('_link','sdtweb elt');
sdtweb('_link','sdtweb pl');sdtweb('_link','sdtweb il');
sdtweb('_link','sdtweb Stack');

model
% - short discussion of DOF, K, Klab
SE
SE.Klab % the matrix labels

% Matrix topology
m=SE.K{1}; k=SE.K{2};figure(1);clf; spy(k);

% When do you have non zero terms ?
nnz(k)/prod(size(k)) % What does this represent ?

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Step2 direct frequency response 

%% Load discuss DOF extraction. sdtweb('mdof')
b=fe_c(CE.DOF,426.03)'  
% How is this related to equations ?
% Why does CE.DOF==426.03 does not give the good result ?
c=b'  % Collocated load (what is it ?)

fecom('colorfacew-edgealpha.1-alpha0');
fecom('shownodemark',CE.DOF(b~=0),'markersize',5,'marker','o','color','r')

% Preallocate structures to store the result
freq=linspace(1390,1420,10)'; % Small range around first resonance
def=struct('def',zeros(size(CE.T,1),length(freq)), ...
       'DOF',SE.DOF, ...
       'data',freq);
C1=struct('X',{{def.data,426.02}}, ...
    'Xlab',{{'Frequency','DOF'}}, ...
    'Y',(b'*def.def).', ...
    'name','Full');

%% Loop on frequencies
%  Comment the different strategies used in this loop
for j1=1:length(freq) 
 % - omega^2 m + K = Z rigidité dynamique
 Z=-(freq(j1)*2*pi)^2*m+(1+.01i)*k; % Why the .01i ?
 % Why never inv(Z) ? 
 if j1==1; fprintf('\nDirect using MATLAB LU \n');
  tic; q=Z\b;toc
 elseif j1==2; fprintf('\nWith steps using MATLAB LU\n');
  tic;[L,U,P,Q,D]=lu(Z);toc % P*(D\Z)*Q = L*U 
  % Z q = D * P' (L * U ) * Q' * q = b
  tic;q= Q*(U\ (L\ (P*(D\b))));toc
  % Stupidly slow : [L,U]=lu(Z);nnz(L)+nnz(U) 
 elseif j1==3;
  fprintf('\nAn other solver : INTEL MKL PARDISO\n')
  ofact('mklserv_utils -silent');
  tic; Zd=ofact(Z);toc
  tic; q=Zd\b;toc
 else; % Other steps using pardiso
   q=ofact(Z,b);
 end
 def.def(:,j1)=CE.T*q;
 C1.Y(j1)=c*q; 
end

%% Another example with eigs/fe_eig 
% Note : profile is better than tic/toc
profile clear;profile on;cc_simul('ScriptProfEigs');profile report


%% Visualize factor sparsity. 
%  What is the difference between the two calls ?
c0 = chol(k + 1e3*m);
[c1,f,p]=chol(k + 1e3*m);
figure(1);
subplot(121);spy(c0);xlabel(sprintf('Density %.2f %%',nnz(c0)/prod(size(c0))*100));
subplot(122);spy(c1);xlabel(sprintf('Density %.2f %%',nnz(c1)/prod(size(c1))*100));

%% Response at input DOF. Just take a look at what I computed
iicom('curveinit',C1);

% GO BACK TO SLIDES NOW
% ------------------------------------------------------------------------
%% Step3 Now subspace size + reduction

%% phase collinearity
figure(1);clf;plot(def.def(:,1),'.');
xlabel('Real part');ylabel('Imaginary part');

[T,s,v]=svd([real(def.def) imag(def.def)],0);
figure(11);clf;
q=v*s;
% what does this frequency evolution represent ? (later slide SVD & variants)
plot(def.data,abs(complex(q(1:size(def.def,2),:),q(size(def.def,2)+1:end,:)))')
xlabel('Frequency');ylabel('Amplitude');
setlines([],{'-','--','-.'},'+ox*sdv^><ph')

% What is the meaning of those shapes ?
s=diag(s)/s(1)
cf.def=struct('def',T,'DOF',def.DOF,'data',s);fecom('ShowFiMdef')

%% Choose a significant norm 
T=[real(def.def) imag(def.def)];[u,s]=svd(T'*m*T);T=T*u';s=diag(s)/s(1)

cf.def=struct('def',T,'DOF',def.DOF,'data',s);fecom(';ShowFiMdef;view3')
% Why are the singular values different ? does this make more sense ?


% ------------------------------------------------------------------------
%% Step4 Rayleigh-Ritz (Galerkin)

% Phase1 : learning here by DFRF
% Phase 2 : basis building : here SVD
T=[real(def.def) imag(def.def)];[u,s]=svd(T'*m*T);T=T*u';s=diag(s)/s(1)
T=T(:,s>1e-8);
% Phase 3 : reduce model
mr=T'*m*T; kr=T'*k*T; br=T'*b; cr=c*T;
RO.T0=T;

%% Phase 4 : exploit reduced model : here reduced DFRF
C2=C1;C2.name='Reduced';C2.X{1}=linspace(freq(1),freq(end),100)';
dr=struct('def',zeros(size(mr,1),length(C2.X{1})),'DOF',(1:size(mr,1))'+.99);
tic;
for j1=1:length(C2.X{1})
 s=2i*pi*C2.X{1}(j1);
 qr= (mr*s^2+(1+.01i)*kr)\br;
 dr.def(:,j1)=qr;
 C2.Y(j1,1)=cr*qr;
end
toc
% Show forced response / full displacement
iicom('curveinit',{'curve','Full',C1;'curve','Red',C2});setlines

% Questions 
%  - what did we gain compared to DFRF ?
%  - how can we gain even more time ?

% GO BACK TO SLIDES NOW
% ------------------------------------------------------------------------
%% Step5 Krylov, orthog, ... Lanczos

%% Base krylov : what is wrong here ?
alpha=-(1300*2*pi)^2;kd=ofact(k+alpha*m); 
T=zeros(size(k,1),10);bc=b;
for j1=1:size(T,2)
  T(:,j1)=kd\bc; % K-1 bc
  bc=m*T(:,j1); % Update load
end
ofact('clear',kd);s1=svd(T);s1=s1/s1(1)


%% Krylov with load orthogonalization
kd=ofact(k+alpha*m);
T=zeros(size(k,1),10);
bc=b;
for j1=1:size(T,2)
  qc=kd\bc; % K-1 bc
  qc=qc/(qc'*m*qc); %normalize
  T(:,j1)=qc;
  bc=m*qc; % Update load
  % Renormalize 
  % Orthogonal to existing  Kq = bc - M T qr with T' K q =0
  ind=1:j1-1;
  if ~isempty(ind)
   MT=m*T(:,ind); 
   % Why is a pinv needed here ?
   bc=bc-MT* (  pinv(T(:,ind)'*MT)*(T(:,ind)'*bc));
   bc=bc/norm(bc);
  end
end
ofact('clear',kd);s2=svd(T);s2=s2/s2(1)

% What did we inprove here ?
disp([s1 s2])

cf.def=fe_eig({m,k,T,SE.DOF},[5 25 1e3]);

% ------------------------------------------------------------------------
%% Step6 sparse reduced model

i1={1:length(SE.DOF)};% Split in multiple segments along length
x=[-167 -89 0 89];
for j1=1:length(x)
 i2=i1{end}(fe_c(SE.DOF(i1{end}),feutil('findnode x<=',SE,x(j1)),'ind'));
 i1{end+1}=setdiff(i1{end},i2);i1{end-1}=i2;
end

T=repmat(RO.T0,1,length(i1));
for j1=1:length(i1);
 T(horzcat(i1{[1:j1-1 j1+1:end]}),j1*size(RO.T0,2)+[-size(RO.T0,2)+1:0])=0;
end
TR=struct('def',T,'DOF',SE.DOF);
figure(1);clf;spy(T'*k*T); % Sparsity of reduced model

% Reduced modes : which is exact
dr=fe_eig({m,k,T,SE.DOF});dr.TR=TR;
cf.def=dr;


%% EndTuto

elseif comstr(Cam,'profeigs')
%% #ProfEigs : illustrate profile within a file

% Get variables from caller
m=[];k=[];b=[];freq=[];CE=[];C1=[];def=[];c=[];
eval(iigui({'freq','m','k','b','CE','C1','def','c'},'GetInCaller')) 

% Direct Matlab
fprintf('\nDirect MATLAB using eigs\n')
[phi,ome]=eigs(k,m,20,'sm');
% Same call with the SDT layer
fprintf('\nIndirect using SDT call to eigs\n')
def=fe_eig({m,k,[]},[6 20 1e3]);

else; error('Script%s',CAM);
end

% --------------------------------------------------------------------------
%% #GUI : unused view for time simulation ------------------------------------
elseif comstr(Cam,'gui'); [CAM,Cam]=comstr(CAM,4);  % Prototype FEPLOT GUI

if comstr(Cam,'use')
  UseGui=1;
elseif isempty(Cam)
    
out=struct('X',{{[],{'u';'y'}}},'Xlab',{{'Time','Channel'}},'Y',[]);

out.Edit=struct('tmin',cingui('ParamEdit','tmin','0','%g','start of computation time range'),...
 'tmax',cingui('ParamEdit','tmax','0','%g','end of computation time range'),...
 'Fs',cingui('ParamEdit','fs','1','%g','Sampling frequency'),...
 'N',cingui('ParamEdit','N','1024','%g','Number of samples'),...
 'Recompute',struct('type','push','value','Recompute','callback', ...
      'cc_simul(''GuiRecompute'');','Prop',{{'Color',[1 0 0]}}));

ci=sdt.urn('iiplot(3)');ci.Stack{'curve','time'}=out;iicom('curtabStack')

elseif comstr(Cam,'recompute')
    keyboard
else error(1);
end

%% #Load : clean load of various tests ----------------------------------------
%% #LoadTest load test configuration from cc_test.mat file -2
elseif comstr(Cam,'loadtest')  

% if you need to reload the test FRFs, just copy these lines
if ~exist('cc_test.mat','file')
  addpath(comgui('cd o:/sdtdata/archive/cc_demo'));
end

% Load test data and clean details in the original file
r1=load('cc_test.mat');
UFS=r1.UFS; if isa(UFS,'sdth');UFS=UFS.vfields{1};end
UFS.header='';
if isfield(r1,'TEST');TEST=r1.TEST;else;TEST=r1.model;end
clear r1

r1=[];
if ishandle(2);
  ci=get(2,'userdata'); 
  if isa(ci,'sdth');r1=ci.Stack{'Test'};else;delete(2);end
  if ~isempty(ci.Stack{'IdFrf'});stack_rm(ci,'curve','IdFrf');end
end
if isfield(r1,'xf')&&isequal(r1.xf,UFS.xf)
 iicom(ci,';submagpha;ylog');%clear UFS;
 iicom('iixonly',{'Test','IdFrf'});
else
 try;[ci,XF]=idcom;ci.Stack{'curve','Test'}=UFS;
 catch;idcom;XF(1)=UFS(1);
 end % SDT 6.0 ou 5.3
end
iicom(ci,';submagpha;ylog');%clear UFS;
ci.IDopt.DataType='vel./force';

if sdtkey('cvsnum','idcom')>1.233&&~isdeployed
  iicom(ci,'dockid'); cf=idcom(ci,'feplot');
else
 cf=idcom(ci,'feplot');
end
if isempty(strfind(Cam,'-nowire'));cf.model=TEST; fecom('undefline');end
out=ci; out1=cf;

% -----------------------------------------------------------------------
%% #LoadCorrel load test configuration from cc_test.mat file and  -2
% FEM model from cc+plaq.mat
elseif comstr(Cam,'loadcorrel')  

% Open Test and wire and return pointers ci to iiplot and cf to feplot
[ci,cf]=cc_simul('loadtest-nowire');

% Provide a first mode estimate if nothing existing yet
r1=ci.Stack{'IdMain'};
if ~isfield(r1,'po')||isempty(r1.po);
 idcom('e .01 1415');idcom('ea'); 
end

% Load FEM model and define topology correlation
cg=feplot(10); 
fname=which('cc+plaq.mat');if isempty(fname);error('file not found');end
cg=fecom(cg,'load',fname);
r1=load('cc_test.mat','correl');
cg.mdl=fe_case(cg.mdl,'SensDof','Test',r1.correl);
cg.sel='reset'; fecom('showFiCevalY');

% Define parameters to be updated
upcom(cg.mdl,'ParStackAddK 50 .1 1000','Cover/base joint','group2');
upcom(cg.mdl,'ParStackAddK 5 .1 10','Screws','group1');

out=ci; out1=cg;

%% #LoadCCAsSDT
elseif comstr(Cam,'loadccassdt')  

fname=which('CCasSDT.mat');
if isempty(fname)
 model=feutil('rmfield',load(which('cc+plaq.mat')),{'mind','ke','me','st', ...
     'ki','mi','wd','DOF','Opt','file','copt','bas'});
 i1=ismember(model.il(:,1),2:7);model.il(i1,7:end)=0;
 model.il(i1,3)=5; % Formulation
 model.il(i1,4)=1e3; % drilling
 model.il(1,3)=30; %beam torsion 
 i1=feutil('findelteltname celas',model);model.Elt(i1,3:4)=-abs(model.Elt(i1,3))*[1 0];
 %model=fe_case(model,'fixdof','drill',[2003.04;2003.05;2003.06;2190.04;2190.05;2190.06]);
 model.il=p_solid(model.il,'dbval 1 d3 -3','dbval 8 D3 -3');
 model=orderfields(model,{'Node','Elt','pl','il','Stack'});
 fname=fullfile(fileparts(which('cc+plaq.mat')),'CCasSDT.mat');
 save(fname,'model','-v7'); % v7 for Python
else;
 load(fname);
 fprintf('Loaded model from ''%s''\n',fname);
end
out=model;

% -----------------------------------------------------------------------
%% #LoadDBFE load drum brake model with cable guide -2
elseif comstr(Cam,'loaddbfe')  

fname=cc_simul('wd','DbCoarseAssembly.mat');

model=[];load(fname,'model')
%model.Node(:,5:7)=model.Node(:,[7 5 6]);
model=fe_case(model,'parAddK','Kv','matid2', ... % Viscoelastic material
   'parAddM','Mg','MatId 2 3'); % Mass of cable guide
[Up,CE]=fe_caseg('assemble NoT -SE -matdes 2 1 -1.1',model);
Up=stack_set(Up,'case','Case 1',CE);

ci=comgui('guiiiplot -reset',2);
cf=comgui('guifeplot -reset',3);
cf.model=Up;

out=ci;out1=cf;
% -----------------------------------------------------------------------
%% #LoadDBTest load drum brake test [ci,cf]=cc_simul('loadDBTest'); -2
elseif comstr(Cam,'loaddbtest')  

% D:\balmes\Dropbox\ENSAM\savoir_TP\VibAM\J92_coarse_assembly.mat
fname=cc_simul('wd','DbTest.mat');load(fname,'FRF','ID','TEST');
% addpath d:/balmes/dropbox/tex/id
cg=feplot(5);
if isempty(cg.mdl.Elt); cg.mdl=TEST;end

ci=iiplot(2);iicom('curveinit','Test',FRF)
r1=ci.Stack{'IdMain'};
if isempty(r1.po) % Initialize result
 ci.Stack{'IdMain'}=struct('po',ID.po(ID.po(:,1)<2500&ID.po(:,1)>1000,:),'res',[]);
 idcom('estlocalpole');
 
 cg.def=ci.Stack{'IdMain'}; % Display mode
 cg.sel={'-Test','colordataEvalA-edgealpha.1'};

end
out=ci;out1=cg;

% -----------------------------------------------------------------------
%% #LoadDBCorrel load drum brake test/correl/FEM -2
%[ci,cf]=cc_simul('loadDBCorrel');
elseif comstr(Cam,'loaddbcorrel')  

% Now display and possibly compute model
fname=cc_simul('wd','DbTest.mat');load(fname,'SE','TEST');
SE=fe_case(SE,'stack_set','SensDof','Test',TEST);
ci=iiplot(2);if isempty(ci.Stack{'Test'});cc_simul('LoadDbTest');end
    
cf=feplot(3);
if sdtdef('cvsnum','feplot')>=1.514 % Use dock now
 fecom(cf,'DockCoShape');
 ii_mac(1,'setMAC',struct('da',ci.Stack{'IdMain'}, ...
     'sela',{{'-Test','ShowFiTestDef'}}));
end
if (isfield(cf.mdl,'il')&&cf.mdl.il(1,3)==-1.2)||isempty(cf.def)|| ...
      isfield(cf.def,'TR');
 MAP=feutil('getnormal node MAP',SE);
 mo1=SE;mo1.Node(:,5:7)=mo1.Node(:,5:7)-MAP.normal*-1.2;mo1.il(3)=-1.2;
 def=fe_eig(mo1,[5 20 1e3]);
 feplot(cf,mo1,fe_def('subdef',def,def.data>1));
 cf.sel={'-linface','colordataEvalA-edgealpha.1'};
end
if sdtdef('cvsnum','feplot')>=1.514 % Use dock now
 ii_mac(1,'SetMac',struct('db',def)) % Set FEM shape
end
%% 

out=ci;out1=cf;
% -----------------------------------------------------------------------
%% #RangeCh : change feplot channel on datatip position
elseif comstr(Cam,'rangech')

if length(dbstack)>20; return; end % Cancel update
gf=3;
if ~ishandle(gf)||~strcmp(get(gf,'tag'),'feplot');
  warning('Expecting model to be displayed with cf=feplot(3);cf.model=model;')
  return
end
cf=get(gf,'userdata'); 
r1=real(cf.def.data); 
if isfield(evt,'pos');pos=evt.pos;else; pos=evt.Position;end
ch=find(r1(:,1)==pos(2)&r1(:,3)==pos(1)); % Left freq/alpha
if isempty(ch)
 ch=find(r1(:,1)==pos(1)&r1(:,2)*100==pos(2)); % right zeta/freq
end
if ~isempty(ch);fecom(cf,sprintf('ch%i',ch));end
figure(1);

% -----------------------------------------------------------------------
%% #RangeDataCursor : allows cursor picking
elseif comstr(Cam,'rangedata')  

figure(1);
subplot(121);
go=findall(gca,'type','line');go=go(1);
RO=struct('X',1,'DataTipFmt','out=sprintf(''alpha=%.2f\n%.2f Hz'',pos);', ...
    'PostFcn',{{@cc_simul,'RangeCh'}});
try;iimouse('datatipnew',go,RO);end

subplot(122);
go=findall(gca,'type','line');go=go(1);
RO=struct('X',1100,'DataTipFmt','out=sprintf(''%.2f Hz\n%.2f %%'',pos.*[1 1]);', ...
    'PostFcn',{{@cc_simul,'RangeCh'}});
try;iimouse('datatipnew',go,RO);end
 
%% #admin -------------------------------------------------------------------
elseif comstr(Cam,'web')
 % cc_simul('WebtpModal');
 fname=fullfile(fileparts(which('cc_simul')),'help',[CAM(4:end) '.html']);
 if ~exist(fname,'file')
 fname=fullfile(fileparts(which('cc_simul')),'html',[CAM(4:end) '.html']);
 end
 if ~exist(fname,'file');fname=strrep(fname,'.html','e.html');end
 fname
 web(fname);
 
elseif comstr(Cam,'wd')
 fname=uf;
 out=which(fname); 
 if isempty(out);out=fullfile(fileparts(which('cc_simul.m')),fname);end
 if isempty(out)||~exist(out,'file');
     out=fullfile(sdeb('wdENSAM\savoir_TP\VibAM\'),fname);
 end
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','cc_simul','CAM',CAM)));
 if nargout==0; clear out; end

elseif strcmp(Cam,'cvs')
    out='$Revision: 1881 $  $Date: 2020-11-04 17:59:16 +0100 (Wed, 04 Nov 2020) $'; return; 
else;error('%s unknown',CAM);
end % simulation ou trace


%% ssEval : ode45 based evaluation
function dx=ssEval(t,x,sys);

dx=sys.a*x+sys.b*sys.u(t);
