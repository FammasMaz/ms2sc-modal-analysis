function [out,out1]=mevib16(varargin);

% Fichier de support TP de mécanique vibratoire

%       Etienne Balmes
%       Copyright (c) 2016-2016 by SDTools, All Rights Reserved.
%       $Revision: 1.13 $  $Date: 2016/11/17 13:51:57 $

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
if nargin<1; CAM='init'; Cam='init';
elseif ischar(varargin{1}); 
  obj=[];evt=[];[CAM,Cam]=comstr(varargin{1},1); carg=2;
else;
 obj=varargin{1};evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1); carg=4;
end



if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);
%% #Script : basic simulations to be edited

if comstr(Cam,'tube')
 %% #ScriptTube

 model=mevib16('MeshTube',struct('lx',1500));
 %% Modes of closed tube
 d1=fe_eig(model,[5 30 1e3]);
 cf=feplot(2);cf.model=model;mevib16('ViewFluidVel',d1);fecom('ch2');
 
 %% Frequency response (test impedence)
 model=stack_set(model,'info','Freq',linspace(50,2000,100)');
 
 % Paroi directe R=16, R=3.2 mélamine 4cm, R=1.57 mélamine 10 cm
 
 model.pl(2,6)=3.2; % Real part of reduced impedance (high <=> reflective)
 FRF=fe_simul('dfrf',model);
 FRF.LabFcn=['sprintf(''%i @ %.1f Hz, R=' sprintf('%.2f',model.pl(2,6)) ...
     ''',ch,def.data(ch,1))'];
 FRF.ID={struct('po',d1.data(:,1)*[1 0])};%Frequencies of closed tube
 
 iicom('curveinit','P',FRF);iicom(';submagpha;ylog'); 
 iicom('ch 10 20 30 40');
 mevib16('ViewFluidVel',FRF);
 fecom('CursorOnIiplot'); % On click change sensor

 
 %% pressure on the central line (use quadradic)
 mo2=model;
 alpha=@(x)4*x./((x+1).^2);
 R1=@(a)(-2*(a-2)-sqrt(4*(a-2).^2-4*a.^2))./2./a;
 R2=@(a)(-2*(a-2)+sqrt(4*(a-2).^2-4*a.^2))./2./a;

 mo2.pl(2,6)=10; % Real part of reduced impedance
 d1=fe_simul('dfrf',stack_set(mo2,'info','Freq',[500 1000 1300]));
 mevib16('ViewFluidVel',d1)
  
 n1=sortrows(feutil('getnode y==0 & z==0 & x>100',mo2),5);
 C2=struct('X',{{d1.data,n1(:,5)}},'Xlab',{{'Freq','Pos'}}, ...
     'Y',(fe_c(d1.DOF,n1(:,1)+.19)*d1.def).');
 figure(1);clf; 
 subplot(211);plot(C2.X{2},real(C2.Y));%mesh(C2.X{1},C2.X{2},abs(C2.Y)')
 subplot(212);plot(C2.X{2},abs(C2.Y));%mesh(C2.X{1},C2.X{2},abs(C2.Y)')
 legend(cellstr(num2str(C2.X{1}(:,1))));axis('tight');
 xlabel('Position');
 
 % TOS = p_min/p_max
 TOS= min(abs(C2.Y),[],2)./max(abs(C2.Y),[],2);
 title(sprintf('Alpha %.3g',mean(alpha(TOS))))

 [{'R1','R2','Alpha'};num2cell([R1(alpha(TOS)) R2(alpha(TOS)) alpha(TOS)])]

 
 % Check velocity of sound by distance between min/max
 ld=330e3/1000;% Wavelength in mm


 %% Periodic computation
 
 mo2=model;mo2.Elt=feutil('selelt withnode{x==}',mo2,max(mo2.Node(:,5)));
 mo2.Node=feutil('getnode groupall',mo2);
 mo2.Node(:,5)=mo2.Node(:,5)-min(mo2.Node(:,5)); % Scale to 1 mm
 mo2.Node(:,5)=mo2.Node(:,5)/max(mo2.Node(:,5))*10; % Scale to 1 mm
 mo2=stack_set(mo2,'info','EigOpt',[5 50 1e3]);
 mo2=fe_cyclic(sprintf('build -1 %.15g 0 0',max(mo2.Node(:,5))),mo2);
 par={'lab "ncx" min 100 max 2.01 scale "ilin" NPoints 30'};
 RD=fe_range('Grid',par); 
 [def,hist]=fe_homo('dftDisp',mo2,RD);
 figure(1);plot(hist.Y'/1000,hist.X{1}(:,1));axis('tight');
 xlabel('Frequency [kHz]');ylabel(hist.Xlab{1});
 set(gca,'xlim',[0 6])
 
elseif comstr(Cam,'beam')
 %% #ScriptBeam (atelier 2-3)
 % Freq, sweep, lissajou, added mass, Proof mass damper, modeshapes
 % Identification ? (-3dB), Impact, RMS responses
 
 model=mevib16('MeshBeam1',struct('lz',5.8)); % Acier inox
 d1=fe_eig(model,[5 10 0]);
 feplot(model,d1);fecom('colordataEvalZ');
 
 % Radiaflex page 61
 % 20 dAN (approx mg) / 3 mm 
 k=(15*2/5e-3) % N/m
 m=.1; sqrt(k/m)/2/pi
 % Buil references 521293/521128/521295

 model=mevib16('MeshBeam1',struct('PM',200e-3,'PS',6000));
 d2=fe_eig(model,[5 20 0]);[d1.data d2.data (d2.data./d1.data-1)*100]
 
 feplot(model,d2);fecom('colordataEvalZ');
 
 
 %%  Now build transfer
 
 model=stack_set(model,'info','Freq',linspace(1,300,500)');
 nor2xf(d1,.01,model,'iiplot "base" -po');
 nor2xf(d2,.01,model,'iiplot "damper" -po');
 iicom('submagpha');
 
 
elseif comstr(Cam,'plate')
 %% #ScriptPlate : forme des modes + masse ajoutée => impact formes

 sdtweb mevib16 MeshPlate % Meshing script
 
 % Start computation from script
 mevib16('SetMv',struct('PNode',9,'PM',.001,'PModes','do'));
 
 mo1=feutil('addelt',model,'mass1',[22 .01 .01 .01]);def=fe_eig(mo1);
 mo1=feutil('addelt',model,'mass1',[40 .005 .005 .005]);def=fe_eig(mo1);

 %model=feutil('addelt',model,'mass1',[22 .01 .01 .01]);def=fe_eig(model);
 def=fe_eig(model);

 feplot(model,def);fecom('colordataEvalZ');
 ii_plp('colormap wcenter thres .05',jet(50))
 %fecom('colormap',jet(7))
 
 
 d3=def;d3.def=d3.def(:,7:8)*[1;1i];d3.data=d3.data(7);cf.def=d3;
 
    
elseif comstr(Cam,'brake')
 %% #ScriptBrake
 %% Focus on cable guide resonance numerical side 
 [ci,cf]=mevib16('loadBrake');

 i1= [23199 23283 5161 5170]';
 nor2xf(def,.01,i1(1)+.03,i1+.03,linspace(10,3000,500)','acc Hz -po iiplot"Test"');
 iicom(';iix:TestOnly;submagpha'); iicom('ch1:4')

 C1=ci.Stack{'Test'};
 imode=[12 14 15 18];
 cp=fe_c(def.DOF,i1+.03)*def.def(:,imode);r2=[0 1 0 0]/cp;%(cp')\[0;1;0];
 r1=[C1.Y*r2(:) C1.Y];
 iw=find(C1.X{1}>500&C1.X{1}<2500);r1=r1*(diag(1./max(abs(r1(iw,:)))));
 figure(1);h=plot(C1.X{1},abs(r1));ii_plp(def.data(imode)*[1 .0]);
 set(h(1),'linewidth',3)
 
%% Modal observation on test
 mevib16('LoadTxt','Plateau/*.txt')
 ci=iiplot; idcom(ci,'ShowTest');
stack_set(ci,'curve','IdMain',struct('po',[
%{'Freq','Damp','Index'}
890.89899 0.007955 % 1 
1160.01668 0.009422 % 3 
1229.81500 0.010880 % 4 
]));
iicom('wmin 850 1293');idcom('estlocalpole');
 
 Id=ci.Stack{'IdMain'}; % Modes
 C1=ci.Stack{'Test'}; % Réponse
 cp=Id.res(1:size(Id.po,1),:).'; 
 r2=[0 0 1]*pinv(cp);r2=r2/norm(r2)*5;
 
 mevib16('ViewBrFilt',r2);
 
 inw=find(C1.X{1}>500 & C1.X{1}<1500);
 xf=C1.Y(inw,:); xf=[xf*r2(:) xf];
 figure(1);h=plot(C1.X{1}(inw),abs(xf)*diag(1./max(abs(xf))));
 set(h(1),'linewidth',3);set(h(2:end),'linestyle',':')
 axis('tight');
 xlabel(C1.Xlab{1});ii_plp(ci.Stack{'IdMain'});
 
 
else; error('Mesh%s unknown',CAM);
end

elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);

if carg>nargin;RO=struct;else;RO=varargin{carg};carg=carg+1;end

%% #Mesh : meshing/loading commands -------------------------------------    
if comstr(Cam,'tube')
 %% #MeshTube : 106 mm inner diameter. L 530 or 1500 thick tube 2 mm

 R1=struct('lx',530,'d',106);
 % Merge non given fields
 RO=feval(sdtroot('@sfield'),'AddMissing',RO,R1);

 node=[0 0 0;RO.lx 0 0;0 0 RO.d/2]; 
 model=feutil('object quad 1 1',node,50,2);
 model=feutil('rev 20 o 0 0 0  360  1 0 0',model);
 
 % End surface to support surface impedance
 elt=feutil('selelt selface & innode {x==}',model,max(model.Node(:,5)));
 model=feutil('addelt',model,feutil('set group 1 matid 2 proid 2',elt));
 
 % Material
 model.pl=m_elastic('dbval 1 Air -unit MM','dbval 2 Air -unit MM');
 model.pl(2,6)=.1; % Default real part for reduced wall impedance
 % Integration (element formulation properties)
 model.il=p_solid('dbval 1 d3 -3', ...
      'dbval 2 wallimp'); % acoustic impedance on a wall
 model.pl(:,4)=340e3; % Velocity at 340 m/s
 model.unit='MM';model.name='Kundt';

 % Add source
 model=fe_case(model,'DofLoad','IN',struct('def',1,'DOF',1.19));
 model.Node(model.Node(:,1)==1,6)=10; % sligtly off axis
 % No damping in fluid
 model=stack_set(model,'info','DefaultZeta',0); 
 ofact('spfmex;');ofact('silent;'); % SDT options

 
elseif comstr(Cam,'beam');[CAM,Cam]=comstr(CAM,5);
%% #MeshBeam
if isfield(RO,'Blx')&&ischar(RO.Blx);RO.Blx=str2num(RO.Blx);end
if ~isfield(RO,'BK');RO.BK=0;elseif ischar(RO.BK);RO.BK=str2num(RO.BK);end % No mass

PA=sdtroot('paramvh');
if isfield(PA,'Mv');Mv=PA.Mv;else;Mv=[];end

if any(RO.Blx==[600 605]) 
 %% Poutre 1 (acier inox) pot 320, spring 480
 R1=struct('ly',40,'Blz',6,'holes',[320 400 480 592],'in',320, ...
     'spring',400,'out',592);
 if ~isfield(RO,'BM')||RO.BM==0;
  RO.BM=.154;
 end
else
 %% Poutre 1 (ancier basique) 530
 R1=struct('ly',30.2,'Blz',4.9,'holes',[232 298 364 430], ...
       'in',298,'out',298,'spring',364);
 if ~isfield(RO,'BM')||RO.BM==0;
  RO.BM=.374; RO.BK=50e3;
 end
 % RO.BM
end
% Merge non given fields
RO=feval(sdtroot('@sfield'),'AddMissing',RO,R1);

node=[0 -RO.ly/2 -RO.Blz;RO.Blx 0 0;0 RO.ly 0; 0 0 RO.Blz];
if remi(RO.Blx,1)==0 % 530 or 600
  r2=feutil('refineline 10',[0 RO.holes RO.Blx])/RO.Blx;
else; % Shift 5 mm
  node(1)=node(1)-5; 
  r2=feutil('refineline 10',[0 RO.holes+5 RO.Blx])/RO.Blx;
end
model=feutil('objecthexa 1 1',node,r2,2,1);
model=feutil('lin2quad',model);

model.unit='MM';
model.pl=m_elastic('dbval 1 Steel -unit MM');model.pl(3)=200e6;
model=p_solid('default;',model);
model=fe_case(model,'FixDof','Base',sprintf('x==%g',min(model.Node(:,5))));

n1=feutil(sprintf('getnode x==%g & y==0 & z==0',RO.spring),model);
n2=feutil(sprintf('getnode x==%g & y==0 & z==0',RO.out),model);
n3=feutil(sprintf('getnode x==%g & y==0 & z==0',RO.in),model);
model=fe_case(model,'DofLoad','In',struct('def',1,'DOF',n3(1)+.03), ...
      'SensDof','Out',n2(1)+.03);

%% Deal with proof mass
if RO.BM<=0 % Nothing
elseif RO.BK==0 % Add mass (but no stiffness since m=0)
  model=feutil('addelt',model,'mass1',[n1(1) RO.BM*[1 1 1 0 0 0]]); 
else % Add proof mass with spring
  model.Node(end+1,[1 5:7])=[max(model.Node(:,1))+1 n1(1,5) 0 2];
  model=feutil('addelt',model,'mass1',[model.Node(end,1) RO.BM*[1 1 1 0 0 0]]); 
  
  elt=n1(1,1);elt(:,2)=model.Node(end,1);
  elt(:,3)=123456;
  elt(:,7)=RO.BK; % Stiffness
  elt(:,10)=RO.BK*RO.Beta; % 
  model=feutil('addelt',model,'celas',elt);
  model=fe_case(model,'fixdof','RotMass',model.Node(end,1)+[1 4:6]'/100);
end

model=fe_case(model,'fixDof','Bend','y==0 -DOF 2');
if nargout==0; 
    d1=fe_eig(model,[5 20 0]);feplot(model,d1);
    fecom('colordataEvalZ');
end
if ~isempty(Mv);
 sdcedit(Mv,'BM',num2str(RO.BM));sdcedit(Mv,'Beta',num2str(RO.Beta));
 sdcedit(Mv,'BNode',num2str(n1(1)));
 sdcedit(Mv,'BK',num2str(RO.BK));
 ua=clean_get_uf('tab;','SDT Root','Mv');ua.JTable.repaint;
end  % Update int


elseif comstr(Cam,'plate')
%% #MeshPlate  [model,def]=mevib16('MeshPlate -cf2')
   
 % Geométrie
 model=femesh('TestQuad4 divide 17 17');
 RO.lx=.25;
 model.Node(:,5:6)=(model.Node(:,5:6)-.5)*RO.lx;
 
 % Propriété plaque
 model.il=p_shell('dbval 110 Kirchoff .0015'); % 
 % Propriété matériau 
 model.pl=m_elastic('dbval 100 Steel');
 
 % Conditions aux limites
 i1=feutil('findnode x>-.012 & x<.012 & y>-.012 & y<.012',model);
 model=fe_case(model,'fixdof','center',i1, ...
     'fixdof','Drill',.06);
 model=stack_set(model,'info','EigOpt',[5 20 0]); % Calcul 20 modes
 model.unit='SI';model.name='SquarePlate';
 
   
else; error('Mesh%s unknown',CAM);
end
out=model;

%% #MeshInitFigures
wd='o:/balmes/scratch/MeVib';if ~exist(wd,'dir');wd=sdtdef('tempdir');end
sdtroot('SetProject;',struct('ProjectWd',wd,'PlotWd',fullfile(wd,'plots')))
cf=comgui('guifeplot -Project"SDT Root"',2);
ci=comgui('guiiiplot -Project"SDT Root"',3);
cingui('plotwd',cf,'@OsDic(SDT Root)',{'FnMI'});

elseif comstr(Cam,'load');[CAM,Cam]=comstr(CAM,5);
%% #Load
 if comstr(Cam,'uff')
  %% #LoadUff
  li=dir(varargin{carg});carg=carg+1;
  if isfield(li,'folder');wd=li(1).folder;
  else;wd=fileparts(varargin{carg-1});
  end
  
  COH=[]; FRF=[];lab=cell(length(li),1);
  for j1=1:length(li)    
    r1=ufread(fullfile(wd,li(j1).name));   
    [un1,lab{j1}]=fileparts(li(j1).name);
    if j1==1;FRF=r1(1); if length(r1)>1;COH=r1(2);end
    else; 
      r2=r1(1);if size(r2.w,1)~=size(FRF.w,1);dbstack;keyboard;end
      FRF=fe_def('curvejoin',FRF,r1(1)); 
      if ~isempty(COH);COH=fe_def('curvejoin',COH,r1(1));end
    end
  end
  FRF=struct('X',{{FRF.w,lab}},'Xlab',{{'Frequency [Hz]','IO'}}, ...
      'Y',FRF.xf,'dof',FRF.dof);
  ci=comgui('guiiiplot -reset',3);
  if ~isempty(COH);stack_set(ci,'curve','COH',COH);end
  iicom('curveinit','Transfer',FRF);iicom('ch',{'IO',1:size(FRF.Y,2)});
  
  
 elseif comstr(Cam,'txt')
  %% #LoadTxt
  if carg>nargin
    [li,wd] = uigetfile('*.txt', 'Pick files','MultiSelect','on');
    if isequal(li,0);return;elseif ischar(li);li={li};end
  else
   li=dir(varargin{carg});carg=carg+1;
   if isfield(li,'folder');wd=li(1).folder;
   else;wd=fileparts(varargin{carg-1});
   end
  end
  
  FRF=[]; COH=[]; lab=cell(length(li),1);
  for j1=1:length(li)    
    if iscell(li);fname=fullfile(wd,li{j1});
    else; fname=fullfile(wd,li(j1).name);
    end
    fid=fopen(fname);st=fscanf(fid,'%c');fclose(fid);
    [un1,lab{j1}]=fileparts(fname);
    st(st==13)=10;st(st(1:end-1)==10&st(2:end)==10)='';
    i1=strfind(st,'Y_Phase');
    if ~isempty(i1)
     i1=i1(1);i2=i1+find(st(i1(1)+(1:200))==10);
     r1=textscan(st(i2(2):end),'%n%n%n');r2=r1{2}.*exp(1i*pi/180*r1{3});
     st='';i1=[];
    end
    if ~isempty(st);i1=strfind(st,'Y_ImagValue');end
    if ~isempty(i1)
     i1=i1(1);i2=i1+find(st(i1(1)+(1:200))==10);
     r1=textscan(st(i2(2):end),'%n%n%n');r2=complex(r1{2},r1{3});
     st='';i1=[];
    end
    if ~isempty(st);i1=strfind(st,'Real');end
    if ~isempty(i1)
     i1=i1(1);i2=i1(1)+find(st(i1(1)+(1:200))==10);
     st=strrep(st,',','.');
     r1=textscan(st(i2(1):end),'%n%n%n');r2=complex(r1{2},r1{3});
     st='';i1=[];
    end
    
    if isempty(FRF);FRF=struct('w',r1{1},'xf',r2);
    else; FRF.xf(:,end+1)=r2;
    end
  end
  FRF.dof=(1:size(FRF.xf,2))'+.01; iw=FRF.w>5;
  FRF=struct('X',{{FRF.w(iw),lab}},'Xlab',{{'Frequency [Hz]','IO'}}, ...
      'Y',FRF.xf(iw,:),'dof',FRF.dof);
  ci=comgui('guiiiplot -reset',3);
  iicom('curveinit','Test',FRF);
  iicom(';iix:TestOnly;w0');iicom('ch',{'IO',1:size(FRF.Y,2)});
  
 elseif comstr(Cam,'brake');
 %% #LoadBrake  -cf2 -ci3
     
  ci=comgui('guiiiplot -reset',3);cf=comgui('guifeplot -reset',2); 
  if ~isempty(strfind(Cam,'cf'))
   % mevib16('loadBrake -cf2'); % load mesh
   fname=mevib16('wd','DbCoarseAssembly.mat');
   model=[];def=[]; load(fname,'model','def')
   feplot(cf,'initmodel-back',model);
   %feplot(cf,'initdef-back',def); 
   cf.sel={'-linface','colordataEvalZ-edgealpha.1'};
   cf.def=def; fecom('ch7');
   ii_plp('colormap',jet(20),cf.ga);
   figure(2);fecom(cf,'renderer opengl');
   eval(iigui({'model','def'},'SetInBaseC')) 
  end
  if ~isempty(strfind(Cam,'ci')) % Load test
   fname=mevib16('wd','DbTest.mat');
   r1=load(fname,'FRF','ID','TEST');
   r2=fe_def('subdef',r1.FRF,r1.FRF.w>200);
   r2=fe_def('subdof',r2,[(1:2)'; fe_c(r2.dof(:,1),.03,'dof')]);
   stack_set(ci,'curve','Test',r2);
   stack_set(ci,'curve','IdMain',struct('res',[],'po',r1.ID.po(4:13,:)));
   iicom(';iix:TestOnly;submagpha;cax2;ylog;w0');idcom(ci,'est');
   cg=feplot(55);
   cg.model=r1.TEST;
   cg.def=ci.Stack{'IdMain'};fecom(cg,'colordataEvalZ');
   
  end
  out=ci;out1=cf;

 else; error('Load%s unknown',CAM);
 end

 
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
%% #View

%% #ViewFluidVel : mevib16('ViewFluidVel')
if comstr(Cam,'fluid')

 cf=feplot; 
 if carg>nargin;d1=cf.def;else;d1=varargin{carg};carg=carg+1;end
 if isempty(fe_c(d1.DOF,.01,'ind')) % Add velocity field
  mo1=cf.mdl.GetData; 
  mo1.Elt=feutil('selelt eltname hexa',mo1);
  C1=fe_stress('stress gstate',mo1,d1);%C1=C1.GroupInfo{1,5};
  d2=fe_stress('expand',mo1,C1);
  d2=struct('def',[reshape(d2.def,3*size(d2.def,1),[]);d1.def], ...
     'DOF',[feutil('getdof',(1:3)'/100,fix(d2.DOF));d1.DOF], ...
     'data',d1.data);
  cf.def=d2;
 end
 cf.sel={'innode {y>=0}','colordata19-edgealpha.1'};
 fecom colorscaleinstanttight;set(cf.ga,'climmode','auto');
%% #ViewMac : mevib16('ViewMac')
elseif comstr(Cam,'mac')
 ci=iiplot(3);
 if isempty(ci.Stack{'Test'});mevib16('LoadBrake-ci3');end  
 cf=feplot(2);if isempty(cf.mdl.Elt);mevib16('LoadBrake-cf2');end  
 model=cf.mdl.GetData;def=cf.def;
 sens=fe_case(model,'sens');
 figure(1);ii_mac(ci.Stack{'IdMain'},fe_def('subdef',def,7:20),'sens',sens,'macplot');
 ii_mac(1,'_plotmactext');ii_mac(1,'macErrorTable');
 
 
%% #ViewBrFilt : mevib16('ViewBrFilt') : fixed modal filter
elseif comstr(Cam,'brfilt')

    
ci=iiplot; C1=ci.Stack{'Test'};
if ~isfield(C1,'Y')||size(C1.Y,2)~=5
   error('On attend 5 mesures, dossier plateau'); 
end

if carg>nargin
 r2=[
  -8.0735e-01 + 1.7924e+00i % 'p1_ref'
  -5.7673e-01 + 8.7476e-01i % p2 support
   1.9819e+00 - 3.0634e-01i % p3 support
  -2.5766e-01 - 2.8997e+00i % p4 fond
   6.5167e-02 - 2.7453e+00i % p5 fond
 ];
else; r2=varargin{carg};carg=carg+1;
end
 inw=find(C1.X{1}>500 & C1.X{1}<1500);
 xf=C1.Y(inw,:); xf=[xf*r2(:) xf];
 figure(1);h=plot(C1.X{1}(inw),abs(xf)*diag(1./max(abs(xf))));
 set(h(1),'linewidth',3);set(h(2:end),'linestyle',':')
 axis('tight');
 xlabel(C1.Xlab{1});ii_plp(ci.Stack{'IdMain'});
 

    
else;error('View%s unknown',CAM);
end
%% #admin
elseif comstr(Cam,'init')
%% #Init mevib16('init')

UI=sdtroot('paramUI');PARAM=sdtroot('paramvh');
if ~isfield(UI.DefBut.Tab,'Mv')||~isfield(PARAM,'Mv')
 UI.DefBut.Tab.Mv=struct('name','Tab.Mv','ParamField','Mv', ...
    'InitFcn',{{'mevib16','Init'}},'SetFcn',{{'mevib16','Set'}});
 sdtroot('paramui',UI);
 mevib16('setMv');
 sdtroot('initMv'); return;
end

  r2=PARAM.Mv;CAM='Mv';
  ua=struct('name',CAM,'IntegerHandle','off','ToolTip',CAM);
  ua=sdt_dialogs('uatable-end0 -tipascol 1',ua,r2);
  % tab data implementation
  ua.NeedClose=1;ua.ParentPannel=UI.gf; ua.ExploName=['PTree.' CAM];
  % store data in java pointer, then display the tab
  [ub,ga]=cinguj('TabbedPaneAdd',UI.gf,ua); 
  set(UI.gf,'visible','on');



elseif comstr(Cam,'set')
%% #Set mevib16('set')

 %% #SetMain - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -2
% mevib16('setMv Reset');mevib16('init')
% if comstr(Cam,'main'); [CAM,Cam]=comstr(CAM,4+1);
%  [r1j,r1,st]=defaultSet(UI,PARAM,RO,'Main',{});
RO=struct;
if ~isempty(evt);[CAM,Cam]=comstr(evt.CAM,1);
   if isfield(evt,'RO');RO=evt.RO;end
else;[CAM,Cam]=comstr(CAM,4);
    if carg<=nargin;RO=varargin{carg};carg=carg+1;end
end
if comstr(Cam,'mv')

 PARAM=sdtroot('paramvh');
 if ~isfield(PARAM,'Mv')||~isempty(strfind(Cam,'reset'))
  %edit(fullfile(fileparts(which('mevib16')),'mevib_en-us.csv'))
  r1=sdt_locale('defcsv','mevib_en-us.csv');
  PARAM.Mv=cinguj('objEditJ',r1.Mv,'SDT Root');
 end
 r1=fe_def('cleanentry',PARAM.Mv);
 
 if isfield(RO,'PModes') 
  %% #PModes compute square plate modes with added mass -2
  RO=feval(sdtroot('@sfield'),'AddMissing',RO,r1);
  
  model=mevib16('MeshPlate'); 
  if ischar(RO.PM);RO.PM=str2double(RO.PM); end
  if ~isfinite(r1.PM);error('Inconsistent mass');end
  if ischar(RO.PNode);RO.PNode=str2double(RO.PNode);end 
  if ~ismember(RO.PNode,model.Node(:,1));
      error('NodeId %i does not exist',RO.PNode);
  end
  mo1=feutil('addelt',model,'mass1',[RO.PNode RO.PM*[1 1 1]]);
  def=fe_eig(mo1);
  
  cf=comgui('guifeplot -project "SDT Root"',2);
  feplot(mo1,def);fecom('ch10');fecom('colordataEvalZ');
  if sdtdef('verm')<806;fecom('renderer zbuffer');end
  ii_plp('colormap wcenter thres .05',jet(50))
  fecom('view3');iimouse('resetview');
  %fecom('colormap',jet(7))
 elseif isfield(RO,'BModes') 
  %% #BModes compute beam modes with added mass -2
  RO=feval(sdtroot('@sfield'),'AddMissing',RO,r1);
  model=mevib16('MeshBeam',RO); 
  d1=fe_eig(model,[5 30 0]);
  
  cf=comgui('guifeplot -project "SDT Root"',2);
  cf.model=model; 
  cf.def=d1; fecom('colordataEvalZ');fecom('renderer zbuffer');
 
 
 elseif isfield(RO,'BDFRF') 
  %% #BDFRF compute square plate transfer without & with damper -2
  RO=feval(sdtroot('@sfield'),'AddMissing',RO,r1);
  model=mevib16('MeshBeam',RO); 
  try;RO.BFreq=eval(RO.BFreq);
  catch;
     RO.BFreq=linspace(5,300,300);
  end
  model=stack_set(model,'info','Freq',RO.BFreq(:));
  % Nominal with no damping
  mo1=model;mo1.Elt=feutil('removeelt eltname celas | eltname mass',model);
  d0=fe_eig(mo1,[5 20 1e3]); 
  nor2xf(d0,.001,model,'iiplot -acc "base" -po');iicom('submagpha');
  
  % Reduced basis with damper
  [SE,CE]=fe_case('assemble -SE -matdes 2 1 4',model);
  d1=fe_eig(SE,[5 20 1e3]);
  damp=feutilb('dtkt',d1.def(fe_c(d1.DOF,CE.DOF,'ind'),:),SE.K);
  damp=damp(:,3)./damp(:,2)+.001;
  nor2xf(d1,.001,model,'iiplot -acc "damper 0.1%" -po');
  nor2xf(d1,damp,model,'iiplot -acc "damper true" -po');
  iicom('iixonly',{'base','damper 0.1%','damper true'})
  
  
  iicom('submagpha');
 %elseif ~isempty(fieldnames(RO)); error('%s nothing done',comstr(RO,-30));
 elseif isfield(RO,'BrModes') 
  %% #BrModes : display brake result -2
   mevib16('loadBrake -cf2');
 elseif isfield(RO,'BrId') 
  %% #BrId : display brake result -2
  mevib16('loadBrake -ci3');
 elseif isfield(RO,'BrFilt') ; mevib16('ViewBrFilt');
 elseif isstruct(RO)&&length(fieldnames(RO))==1
  st=fieldnames(RO);st=st{1};
  switch st
  case {'BrTxt','BTxt'};mevib16('LoadTxt');
  case 'BrMAC';mevib16('viewmac');
  end
   
 end
end
out=[];

%% wd : path utilities -------------------------------------------------------
elseif comstr(Cam,'wd')
 fname=varargin{carg};carg=carg+1;
 out=which(fname); 
 if isempty(out);out=fullfile(fileparts(which('mevib16.m')),fname);end
 if isempty(out)||~exist(out,'file');
     out=fullfile(sdeb('wd/ENSAM/savoir_TP/Mevib16'),fname);
 end
 %load(fullfile(which D:\balmes\Dropbox\ENSAM\savoir_TP\VibAM\DbTest.mat FRF ID SE TEST

elseif strcmp(Cam,'cvs')
    out='$Revision: 1.13 $  $Date: 2016/11/17 13:51:57 $'; return; 

else; error('%s unknown',CAM);
end
