
pw0=pwd;
wd=fullfile(pw0,'s68mevib');
if ~exist(wd,'dir');wd=fileparts(pw0);end
cd(wd);sdtcheck('path');
addpath(pw0);
try;mevib16('init');mevib16('init');end
cd(pw0);
try;cd(fullfile(getenv('USERPROFILE'),'Desktop'));end
web(fullfile(pw0,'html','tpModal.html'))

 mevib16('loadBrake -cf2'); % Display FEM modes button
 model
 def
 