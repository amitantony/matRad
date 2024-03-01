% test geometry

%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
%%
load("BOXPHANTOM.mat")

%% Create a CT image series

ctDim = [160,160,130]; % x,y,z dimensions
ctResolution = [3,3,3]; % x,y,z the same here!

% This uses the phantombuilder class, which helps to easily implement a 
% water phantom containing geometrical 3D shapes as targets and organs
builder = matRad_PhantomBuilder(ctDim,ctResolution,1);

%% Create the VOI data for the phantom
% To add a VOI there are (so far) two different options
% either a Box or a spherical Volume (either OAR or Target)
% NOTE: The order in which the objectives are initialized matters!
% In case of overlaps in the objectives, the firstly created objectives have
% a higher priority! This means that if two VOI have an overlap with
% different HU, then the value of the firstly initialized objective will be
% set in the overlap region


%define objectives for the VOI

objective1 = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(100,30));

builder.addSphericalTarget('PTV',10,'offset',[-5 -5 0],'objectives',objective1,'HU',0);
builder.addSphericalOAR('OAR1',12,'offset',[19 -7 0],'objectives',objective2);
builder.addSphericalOAR('OAR2',7,'offset',[-5 7 0],'objectives',objective2);
builder.addSphericalOAR('OAR3',7,'offset',[-5 -16 0],'objectives',objective2);
builder.addBoxOAR('Body',[60,70,60],'objectives',objective2);

% This adds two Box OAR and one Spherical Target in the middle
% For Box Volumes a name (here Volume2 and Volume3) has to be specified,
% as well as the dimension of the Volume.
% For spherical Volumes a name has to be specified, as well as the radius
% of the sphere
%(Optional) To move the VOI from the center of the ct an offset can be set.
%(Optional) The objectives can already be set here, however this can also
%be done later on
%(Optional) The HU of the VOI can be set (normal value: 0 (water)) 
%% Add Boxphantom
objective1 = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,30));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));

builder.addBoxOAR('OAR2',[10,20,80],'offset',[-12,0,0],'objectives',objective2);
builder.addBoxOAR('OAR3',[20,10,80],'offset',[0,-16,0],'objectives',objective2);
% builder.addBoxOAR('OAR3',[40,-40,40])

%% Get the ct and cst (stored as properties of the phantomBuilder class)

[ct2,cst2] = builder.getctcst();

%% Add to Boxphantom
cst{3,1}    = 2;
cst{3,2}    = cst2{1,2};
cst{3,3}    = 'OAR';
cst{3,4}{1} = cst2{1,4}{1};
cst{3,5}    = cst2{1,5};

cst{3,6}{1} = cst2{1,6}{1}; 

cst{4,1}    = 3;
cst{4,2}    = cst2{2,2};
cst{4,3}    = 'OAR';
cst{4,4}{1} = cst2{2,4}{1};
cst{4,5}    = cst2{2,5};

cst{4,6}{1} = cst2{2,6}{1}; 

%% Treatment Plan

pln.radiationMode = 'protons';            
pln.machine       = 'Generic';

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
% 'none':     physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model 
% and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
modelName    = 'none';
quantityOpt  = 'RBExD';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0:70:355];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
%%
matRadGUI;

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Export dij matrix
%matRad_exportDij('dij.bin',dij,stf);
%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%%
matRadGUI
%% Plot the resulting dose slice
plane      = 3;
slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([resultGUI.physicalDose(:)])];

figure,title('phantom plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);


%% 
% We export the the created phantom & dose as dicom. This is handled by the 
% class matRad_DicomExporter. When no arguments are given, the exporter searches
% the workspace itself for matRad-structures. The output directory can be set by
% the property dicomDir. While the different DICOM datasets (ct, RTStruct, etc) 
% can be exported individually, we call the wrapper to do all possible exports.
dcmExport = matRad_DicomExporter();
dcmExport.dicomDir = [pwd filesep 'dicomExport'];
dcmExport.matRad_exportDicom();

