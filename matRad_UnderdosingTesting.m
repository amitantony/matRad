%% Proton Optimization 2 beam without dirty Dose
matRad_rc;

load("BOXPHANTOM.mat")

% cube = zeros(183,183,90);
% cube(cst{6,4}{1}) = 1;
% vResolution = ct.resolution;
% vMargin = [];
% vMargin.x = 5;
% vMargin.y = 5;
% vMargin.z = 5;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);
% 
% cst{11,1}    = 10;
% cst{11,2}    = 'Margin';
% cst{11,3}    = 'OAR';
% cst{11,4}{1} = find(mVOIEnlarged);
% cst{11,5}    = cst{3,5};
% 
% cst{11,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,60)); 
% 
% cube = zeros(183,183,90);
% cube(cst{7,4}{1}) = 1;
% vResolution = ct.resolution;
% vMargin = [];
% vMargin.x = 5;
% vMargin.y = 5;
% vMargin.z = 5;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);
% 
% cst{12,1}    = 11;
% cst{12,2}    = 'Margin';
% cst{12,3}    = 'OAR';
% cst{12,4}{1} = find(mVOIEnlarged);
% cst{12,5}    = cst{3,5};
% 
% cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,60)); 

%cst{8,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,0));

cst{1,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,30));


%% Define radiation modality
plnO5.radiationMode = 'protons';        
plnO5.machine       = 'Generic';

% Calculate LET
plnO5.propDoseCalc.calcLET = 1;

% Set some beam parameters
plnO5.numOfFractions        = 5;
plnO5.propStf.gantryAngles  = [90];
plnO5.propStf.couchAngles   = zeros(numel(plnO5.propStf.gantryAngles),1);
plnO5.propStf.bixelWidth    = 5;
plnO5.propStf.numOfBeams    = numel(plnO5.propStf.gantryAngles);
plnO5.propStf.isoCenter     = ones(plnO5.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
plnO5.propOpt.runDAO        = 0;
plnO5.propSeq.runSequencing = 0;


% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

plnO5.bioParam = matRad_bioModel(plnO5.radiationMode,quantityOpt, modelName);
plnO5.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
plnO5.propDoseCalc.doseGrid.resolution.x = 6; % [mm]
plnO5.propDoseCalc.doseGrid.resolution.y = 6; % [mm]
plnO5.propDoseCalc.doseGrid.resolution.z = 6; % [mm]

% Generate Beam Geometry STF
stfO5 = matRad_generateStf(ct,cst,plnO5);

% Dose Calculation
dijO5 = matRad_calcParticleDose(ct,stfO5,plnO5,cst);

%% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dijO5 = matRad_calcDirtyDose_O(2,dijO5);

% Inverse Optimization for IMPT
resultGUI_O5 = matRad_fluenceOptimization(dijO5,cst,plnO5);

% resultGUIRefNew = resultGUI;
%%
phys_P10_L4 = matRad_calcQualityIndicators(cst,pln,resultGUI_p10_L4.physicalDose);
RBExD_P10_L4 = matRad_calcQualityIndicators(cst,pln,resultGUI_p10_L4.RBExD);
LET_P10_L4 = matRad_calcQualityIndicators(cst,pln,resultGUI_p10_L4.LETd);

%%
cube = result_O_DD{1,1}.RBExD;
plane = 3;
slice = 80;
doseWindow = [0 max(cube(:))];

figure,
subplot(2,2,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('Old branch, RBExD')
zoom(1.3)

cube = result_O_mL;
plane = 3;
slice = 80;
doseWindow = [0 max(cube(:))];

figure,
subplot(2,2,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('P10 L4 plan')
zoom(1.3)

cube = resultGUI_p80_L5.RBExD;
plane = 3;
slice = 80;
doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('P80 L5 plan')
zoom(1.3)

cube = resultGUI_p70_L5.RBExD;
plane = 3;
slice = 80;
doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('P70 L5 plan')
zoom(1.3)

cube = resultGUI_p70_L4.RBExD - resultGUI_p10_L4.RBExD;
plane = 3;
slice = 80;
doseWindow = [-0.2 0.2];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan P70 minus P10')
zoom(1.3)

cube = resultGUI_p100_L5.RBExD - resultGUI_p80_L5.RBExD;
plane = 3;
slice = 80;
doseWindow = [-0.2 0.2];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan P100 minus P80')
zoom(1.3)

cube = resultGUI_p100_L5.RBExD - resultGUI_p70_L5.RBExD;
plane = 3;
slice = 80;
doseWindow = [-0.2 0.2];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan P100 minus P70')
zoom(1.3)


% maximumRef = max(resultGUIRef.RBExD(cst{9,4}{1}));
% minimumRef = min(resultGUIRef.RBExD(cst{9,4}{1}));

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{11,1}    = 10;
cst{11,2}    = 'Margin';
cst{11,3}    = 'OAR';
cst{11,4}{1} = find(mVOIEnlarged);
cst{11,5}    = cst{3,5};

cst{11,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30)); 

cube = zeros(183,183,90);
cube(cst{4,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{12,1} = 11;
cst{12,2} = 'Margin';
cst{12,3} = 'OAR';
cst{12,4}{1} = find(mVOIEnlarged);
cst{12,5} = cst{11,5};

cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cube = zeros(183,183,90);
cube(cst{10,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{13,1} = 12;
cst{13,2} = 'Margin';
cst{13,3} = 'OAR';
cst{13,4}{1} = find(mVOIEnlarged);
cst{13,5} = cst{11,5};

cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
cube(cst{1,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{14,1} = 13;
cst{14,2} = 'Margin';
cst{14,3} = 'OAR';
cst{14,4}{1} = find(mVOIEnlarged);
cst{14,5} = cst{11,5};

cst{14,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cst{9,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(300,maximumRef));
cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(10,10));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 240 270];
pln.propStf.couchAngles   = [0 0 0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultTarget = resultGUI;

cube = resultTarget.RBExD;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
zoom(1.3)

cube = resultTarget.dirtyDose;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
zoom(1.3)

%% Differences

cube = resultGUIRef.RBExD - resultTarget.RBExD;
plane = 3;
slice = 34;
doseWindow = [-1 1];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')
zoom(1.3)

cube = resultGUIRef.dirtyDose - resultTarget.dirtyDose;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')
zoom(1.3)
