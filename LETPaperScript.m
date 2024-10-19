%% My processing script 
matRad_rc;
% clear
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 500000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-06;
% load TG119.mat
 %%
% %% add core
% cube = zeros(ct.cubeDim);
% cube(cst{1,4}{1}) = 1;
% vResolution = ct.resolution;
% vMargin = [];
% vMargin.x = 5;
% vMargin.y = 5;
% vMargin.z = 5;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin,1);
% 
% cst{4,1}    = 3;
% cst{4,2}    = 'Core_Big';
% cst{4,3}    = 'OAR';
% cst{4,4}{1} = find(mVOIEnlarged);
% cst{4,5}    = cst{1,5};
% 
% %% add Target margin
% cube = zeros(ct.cubeDim);
% cube(cst{2,4}{1}) = 1;
% vResolution = ct.resolution;
% vMargin = [];
% vMargin.x = 5;
% vMargin.y = 5;
% vMargin.z = 5;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin,1);
% 
% cst{5,1}    = 3;
% cst{5,2}    = 'Target_Conform';
% cst{5,3}    = 'OAR';
% cst{5,4}{1} = find(mVOIEnlarged);
% cst{5,5}    = cst{1,5};
% 
% cst{2,5}.alphaX  = 0.5;
% cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
%
cst{3,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30)); 
% cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,40)); 
% cst{3,6}{2} = struct(DoseObjectives.matRad_MeanDose(100,0));
%%

% meta information for treatment plan (1) 
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [250 290]; % [?] ;
%pln(1).propStf.gantryAngles    = [90];
pln(1).propStf.couchAngles     = zeros(numel(pln(1).propStf.gantryAngles),1); % [?] ; 
pln(1).propStf.numOfBeams      = numel(pln(1).propStf.gantryAngles);
pln(1).propStf.isoCenter       = ones(pln(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(1).propDoseCalc.calcLET = 1;

pln(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).propOpt.spatioTemp      = 0;
pln(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'MCN';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);

%pln = pln(1);

% meta information for treatment plan (2) 
pln(2).numOfFractions  = 25;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:40:359]; % [?] ;
pln(2).propStf.couchAngles     = zeros(numel(pln(2).propStf.gantryAngles),1);  % [?] ; 
pln(2).propStf.numOfBeams      = numel(pln(2).propStf.gantryAngles);
pln(2).propStf.isoCenter       = ones(pln(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).propOpt.spatioTemp      = 0;
pln(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = ['effect'];     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(2).bioParam = matRad_bioModel(pln(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);

% Plan Wrapper
plnJO = matRad_plnWrapper(pln);
% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);

%% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij,pln);
dij = matRad_calcmLETDose(dij,pln);

%% compress optimized results into one cube 
%
dij.precon = 1;
% dij.wInit = [result_pre{1}.w; result_pre{2}.w];
[result_pre,optimizer_U300_preCon] = matRad_fluenceOptimizationJO(dij,cst,plnJO);
% Visualization
slice = 65;

photon_plan = result_pre{2};
proton_plan = result_pre{1};
totalPlan = pln(1).numOfFractions.*proton_plan.(quantityOpt) + pln(2).numOfFractions.*photon_plan.(quantityOpt);

f = figure;
subplot(1,3,1);
    imagesc(proton_plan.(quantityOpt)(:,:,slice).*pln(1).numOfFractions);
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Proton Plan');
subplot(1,3,2);
    imagesc(photon_plan.(quantityOpt)(:,:,slice).*pln(1).numOfFractions);
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Photon Plan');
subplot(1,3,3);
    imagesc(totalPlan(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Total Plan');
%%
if iscell(result_pre)
[resultGUI, result] = matRad_accumulateCubesMixMod(result_pre, pln,ct);
end
% pln = pln(1) ;
% resultGUI.totalEffect = resultGUI.mod1effect + resultGUI.mod2effect;

%% DVH

%% plotting
color = ['b', 'm', 'g' , 'r'];
figure
subplot(3,1,1)
structure = [2 3 4];
for i = [1 2 3]
    plot(dvh_Pwithout_p(structure(i)).doseGrid,dvh_Pwithout_p(structure(i)).volumePoints,'Color',color(i),'LineStyle','-','LineWidth',2)
    hold on
end

%% Visualization
slice = 65;
ResultCell = result_DD_Full;
photon_plan = ResultCell{2};
proton_plan = ResultCell{1};
quantityOpt = 'effect';
totalPlan = pln(1).numOfFractions.*proton_plan.(quantityOpt) + pln(2).numOfFractions.*photon_plan.(quantityOpt);
% matRad_calcQualityIndicators(cst,pln,totalPlan)
%%



f = figure;
subplot(1,3,1);
    imagesc(proton_plan.(quantityOpt)(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Proton Plan');
subplot(1,3,2);
    imagesc(photon_plan.(quantityOpt)(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Photon Plan');
subplot(1,3,3);
    imagesc(totalPlan(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,slice);
    title('Total Plan');


%% %% ficures slice 

plane = 3;
slice = 65;
cube = proton_plan.(quantityOpt);
doseWindow = [0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
figure,
subplot(1,3,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title(['Referenz Proton ' quantityOpt])
zoom(1.5)
subplot(1,3,2)
cube = photon_plan.(quantityOpt);
doseWindow = [0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,photon_plan.(quantityOpt),plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title(['Referenz Photon ' quantityOpt])
zoom(1.5)
subplot(1,3,3)
cube = totalPlan;
doseWindow = [0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Referenz Proton BED')
zoom(1.5)
figure
cube = proton_plan.dirtyDose;
doseWindow = [0 max(cube(:))];
isoStep = [0:0.1*doseWindow(2):doseWindow(2)];
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,isoStep);
title('Referenz Proton BED')
zoom(1.5)



%% calc DVH for the current cube 
resultDVH  = matRad_calcDVH(cst,totalPlan,'cum');
% LET_DVH = resultDVH;
figure,
vois = [3,4,5,12,13,14,15];
for i  = vois
    plot(resultDVH(i).doseGrid,resultDVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-');
    hold on
    plot(LET_DVH(i).doseGrid,LET_DVH(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
end 
hold on 
c =1;

%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end

legend ([leg{:}],names ) 




