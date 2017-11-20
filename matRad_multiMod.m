load BOXPHANTOM.mat

 
% meta information for treatment plan
% pln as arr struct makes more sense  ~(O_o)~

% photon bit pln(1)
pln(1).bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).gantryAngles    = [0:72:359]; % [?]
pln(1).couchAngles     = [0 0 0 0 0]; % [?]
pln(1).numOfBeams      = numel(pln(1).gantryAngles);
pln(1).numOfVoxels     = prod(ct.cubeDim);
pln(1).isoCenter       = ones(pln(1).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(1).voxelDimensions = ct.cubeDim;
pln(1).radiationMode   = 'photons';     % either photons / protons / carbon
%pln.radiationMode   = {'carbon','carbon','photons'};     % either photons / protons / carbon
pln(1).bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln(1).numOfFractions  = 30;
pln(1).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).machine         = 'Generic';
% proton pln(2)
pln(2).bixelWidth      = 7; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).gantryAngles    = [0]; % [?]
pln(2).couchAngles     = [0]; % [?]
pln(2).numOfBeams      = numel(pln(2).gantryAngles);
pln(2).numOfVoxels     = prod(ct.cubeDim);
pln(2).isoCenter       = ones(pln(2).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(2).voxelDimensions = ct.cubeDim;
pln(2).radiationMode   = 'protons';     % either photons / protons / carbon
%pln.radiationMode   = {'carbon','carbon','photons'};     % either photons / protons / carbon
pln(2).bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln(2).numOfFractions  = 5;
pln(2).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).machine         = 'Generic';




%%

%stf = matRad_generateStf(ct,cst,pln); % do not foget this also need adjustment!!!
% for pln2 
   stf = matRad_generateStf(ct,cst,pln);

%% 
%todo write wrapper function
%dij = matRad_multiModDoseCalcWrapper(ct,stf,pln,cst);

dij = matRad_calcCombiDose(ct,stf,pln,cst);

%% fluence opti uses some parts of pln for options 
%  change basic pipeline for it 

resultGUI = matRad_fluenceOptimization(dij,cst,pln);






