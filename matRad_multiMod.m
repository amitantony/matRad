% might need a circular phantom
%load BOXPHANTOM.mat

 %%
% meta information for treatment plan
% pln as arr struct makes more sense  ~(O_o)~

% photon bit pln(1)
pln(1).bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).gantryAngles    = [0]; % [?]
pln(1).couchAngles     = [0]; % [?]
pln(1).numOfBeams      = numel(pln(1).gantryAngles);
pln(1).numOfVoxels     = prod(ct.cubeDim);
pln(1).isoCenter       = ones(pln(1).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(1).voxelDimensions = ct.cubeDim;
pln(1).radiationMode   = 'photons';     % either photons / protons / carbon
%pln.radiationMode   = {'carbon','carbon','photons'};     % either photons / protons / carbon
pln(1).bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln(1).numOfFractions  = 1;
pln(1).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).machine         = 'Generic';

% proton pln(2)
pln(2).bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).gantryAngles    = [ 180]; % [?]
pln(2).couchAngles     = [ 0]; % [?]
pln(2).numOfBeams      = numel(pln(2).gantryAngles);
pln(2).numOfVoxels     = prod(ct.cubeDim);
pln(2).isoCenter       = ones(pln(2).numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln(2).voxelDimensions = ct.cubeDim;
pln(2).radiationMode   = 'protons';     % either photons / protons / carbon
%pln.radiationMode   = {'carbon','carbon','photons'};     % either photons / protons / carbon
pln(2).bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln(2).numOfFractions  = 1;
pln(2).runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).machine         = 'Generic';




%%

%stf = matRad_generateStf(ct,cst,pln); % do not foget this also need adjustment!!!
% for pln4 
   stf = matRad_generateStf(ct,cst,pln);

%% 
%todo write wrapper function
%dij = matRad_multiModDoseCalcWrapper(ct,stf,pln,cst);

dij = matRad_calcCombiDose(ct,stf,pln,cst);

%% fluence opti uses some parts of pln for options 
%  change basic pipeline for it 

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);    
%%
wInit = ones(sum([stf.totalNumOfBixels]),1);
wInit = resultGUI.w;

wInitxRay = wInit;
wInitProt = wInit;
wTot      = wInit;
off=1;
for i=1:numel(stf)
   if strcmp(stf(i).radiationMode, 'photons')
       off=off+stf(i).numOfRays;
   end
end

wInitxRay(off:end) = 0;
wInitProt(1:off-1)   = 0;

doseXRay   = reshape(dij.physicalDose{1} * wInitxRay,[ct.cubeDim]);
doseProton = reshape(dij.physicalDose{1} * wInitProt,[ct.cubeDim]);
doseTot    = reshape(dij.physicalDose{1} * wTot,[ct.cubeDim]);

slice = round(pln(1).isoCenter(1,3)/ct.resolution.z);

figure,
subplot(131),imagesc(doseXRay(:,:,slice)),colorbar
subplot(132),imagesc(doseProton(:,:,slice)),colorbar
subplot(133),imagesc(doseTot(:,:,slice)),colorbar
%%

   figure;
hold on,
 plot(resultGUI.physicalDose(:,80,80));
 plot(doseXRay(:,80,80));
 plot(doseProton(:,80,80));
 legend('total', 'Xray', 'Carbon');
 title ('Default constraints, only Physical dose');
 hold off;



%%
% doseIsoLevels = [0.2:0.1:1.2]*60;
% 
% figure,
% [hCMap,hDose,hCt,hContour,hIsoDose] =...
%     matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,3,slice,[],[],[],[],[],doseIsoLevels);
% 
% 
% 
% 





