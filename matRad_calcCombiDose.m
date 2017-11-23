function dij = matRad_calCombiDose(ct,stf,pln,cst)

mod = numel(pln);
pntr=1;
 % over each beam

dij.numOfBeams =0;
dij.numOfVoxels = 0;
dij.resolution = 0;
dij.dimensions = 0;
dij.numOfScenarios = 0;   
dij.numOfRaysPerBeam = [];
dij.totalNumOfBixels = 0;
dij.totalNumOfRays = 0;
dij.bixelNum = [];
dij.rayNum = [];
dij.beamNum = [];
dij.physicalDose = {[]};
 
 
for i = 1:mod
    
% pick the correct pln for stf pls    
% for j = 1 : numel(pln)
% if strcmp( stf(i).radiationMode, pln(j).radiationMode)
% currPlan=pln(j);
% end
% end


currStf=stf(pntr : pntr + pln(i).numOfBeams -1);
    
if strcmp(pln(i).radiationMode,'photons')
    dijt = matRad_calcPhotonDose(ct,currStf,pln(i),cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln(i).radiationMode,'protons') || strcmp(pln(i).radiationMode,'carbon')
    dijt = matRad_calcParticleDose(ct,currStf,pln(i),cst);
    % precondition correction for the Dij for particle doses ( right now
    % arbitrary) 
    dijt.physicalDose{1} = dijt.physicalDose{1} .* 10;
end

% put the two dijs together ?? just like that ? 
dij.numOfBeams = dij.numOfBeams + dijt.numOfBeams;
dij.numOfVoxels = dijt.numOfVoxels;
dij.resolution = dijt.resolution;
dij.dimensions = dijt.dimensions;
dij.numOfScenarios = dijt.numOfScenarios;   % dunno wat to do 
dij.numOfRaysPerBeam = [dij.numOfRaysPerBeam dijt.numOfRaysPerBeam];
dij.totalNumOfBixels = dij.totalNumOfBixels + dijt.totalNumOfBixels;
dij.totalNumOfRays = dij.totalNumOfRays + dijt.totalNumOfRays;
dij.bixelNum = [dij.bixelNum; dijt.bixelNum];
dij.rayNum = [dij.rayNum; dijt.rayNum];
dij.beamNum = [dij.beamNum; dijt.beamNum];

dij.physicalDose = {[dij.physicalDose{1} dijt.physicalDose{1}]};

pntr = pntr+pln(i).numOfBeams;

end

fprintf('SUCCESS !! \n');



end

