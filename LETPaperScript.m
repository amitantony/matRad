%% My processing script 

% compress optimized results into one cube 
% this is a super shot doable script from the other branch 

if iscell(result_pre)
[resultGUI, result] = matRad_accumulateCubesMixMod(result_pre, pln,ct);
end

resultGUI.totalEffect = resultGUI.mod1effect + resultGUI.mod2effect;

%% DVH

%% ficures slice 

%% Visualization
slice = 59;

photon_plan = resultGUI{2};
proton_plan = resultGUI{1};
totalPlan = pln(1).numOfFractions.*proton_plan.(quantityOpt) + pln(2).numOfFractions.*photon_plan.(quantityOpt);

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


%% profiles ? 

