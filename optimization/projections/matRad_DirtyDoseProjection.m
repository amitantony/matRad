classdef matRad_DirtyDoseProjection < matRad_BackProjection
% matRad_DirtyDoseProjection class to compute dirty dose during optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    methods
        function obj = matRad_DirtyDoseProjection()
            
        end
    end
    
    methods 
        function d = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.dirtyDose{scen})
                d = dij.dirtyDose{scen}*w;
            else
                d = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [dExp,dOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.dirtyDoseExp{scen})
                dExp = dij.dirtyDoseExp{scen}*w;

                for i = 1:size(dij.dirtyDoseOmega,2)
                   dOmegaV{scen,i} = dij.dirtyDoseOmega{scen,i} * w;
                end 
            else
                dExp = [];
                dOmegaV = [];
            end             
        end
     
        function wGrad = projectSingleScenarioGradient(~,dij,dirtyDoseGrad,scen,~)
            if ~isempty(dij.dirtyDose{scen})
               wGrad = (dirtyDoseGrad{scen}' * dij.dirtyDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dExpGrad,dOmegaVgrad,scen,~)
            if ~isempty(dij.dirtyDoseExp{scen})
                wGrad = (dExpGrad{scen}' * dij.dirtyDoseExp{scen})';
                wGrad = wGrad + 2 * dOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end