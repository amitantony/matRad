classdef matRad_SquaredOverdosingDirtyDose < DirtyDoseObjectives.matRad_DirtyDoseObjective
% matRad_SquaredOverdosingDirtyDose implements a penalized dirty dose 
%   See matRad_DirtyDoseObjective for interface description
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        name = 'Squared Overdosing Dirty Dose';
        parameterNames = {'d^{max}'}; 
        parameterTypes = {'dirtyDose'}; 
    end
    
    properties
        parameters = {30}; 
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredOverdosingDirtyDose(penalty,dMax) 
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DirtyDoseObjectives.matRad_DirtyDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 2 && isscalar(dMax)
                    obj.parameters{1} = dMax;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
       %% Calculates the Objective Function value
        function fDirtyDose = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            % overdose : dirtyDose minus prefered dose
            overdose = dirtyDose - obj.parameters{1};
            
            % apply positive operator
            overdose(overdose<0) = 0;
            
            % calculate objective function
            fDirtyDose = 1/numel(dirtyDose) * (overdose'*overdose);
        end
        
        %% Calculates the Objective Function gradient
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            % overdose : dirtyDose minus prefered dose
            overdose = dirtyDose - obj.parameters{1};
            
            % apply positive operator
            overdose(overdose<0) = 0;
            
            % calculate delta
            fDirtyDoseGrad = 2 * 1/numel(dirtyDose) * overdose;
        end
    end
    
end
