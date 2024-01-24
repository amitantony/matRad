function dij = matRad_calcDirtyDose(LET_thres,dij,pln)
% Calculates Dirty and Clean Dose by using LET threshold
%
% call
%   dij = matRad_calcDirtyDose(LET_thres,dij)
%
% input
%   LET_thres:  LET threshold, above: dirty dose, below: clean dose
%   dij:        matRad dij struct
%
% output
%   dij:        matRad dij struct with dirty and clean dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('dij','var') || isempty(dij)
    disp('Not enough input arguments! Calculation is not working.')
else
    if ~exist('LET_thres','var') || isempty(LET_thres)
        disp('Not enough input arguments! Calculation is not working.')
    else
        dij.dirtyDoseThreshold                       = LET_thres;
        [dij.LETmaskDirty,dij.LETmaskClean,dij.mLET] = matRad_calcLETmask(dij);

        for k = 1:dij.numOfModalities
            dij.dirtyDose{1,k}                       = dij.LETmaskDirty{1,k} .* dij.original_Dijs{1,k}.physicalDose{1};
            dij.cleanDose{1,k}                       = dij.LETmaskClean{1,k} .* dij.original_Dijs{1,k}.physicalDose{1};
        end
        
        m = 2;
        rad = ["protons","photons","carbon","helium"];

        if pln(2).radiationMode == rad(1,2)
            DirtySize = zeros(size(dij.original_Dijs{1,m-1}.physicalDose{1}));
            LETDirty_Size = size(dij.dirtyDose{1,m});
            row = LETDirty_Size(1,m-1);
            column = LETDirty_Size(1,m);
            DirtySize(1:row,1:column) = dij.dirtyDose{1,m};
            photonDD = sparse(DirtySize);
            dij.dirtyDose{1,3} = dij.dirtyDose{1,m-1} + photonDD;
            DirtySize(1:row,1:column) = dij.cleanDose{1,m};
            photonCD = sparse(DirtySize);
            dij.cleanDose{1,3} = dij.cleanDose{1,m-1} + photonCD;
        elseif pln(2).radiationMode == rad(1,3)
            DirtySize = zeros(size(dij.original_Dijs{1,m}.physicalDose{1}));
            LETDirty_Size = size(dij.dirtyDose{1,m-1});
            row = LETDirty_Size(1,m-1);
            column = LETDirty_Size(1,m);
            DirtySize(1:row,1:column) = dij.dirtyDose{1,m-1};
            protonDD = sparse(DirtySize);
            dij.dirtyDose{1,3} = dij.dirtyDose{1,m} + protonDD;
            DirtySize(1:row,1:column) = dij.cleanDose{1,m-1};
            protonCD = sparse(DirtySize);
            dij.cleanDose{1,3} = dij.cleanDose{1,m} + protonCD;
        end

    end
end

end