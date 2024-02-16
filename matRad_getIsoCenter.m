function isoCenter = matRad_getIsoCenter(cst,ct,visBool)
% computes the isocenter [mm] as the joint center of gravity 
% of all volumes of interest that are labeled as target within the cst 
% struct
% 
% call
%   isoCenter = matRad_getIsoCenter(cst,ct,visBool)
%
% input
%   cst:        matRad cst struct
%   ct:         ct cube
%   visBool:    toggle on/off visualization (optional)
%
% output
%   isoCenter:  isocenter in [mm]   
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if visBool not set toogle off visualization
if nargin < 3
    visBool = 0;
end

% Initializes V variable.
V = [];

%Check if any constraints/Objectives have been defined yet
noObjOrConst = all(cellfun(@isempty,cst(:,6)));

% Save target indices in V variable.
for i = 1:size(cst,1)
    % We only let a target contribute if it has an objective/constraint or
    % if we do not have specified objectives/constraints at all so far
    if isequal(cst{i,3},'TARGET') && (~isempty(cst{i,6}) || noObjOrConst)
        V = [V; cst{i,4}{1}];
    end
end

% Delete repeated indices, one voxel can belong to two VOIs, because
% VOIs can be overlapping.
V = unique(V);

% throw error message if no target is found
if isempty(V)
    error('Could not find target');
end

% Transform subcripts from linear indices 
[coordV(:,1), coordV(:,2), coordV(:,3)] = ind2sub(ct.cubeDim,V);

% Transform to [mm]
coord = matRad_cubeToWorldCoordinates(coordV, ct);

% Calculated isocenter.
isoCenter = mean(coord);

% Visualization
if visBool

    clf
    hold on
    
    xCoordsV = coord (:,1); %xCoordsV * ct.resolution.x;
    yCoordsV = coord (:,2); %yCoordsV * ct.resolution.y;
    zCoordsV = coord (:,3); %zCoordsV * ct.resolution.z;

    % Plot target
    plot3(yCoordsV,xCoordsV,zCoordsV,'kx')
    
    % Show isocenter: red point
    plot3(isoCenter(2),isoCenter(1),isoCenter(3),'r.','MarkerSize',30)
    
    xlabel('y [mm]')
    ylabel('x [mm]')
    zlabel('z [mm]')
    
end
