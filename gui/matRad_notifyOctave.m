% Copyright (C) 2020 nikla
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.

% -*- texinfo -*- 
% @deftypefn {} {@var{retval} =} matRad_notifyOctave (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: nikla <nikla@LAPTOP-NIKLAS>
% Created: 2020-06-04

function matRad_notifyOctave(hObject,eventName,evt)
    global eventMap;
    if isempty(eventMap)
        return;
    end
        
    persistent warningprinted;
    if isempty(warningprinted)
        hObject.showWarning('You are using an experimental implementation of ''notify'' concept for Octave');
        warningprinted=2;
    end
    
    objHandles = {eventMap(:).src};
    
    %clean Handles
    validH = cellfun(@isa,objHandles,repmat({'handle'},size(objHandles))); %cellfun(@ishandle,objHandles);
    objHandles = objHandles(validH);
    eventMap = eventMap(validH);
    
    %Find the object    
    objIx = cellfun(@(h) isequal(hObject,h), objHandles);
    
    objEvents = eventMap(objIx);
    
    if isempty(objEvents)
        return;
    end
    
    allNames = {objEvents(:).event};
    eventIx = find(strcmp(eventName,allNames)); 
    
    for runIx = 1:numel(eventIx)
        runEventIx = eventIx(runIx);
        runEvent = objEvents(runEventIx);
        if nargin < 3
            runEvent.callback(hObject);
        else
            runEvent.callback(hObject,evt);
    end 
end
