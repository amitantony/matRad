function [nameList, classList, handleList] = getAvailableEnginesOctave(pln,optionalPath)
    % Helper function which does exactly the same as getAvailableEngines
    % only Octave compatible, could also be solved with a lot of if
    % conditions in the main function, but is cleaner and better to debug
    %
    % call:
    %   [nameList, handleList] = DoseEngines.matRad_DoseEngine.getAvailableEngines(pln,optional_path)
    %
    % input:
    %   pln: containing proposed dose calc and machine file informations
    %   optionalPath: searches for dose calc engines in given
    %
    % returns:
    %   nameList: cell-array conatining readable names for engines
    %   classList: cell-array conatining full classnamens for available engines
    %   handleList: cell-array containing function-handles to
    %                  available engines constructor (call the included handle by adding Parentheses e.g. handleList{1}())

    matRad_cfg = MatRad_Config.instance();

    nameList = {};
    classList = {};
    handleList = {};

    switch nargin

        case 0
            % meta.package works perfectly in matlab but won't work in Octave
            % because Octave doesn't fully support package folder and class definition,
            % so we use a somewhat unclean work around in Octave which propably won't work withg class folder.
            % So when using engine classes in octave it's better to initionlize them from hand
            % get all files with engine in name and ending with
            packageContent = dir([matRad_cfg.matRadRoot filesep 'doseCalc' filesep '+DoseEngines' filesep '*.m']);
            subfolderContent = dir([matRad_cfg.matRadRoot filesep 'doseCalc' filesep '+DoseEngines' filesep '*' filesep '*.m']);

            % concatenate all .m files and class folder
            files = vertcat(packageContent, subfolderContent);

            % itterate through the meta classes in the package
            for i = 1:length(files)

                [~,className] = fileparts(files(i).name);
                % get meta class with package name infront
                mc = meta.class.fromName(['DoseEngines.' className]);
                
                % Check if we found a class
                if isempty(mc)
                    continue;
                end

                % check for the isCalcEngine property,
                % could be done cleaner with the superclasses method
                % which sadly isn't available in octave
                propNames = cellfun(@(x) x.Name, mc.PropertyList, 'UniformOutput', false);
                mc_isEngineIdx = find(strcmp('isDoseEngine', propNames));
                if(mc_isEngineIdx && mc.PropertyList{mc_isEngineIdx}.DefaultValue)

                    handleList{end+1} = str2func(mc.Name);
                    classList{end+1} = mc.Name;

                    %get readable name from metaclass properties
                    nameProp = mc.PropertyList(strcmp(propNames, 'name'));
                    nameProp = nameProp{1};

                    if (nameProp.Abstract)
                        % ND -> not defined meaning abstract class
                        % without a name
                        nameList{end+1} = 'ND';
                    else
                        nameList{end+1} = nameProp.DefaultValue;
                    end


                end

            end
        case 1
            % meta.package works perfectly in matlab but won't work in Octave
            % because Octave doesn't fully support package folder and class definition,
            % so we use a somewhat unclean work around in Octave which propably won't work withg class folder.
            % So when using engine classes in octave it's better to initionlize them from hand
            % get all files with engine in name and ending with
            packageContent = dir([matRad_cfg.matRadRoot filesep 'doseCalc' filesep '+DoseEngines' filesep '*.m']);
            subfolderContent = dir([matRad_cfg.matRadRoot filesep 'doseCalc' filesep '+DoseEngines' filesep '*' filesep '*.m']);

            % concatenate all .m files and class folder
            files = vertcat(packageContent, subfolderContent);

            % itterate through the meta classes in the package
            for i = 1:length(files)

                [~,className] = fileparts(files(i).name);
                % get meta class with package name infront
                mc = meta.class.fromName(['DoseEngines.' className]);
                % skip class if abstract
                if ~isempty(mc) && ~(mc.Abstract)
                    % check for the isCalcEngine property,
                    % could be done cleaner with the superclasses method
                    % which sadly isn't available in octave
                    propNames = cellfun(@(x) x.Name, mc.PropertyList, 'UniformOutput', false);
                    mc_isEngineIdx = find(strcmp('isDoseEngine', propNames));

                    if (mc_isEngineIdx && mc.PropertyList{mc_isEngineIdx}.DefaultValue)

                        % get radiation mode from meta class property
                        radModeIdx = find(strcmp('possibleRadiationModes', propNames));
                        propValue = mc.PropertyList{radModeIdx}.DefaultValue;

                        if(any(strcmp(propValue, pln.radiationMode)))
                            % get radiation mode from the in pln proposed basedata machine file
                            machineMode = DoseEngines.matRad_DoseEngine.loadMachine(pln).meta.radiationMode;

                            % add current class to return lists if the
                            % radiation mode is compatible
                            if(any(strcmp(propValue, machineMode)))
                                handleList{end+1} = str2func(mc.Name);
                                classList{end+1} = mc.Name;

                                %get readable name from metaclass properties
                                nameProp = mc.PropertyList(strcmp(propNames, 'name'));
                                nameProp = nameProp{1};

                                if (nameProp.Abstract)
                                    % ND -> not defined meaning abstract class
                                    % without a name
                                    nameList{end+1} = 'ND';
                                else
                                    nameList{end+1} = nameProp.DefaultValue;
                                end

                            end

                        end

                    end

                end

            end

        case 2
            % check if path is valid and add it to the current
            % matlab path
            if(isfolder(optionalPath))
                addpath(optionalPath);
            end

            % in matlab what work pretty good but in octave it won't read class folders,
            % so we have to get the files ourself
            folderContent = dir([optionalPath '*.m']);
            subfolderContent = dir([optionalPath filesep '**' filesep '*.m']);

            % concatenate all .m files and class folder
            files = vertcat(folderContent, subfolderContent);

            for i = 1:length(files)
                [~,className] = fileparts(files(i).name);
                mc = meta.class.fromName(className);
                if (~isempty(mc) && ~(mc.Abstract))

                    % check for the isCalcEngine property,
                    % could be done cleaner with the superclasses method
                    % which sadly isn't available in octave
                    propNames = cellfun(@(x) x.Name, mc.PropertyList, 'UniformOutput', false);
                    mc_isEngineIdx = find(strcmp('isDoseEngine', propNames));

                    if (mc_isEngineIdx && mc.PropertyList{mc_isEngineIdx}.DefaultValue)

                        % get radiation mode from meta class property
                        radModeIdx = find(strcmp('possibleRadiationModes', propNames));
                        propValue = mc.PropertyList{radModeIdx}.DefaultValue;

                        if(any(strcmp(propValue, pln.radiationMode)))
                            % get radiation mode from the in pln proposed basedata machine file
                            machineMode = DoseEngines.matRad_DoseEngine.loadMachine(pln).meta.radiationMode;

                            % add current class to return lists if the
                            % radiation mode is compatible
                            if(any(strcmp(propValue, machineMode)))
                                handleList{end+1} = str2func(mc.Name);
                                classList{end+1} = mc.Name;

                                %get readable name from metaclass properties
                                nameProp = mc.PropertyList(strcmp(propNames, 'name'));
                                nameProp = nameProp{1};

                                if (nameProp.Abstract)
                                    % ND -> not defined meaning abstract class
                                    % without a name
                                    nameList{end+1} = 'ND';
                                else
                                    nameList{end+1} = nameProp.DefaultValue;
                                end

                            end
                        end
                    end

                end
            end
    end
end