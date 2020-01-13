function T = combPhaseStats(dataFolder, exportFile)

% Handle optional arguments
if (nargin < 2) exportFile = []; end
if (nargin < 1) dataFolder = []; end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

parentPath = [];

% If not supplied, prompt for data folder to analyze
if isempty(dataFolder)
  dataFolder = uigetdir(pwd, 'Select folder containing analyzed Spike recordings');
end
if (dataFolder == 0) return; end

% Parse dataFolder to determine default save name
[parentPath, dataFolderName, ~] = parsePath(dataFolder);

% Select folders to save analyzed matlab files
if isempty(exportFile)
  defaultPath = [parentPath 'phStats_' dataFolderName '.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of phase stats', defaultPath);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('No file to be exported - no file selected'); end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Initialize output data structure:
T = struct;

%% Import phase data
% Extract file names
cd (dataFolder);
dir_temp = dir('*.mat'); % Find only matlab files
names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nDataFiles = length(file);
cd (curPath);

% Determine individual file names:
dataFile{nDataFiles}    = [];

for i = 1:nDataFiles
  dataFile{i}  = [dataFolder slash file{i}];
end

fileName{nDataFiles, 1} = [];
gPhaseAve = NaN * ones(nDataFiles, 1);
gPhaseR   = NaN * ones(nDataFiles, 1);
gPhaseP   = NaN * ones(nDataFiles, 1);
gPhaseZ   = NaN * ones(nDataFiles, 1);
rPhaseAve = NaN * ones(nDataFiles, 1);
rPhaseR   = NaN * ones(nDataFiles, 1);
rPhaseP   = NaN * ones(nDataFiles, 1);
rPhaseZ   = NaN * ones(nDataFiles, 1);

for i = 1:nDataFiles
  data = load(dataFile{i});
  fileName{i} = data.saveName;
  
  if isfield(data.C.spike.gamma, 'phaseAve') 
    gPhaseAve(i) = data.C.spike.gamma.phaseAve;
    if (gPhaseAve(i) < 0) gPhaseAve(i) = gPhaseAve(i) + 2*pi; end
    gPhaseR(i) = data.C.spike.gamma.phaseR;
  end
  
  if isfield(data.C.spike.gamma, 'phaseP')
    gPhaseP(i) = data.C.spike.gamma.phaseP;
    gPhaseZ(i) = data.C.spike.gamma.phaseZ;
  end
  
  if isfield(data.C.spike.R, 'phaseAve')
    rPhaseAve(i) = data.C.spike.R.phaseAve;
    if (rPhaseAve(i) < 0) rPhaseAve(i) = rPhaseAve(i) + 2*pi; end
    rPhaseR(i) = data.C.spike.R.phaseR;
  end
  
  if isfield(data.C.spike.R, 'phaseP')
    rPhaseP(i) = data.C.spike.R.phaseP;
    rPhaseZ(i) = data.C.spike.R.phaseZ;
  end
  
end

varNames = {'gPhaseAve', 'gPhaseR', 'gPhaseP', 'gPhaseZ', 'rPhaseAve', 'rPhaseR', 'rPhaseP', 'rPhaseZ'};
T = table(gPhaseAve, gPhaseR, gPhaseP, gPhaseZ, rPhaseAve, rPhaseR, rPhaseP, rPhaseZ, 'VariableNames', varNames);

% Replace NaN values with blanks
tmp = table2cell(T);
tmp(isnan(T.Variables)) = {[]};
T = array2table(tmp,'VariableNames', T.Properties.VariableNames);

T = [table(fileName, 'VariableNames', {'fileName'}) T];
writetable(T, exportFile, 'Delimiter', ',');

fprintf('complete\n');

end
