function S = combCaCorr(param, dataFolder, saveFile)

% Handle optional arguments
if (nargin < 3) saveFile   = []; end
if (nargin < 2) dataFolder = []; end
if (nargin < 1) param      = struct; end

% Handle case in which empty variable is supplied:
if isempty(param) param    = struct; end

% Set default parameters if not specified
if ~isfield(param,'nBins') param.nBins = 1000; end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

parentPath = [];

% If not supplied, prompt for folders to analyze
if isempty(dataFolder)
  dataFolder = uigetdir(pwd, 'Select folder containing analyzed Calcium recordings');
end
if (dataFolder == 0) return; end

% Parse dataFolder to determine default save name
[parentPath, dataFolderName, ~] = parsePath(dataFolder);

% Select folders to save analyzed matlab files
if isempty(saveFile)
  defaultPath = [parentPath dataFolderName '.mat'];
  [saveName, savePath] = uiputfile('.mat','Select file to save output matlab file', defaultPath);
  saveFile = [savePath saveName];
  if ~all(saveFile)
    warning('No file to be saved - no file selected');
  else
    [parentPath, saveFileName, ~] = parsePath(saveFile);
  end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract file names
cd (dataFolder);
dir_temp = dir('*.mat'); % Find only matlab files
names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nFiles = length(file);
cd (curPath);

% Determine individual file names:
dataFile{nFiles}    = [];

for i = 1:nFiles
  dataFile{i}  = [dataFolder slash file{i}];
end

% Initialize output data structure:
S = struct;

S.cdfSWRX{nFiles}  = [];
S.cdfSWRF{nFiles}  = [];
S.cdfCellX{nFiles} = [];
S.cdfCellF{nFiles} = [];

j = 1;
for i = 1:nFiles
  data = load(dataFile{i});
  
  if data.Ca.nChannels >= 5
    S.cdfSWRX{j} = data.SWR.Ca.cdfX;
    S.cdfSWRF{j} = data.SWR.Ca.cdfF;
    S.cdfCellX{j} = data.Ca.SWR.cdfX;
    S.cdfCellF{j} = data.Ca.SWR.cdfF;
    j = j + 1;
  else
    S.cdfSWRX(j)  = [];
    S.cdfSWRF(j)  = [];
    S.cdfCellX(j) = [];
    S.cdfCellF(j) = [];
  end
end
nFiles = length(S.cdfSWRX);

% Initialize standardized CDF arrays
S.bins = linspace(0, 1, param.nBins)';
S.cdfSWR  = zeros(param.nBins, nFiles);
S.cdfCell = zeros(param.nBins, nFiles);

for i = 1:nFiles
  
  % Interpolate SWR-SWR correlation CFD
  % Trim duplicate X-values (would cause interpolation to crash)
  [cdfX, indX] = unique(S.cdfSWRX{i}, 'last');
  cdfF = S.cdfSWRF{i}(indX);
  
  % Cap CDF with X,Y = (0,0) and (1,1) if not present
  if cdfX(1) ~= 0
    cdfX = vertcat(0, cdfX);
    cdfF = vertcat(0, cdfF);
  end
  
  if cdfX(end) ~= 1
    cdfX = vertcat(cdfX, 1);
    cdfF = vertcat(cdfF, 1);
  end
  
  S.cdfSWR(:, i) = interp1(cdfX, cdfF, S.bins, 'previous');
  
  % Interpolate Cell-Cell correlation CFD
  % Trim duplicate X-values (will crash interpolation)
  [cdfX, indX] = unique(S.cdfCellX{i}, 'last');
  cdfF = S.cdfCellF{i}(indX);
  
  % Cap CDF with X,Y = (0,0) and (1,1) if not present
  if cdfX(1) ~= 0
    cdfX = vertcat(0, cdfX);
    cdfF = vertcat(0, cdfF);
  end
  
  if cdfX(end) ~= 1
    cdfX = vertcat(cdfX, 1);
    cdfF = vertcat(cdfF, 1);
  end
  
  S.cdfCell(:, i) = interp1(cdfX, cdfF, S.bins, 'previous');
  
end

S.cdfSWRAve = mean(S.cdfSWR, 2);
S.cdfSWRSEM = std(S.cdfSWR, 0, 2) / sqrt(size(S.cdfSWR, 2));

S.cdfCellAve = mean(S.cdfCell, 2);
S.cdfCellSEM = std(S.cdfCell, 0, 2) / sqrt(size(S.cdfCell, 2));

%% Save file
if all(saveFileName)
  fprintf(['saving file ' saveFileName '... ']);
  save(saveFile,'-struct','S');
  fprintf('done\n');
end

fprintf('complete\n');

end
