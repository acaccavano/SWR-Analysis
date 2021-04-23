function S = combCaCorr(param, dataFolder, saveFile)
%% S = combCaCorr(param, dataFolder, saveFile)
%
%  Function to combine folder of *.mat files with SWR-SWR and Cell-Cell correlations produced from analyzeCaFile

% Handle optional arguments
if (nargin < 3); saveFile   = []; end
if (nargin < 2); dataFolder = []; end
if (nargin < 1); param      = struct; end

% Handle case in which empty variable is supplied:
if isempty(param); param    = struct; end

% Set default parameters if not specified
if ~isfield(param,'nBins');          param.nBins          = 1000;   end
if ~isfield(param,'cellTypeOption'); param.cellTypeOption = 1;      end
if ~isfield(param,'cellType1');      param.cellType1      = 'Deep'; end
if ~isfield(param,'cellType2');      param.cellType2      = 'Supe'; end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% If not supplied, prompt for folders to analyze
if isempty(dataFolder)
  dataFolder = uigetdir(pwd, 'Select folder containing analyzed Calcium recordings');
end
if (dataFolder == 0); return; end

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
    [~, saveFileName, ~] = parsePath(saveFile);
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
S.corrVectorSWR    = [];
S.corrVectorCell   = [];

if param.cellTypeOption
  S.cdfC1X{nFiles} = [];
  S.cdfC1F{nFiles} = [];
  S.cdfC2X{nFiles} = [];
  S.cdfC2F{nFiles} = [];
  S.corrVectorC1   = [];
  S.corrVectorC2   = [];
end

j = 1;
k = 1;
l = 1;
for i = 1:nFiles
  data = load(dataFile{i});
  
  if data.Ca.nChannels >= 5
    S.cdfSWRX{j} = data.SWR.Ca.cdfX;
    S.cdfSWRF{j} = data.SWR.Ca.cdfF;
    S.cdfCellX{j} = data.Ca.SWR.cdfX;
    S.cdfCellF{j} = data.Ca.SWR.cdfF;
    S.corrVectorSWR = vertcat(S.corrVectorSWR, data.SWR.Ca.corrVector);
    S.corrVectorCell = vertcat(S.corrVectorCell, data.Ca.SWR.corrVector);
    j = j + 1;
    
    if param.cellTypeOption
      indC1  = strcmp(data.Ca.cellType, param.cellType1);
      indC2  = strcmp(data.Ca.cellType, param.cellType1);
      if sum(indC1) >= 5
        S.cdfC1X{k} = data.Ca.SWR.cdfC1X;
        S.cdfC1F{k} = data.Ca.SWR.cdfC1F;
        S.corrVectorC1 = vertcat(S.corrVectorC1, data.Ca.SWR.corrVectorC1);
        k = k + 1;
      else
        S.cdfC1X(k) = [];
        S.cdfC1F(k) = [];
      end
      
      if sum(indC2) >= 5
        S.cdfC2X{l} = data.Ca.SWR.cdfC2X;
        S.cdfC2F{l} = data.Ca.SWR.cdfC2F;
        S.corrVectorC2 = vertcat(S.corrVectorC2, data.Ca.SWR.corrVectorC2);
        l = l + 1;
      else
        S.cdfC2X(l) = [];
        S.cdfC2F(l) = [];
      end
    end
  else
    S.cdfSWRX(j)  = [];
    S.cdfSWRF(j)  = [];
    S.cdfCellX(j) = [];
    S.cdfCellF(j) = [];
    if param.cellTypeOption
      S.cdfC1X(k) = [];
      S.cdfC1F(k) = [];
      S.cdfC2X(l) = [];
      S.cdfC2F(l) = [];
    end
  end
end
nFiles   = length(S.cdfSWRX);
nFilesC1 = length(S.cdfC1X);
nFilesC2 = length(S.cdfC2X);

[cdfSWRCombF, cdfSWRCombX]   = ecdf(S.corrVectorSWR);
[cdfCellCombF, cdfCellCombX] = ecdf(S.corrVectorCell);
if param.cellTypeOption
  [cdfC1CombF, cdfC1CombX] = ecdf(S.corrVectorC1);
  [cdfC2CombF, cdfC2CombX] = ecdf(S.corrVectorC2);
end

% Initialize standardized CDF arrays
S.bins = linspace(0, 1, param.nBins)';
S.cdfSWRComb  = zeros(param.nBins, 1);
S.cdfCellComb = zeros(param.nBins, 1);
S.cdfSWR  = zeros(param.nBins, nFiles);
S.cdfCell = zeros(param.nBins, nFiles);

if param.cellTypeOption
  S.cdfC1Comb = zeros(param.nBins, 1);
  S.cdfC1     = zeros(param.nBins, nFilesC1);
  S.cdfC2Comb = zeros(param.nBins, 1);
  S.cdfC2     = zeros(param.nBins, nFilesC2);
end

% Standardize combined CDFs:
S.cdfSWRComb  = standardizeCDF(cdfSWRCombX, cdfSWRCombF, S.bins);
S.cdfCellComb = standardizeCDF(cdfCellCombX, cdfCellCombF, S.bins);
if param.cellTypeOption
  S.cdfC1Comb = standardizeCDF(cdfC1CombX, cdfC1CombF, S.bins);
  S.cdfC2Comb = standardizeCDF(cdfC2CombX, cdfC2CombF, S.bins);
end

% Standardize individual file CDFs:
for i = 1:nFiles
  S.cdfSWR(:, i)  = standardizeCDF(S.cdfSWRX{i}, S.cdfSWRF{i}, S.bins);
  S.cdfCell(:, i) = standardizeCDF(S.cdfCellX{i}, S.cdfCellF{i}, S.bins);
end

if param.cellTypeOption
  for i = 1:nFilesC1
    S.cdfC1(:, i) = standardizeCDF(S.cdfC1X{i}, S.cdfC1F{i}, S.bins);
  end
  for i = 1:nFilesC2
    S.cdfC2(:, i) = standardizeCDF(S.cdfC2X{i}, S.cdfC2F{i}, S.bins);
  end
end

% Average individual file CDFs:
S.cdfSWRAve = mean(S.cdfSWR, 2);
S.cdfSWRSEM = std(S.cdfSWR, 0, 2) / sqrt(size(S.cdfSWR, 2));
S.cdfCellAve = mean(S.cdfCell, 2);
S.cdfCellSEM = std(S.cdfCell, 0, 2) / sqrt(size(S.cdfCell, 2));
if param.cellTypeOption
  S.cdfC1Ave = mean(S.cdfC1, 2);
  S.cdfC1SEM = std(S.cdfC1, 0, 2) / sqrt(size(S.cdfC1, 2));
  S.cdfC2Ave = mean(S.cdfC2, 2);
  S.cdfC2SEM = std(S.cdfC2, 0, 2) / sqrt(size(S.cdfC2, 2));
end

%% Save file
if all(saveFileName)
  fprintf(['saving file ' saveFileName '... ']);
  save(saveFile,'-struct','S');
  fprintf('done\n');
end

fprintf('complete\n');

end


function cdfOut = standardizeCDF(cdfXIn, cdfFIn, bins)
% function to interpolate CFD to standardized nBins

  % Trim duplicate X-values (would cause interpolation to crash)
  [cdfX, indX] = unique(cdfXIn, 'last');
  cdfF = cdfFIn(indX);
  
  % Cap CDF with X,Y = (0,0) and (1,1) if not present
  if cdfX(1) ~= 0
    cdfX = vertcat(0, cdfX);
    cdfF = vertcat(0, cdfF);
  end
  
  if cdfX(end) ~= 1
    cdfX = vertcat(cdfX, 1);
    cdfF = vertcat(cdfF, 1);
  end
  
  cdfOut = interp1(cdfX, cdfF, bins, 'previous');
  
end