function analyzePSCBatch(param, dataFolder, saveFolder, pscFolder, expPSCFolder, expSWRFolder)

%% Handle input arguments
if (nargin < 6) expSWRFolder = []; end
if (nargin < 5) expPSCFolder = []; end
if (nargin < 4) pscFolder    = []; end
if (nargin < 3) saveFolder   = []; end
if (nargin < 2) dataFolder   = []; end
if (nargin < 1) param        = struct; end

% Handle case in which empty variable is supplied:
if isempty(param) param      = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum')              param.fileNum              = 2;   end
if ~isfield(param,'importPSCOption')      param.importPSCOption      = 1;   end
if ~isfield(param,'swrPSCOption')         param.swrPSCOption         = 1;   end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 0;   end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption   = 1;   end
if ~isfield(param,'swrWindow')            param.swrWindow            = 100; end
if ~isfield(param,'parsePSCOption')       param.parsePSCOption       = 1;   end
if ~isfield(param,'calcEvMatrixOption')   param.calcEvMatrixOption   = 1;   end
if ~isfield(param,'expPSCEvOption')       param.expPSCEvOption       = 1;   end
if ~isfield(param,'expSWREvOption')       param.expSWREvOption       = 1;   end
if ~isfield(param,'gammaOption')          param.gammaOption          = 1;   end
if ~isfield(param,'gammaLim1')            param.gammaLim1            = 20;  end
if ~isfield(param,'gammaLim2')            param.gammaLim2            = 50;  end
if ~isfield(param,'rOption')              param.rOption              = 1;   end
if ~isfield(param,'rLim1')                param.rLim1                = 120; end
if ~isfield(param,'rLim2')                param.rLim2                = 220; end
if ~isfield(param,'spectOption')          param.spectOption          = 1;   end
if ~isfield(param,'spectLim1')            param.spectLim1            = 1;   end
if ~isfield(param,'spectLim2')            param.spectLim2            = 600; end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption      = 0;   end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% If not supplied, prompt for folder to analyze
if isempty(dataFolder)
  if param.reAnalyzeOption
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed *.mat files of LFP + PSCs');
  else
    dataFolder = uigetdir(pwd, 'Select folder containing *.mat files of analyzed LFP + imported cell channel');
  end
  if (dataFolder == 0) return; end
end
[parentPath, ~, ~] = parsePath(dataFolder);

% Select folder to save analyzed matlab files
if isempty(saveFolder)
  saveFolder = uigetdir(parentPath, 'Select same or alternate folder to save *.mat files');
  if (saveFolder == 0) error('No save folder selected'); end
  [parentPath, ~, ~] = parsePath(saveFolder);
end

% Select folder of PSC events exported from pClamp
if isempty(pscFolder) && param.importPSCOption && ~param.reAnalyzeOption
  pscFolder = uigetdir(parentPath, 'Select folder of PSC event *.txt files exported from pClamp');
  if (pscFolder == 0) error('No PSC event folder selected'); end
  [parentPath, ~, ~] = parsePath(pscFolder);
end

% Select folder to export PSC event files (if selected)
if isempty(expPSCFolder) && param.expPSCEvOption
  expPSCFolder = uigetdir(parentPath, 'Select folder to export PSC event *.csv files');
  if (expPSCFolder == 0) 
    warning('No PSC event files to be exported - folder not selected'); 
  else
    [parentPath, ~, ~] = parsePath(expPSCFolder);
  end
end

% Select export folder of SWR events (if selected)
if isempty(expSWRFolder) && param.expSWREvOption
  expSWRFolder = uigetdir(parentPath, 'Select folder to export updated SWR event *.csv files');
  if (expSWRFolder == 0) warning('No updated SWR event files to be exported - folder not selected'); end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract input file names
cd (dataFolder);
dir_temp   = dir('*.mat'); % Find only *.mat files
names      = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
dataFiles  = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nDataFiles = length(dataFiles);

if param.importPSCOption && ~param.reAnalyzeOption
  cd (pscFolder);
  dir_temp   = dir('*.txt'); % Find only *.txt files
  names      = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  pscFiles   = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nPSCFiles  = length(pscFiles);
end

cd (curPath);

% Error handling - file number mismatch:
if param.importPSCOption && ~param.reAnalyzeOption
  if (nDataFiles ~= nPSCFiles)
    error('Unequal number of files in data and PSC folders - analysis will be mismatched');
  end
end

% Determine individual file names:
dataFile{nDataFiles}   = [];
saveFile{nDataFiles}   = [];
pscFile{nDataFiles}    = [];
expPSCFile{nDataFiles} = [];
expSWRFile{nDataFiles} = [];

for i = 1:nDataFiles
  % Determine individual import files/folders
  dataFile{i} = [dataFolder slash dataFiles{i}];
  [~, dataFileName, ~] = parsePath(dataFile{i});
  
  if param.importPSCOption && ~param.reAnalyzeOption
    pscFile{i}  = [pscFolder slash pscFiles{i}];
  end
  
  % Determine individual output *.mat file names
  saveFile{i} = [saveFolder slash dataFileName '.mat'];
  
  % Determine individual exported PSC *.txt file names (if selected)
  if ~isempty(expPSCFolder) && param.expPSCEvOption
    if (expPSCFolder ~= 0)
      expPSCFile{i} = [expPSCFolder slash dataFileName '_pscEvents.csv'];
    end
  end
  
  % Determine individual exported SWR event *.csv file names (if selected)
  if ~isempty(expSWRFolder) && param.expSWREvOption
    if (expSWRFolder ~= 0)
      expSWRFile{i} = [expSWRFolder slash dataFileName '_swrEvents.csv'];
    end
  end
  
end

parfor i = 1:nDataFiles
  data = load(dataFile{i});
  analyzePSCFile(data, [], param, saveFile{i}, pscFile{i}, expPSCFile{i}, expSWRFile{i});
end
fprintf('complete\n');
end
