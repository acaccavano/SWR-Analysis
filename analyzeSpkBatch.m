function analyzeSpkBatch(param, dataFolder, saveFolder, spkFolder, bstFolder, expSpkFolder, expBstFolder, expSWRFolder, expAveFolder)
%% analyzeSpkBatch(param, dataFolder, saveFolder, spkFolder, bstFolder, expSpkFolder, expBstFolder, expSWRFolder)
%
%  Function to run batch analyzeSpkFile on batch of files
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   param      = structure containing all parameters including:
%     param.fileNum              = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.importSpkOption      = boolean flag to import spike file (needed unless reanalyzing) (default = 1)
%     param.swrSpkOption         = boolean flag to calculate coincidence of SWRs and spikes (default = 1)
%     param.swrBstOption         = boolean flag to calculate coincidence of SWRs and bursts (default = 1)
%     param.lfpSpkOption         = boolean flag to calculate spike-phase over duration of recording (instead of swrSpkOption - only 1 selectable, default = 0)
%     param.sdMultPhase          = for spike-phase coupling, the SD multiple the peak-trough amplitude must exceed in the oscillation surrounding the spike to be be considered
%     param.useSWRDurationOption = boolean flag to use detected SWR detection for coincidence detection (default = 1)
%     param.useSWRWindowOption   = boolean flag to use standard swrWindow for coincidence detection (default = 0)
%     param.swrWindow            = +/- window around SWR peak events (default = 100 ms)
%     param.parseSpkOption       = currently mandatory boolean flag to parse spikes into standard window (default = 1)
%     param.calcEvMatrixOption   = currently mandatory boolean flag to calc standard event SWR-Spk event matrix (default = 1)
%     param.calcPhaseStats       = Calculates statistics, WARNING: requires circ_stat toolbox (default = 1)
%     param.expSpkEvOption       = boolean flag to determine whether to export csv table of Spk events (default = 1)
%     param.expBstEvOption       = boolean flag to determine whether to export csv table of Bst events (default = 1)
%     param.expSWREvOption       = boolean flag to determine whether to export csv table of SWR events (default = 1)
%     param.reAnalyzeOption      = option to re-analyze file (default = 0)
%     param.expAveOption         = boolean flag to determine whether to export csv table of average statistics
%     param.nBins                = For spike histogram, not currently selectable from UI (default = 100)
%   dataFolder    = full path to matlab folder to import (if not set, will prompt)
%   saveFolder    = full path to matlab folder to save (can be same, if not set, will prompt)
%   spkFolder     = full path to pClamp spike event folder to import (if not set, will prompt)
%   bstFolder     = full path to pClamp burst event folder to import (if not set, will prompt)
%   expSpkFolder  = full path to spike event csv folder to export (if not set, will prompt)
%   expBstFolder  = full path to burst event csv folder to export (if not set, will prompt)
%   expSWRFolder  = full path to SWR event csv folder to export (if not set, will prompt)
%   expAveFolder  = full path to folder of exported csv table of averages (if not set and expAveOption = 1, will prompt)

%% Handle input arguments
if (nargin < 9) expAveFolder = []; end
if (nargin < 8) expSWRFolder = []; end
if (nargin < 7) expBstFolder = []; end
if (nargin < 6) expSpkFolder = []; end
if (nargin < 5) bstFolder    = []; end
if (nargin < 4) spkFolder    = []; end
if (nargin < 3) saveFolder   = []; end
if (nargin < 2) dataFolder   = []; end
if (nargin < 1) param        = struct; end

% Handle case in which empty variable is supplied:
if isempty(param) param      = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum')              param.fileNum              = 2;   end
if ~isfield(param,'importSpkOption')      param.importSpkOption      = 1;   end
if ~isfield(param,'swrSpkOption')         param.swrSpkOption         = 1;   end
if ~isfield(param,'swrBstOption')         param.swrBstOption         = 1;   end
if ~isfield(param,'lfpSpkOption')         param.lfpSpkOption         = 0;   end
if ~isfield(param,'sdMultPhase')          param.sdMultPhase          = 4;   end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 1;   end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption   = 0;   end
if ~isfield(param,'swrWindow')            param.swrWindow            = 100; end
if ~isfield(param,'parseSpkOption')       param.parseSpkOption       = 1;   end
if ~isfield(param,'calcEvMatrixOption')   param.calcEvMatrixOption   = 1;   end
if ~isfield(param,'calcPhaseStats')       param.calcPhaseStats       = 1;   end
if ~isfield(param,'expSpkEvOption')       param.expSpkEvOption       = 1;   end
if ~isfield(param,'expBstEvOption')       param.expBstEvOption       = 1;   end
if ~isfield(param,'expSWREvOption')       param.expSWREvOption       = 1;   end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption      = 0;   end
if ~isfield(param,'expAveOption')         param.expAveOption         = 1;   end
if ~isfield(param,'nBins')                param.nBins                = 100; end % Not yet selectable in UI

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% If not supplied, prompt for folder to analyze
if isempty(dataFolder)
  if param.reAnalyzeOption
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed *.mat files of LFP + Spikes');
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

% Select folder of spike events exported from pClamp
if isempty(spkFolder) && param.importSpkOption && ~param.reAnalyzeOption
  spkFolder = uigetdir(parentPath, 'Select folder of spike event *.txt files exported from pClamp');
  if (spkFolder == 0) error('No spike event folder selected'); end
  [parentPath, ~, ~] = parsePath(spkFolder);
end

% Select folder of burst events exported from pClamp
if isempty(bstFolder) && param.importSpkOption && param.swrBstOption && ~param.reAnalyzeOption
  bstFolder = uigetdir(parentPath, 'Select folder of burst event *.txt files exported from pClamp');
  if (bstFolder == 0) error('No burst event folder selected'); end
  [parentPath, ~, ~] = parsePath(bstFolder);
end

% Select folder to export spike event files (if  selected)
if isempty(expSpkFolder) && param.expSpkEvOption
  expSpkFolder = uigetdir(parentPath, 'Select folder to export spike event *.csv files');
  if (expSpkFolder == 0) 
    warning('No spike event files to be exported - folder not selected'); 
  else
    [parentPath, ~, ~] = parsePath(expSpkFolder);
  end
end

% Select folder to export burst event files (if selected)
if isempty(expBstFolder) && param.expBstEvOption
  if (expSpkFolder == 0)
    expBstFolder = uigetdir(parentPath, 'Select folder to export burst event *.csv files');
    if (expBstFolder == 0) 
      warning('No burst event files to be exported - folder not selected'); 
    else
      [parentPath, ~, ~] = parsePath(expBstFolder);
    end
  else
    expBstFolder = expSpkFolder;
  end
end

% Select export folder of SWR events (if selected)
if isempty(expSWRFolder) && param.expSWREvOption
  expSWRFolder = uigetdir(parentPath, 'Select folder to export updated SWR event *.csv files');
  if (expSWRFolder == 0) 
    warning('No updated SWR event files to be exported - folder not selected'); 
  else
    [parentPath, ~, ~] = parsePath(expSWRFolder);
  end
end

% Select folder to export average files, if option selected
if isempty(expAveFolder) && param.expAveOption
  expAveFolder = uigetdir(parentPath, 'Select folder to export average statistics *.csv files');
  if (expAveFolder == 0) warning('No files to be exported - average folder not selected'); end
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

if param.importSpkOption && ~param.reAnalyzeOption
  cd (spkFolder);
  dir_temp   = dir('*.txt'); % Find only *.txt files
  names      = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  spkFiles   = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nSpkFiles  = length(spkFiles);
end

if param.importSpkOption && param.swrBstOption && ~param.reAnalyzeOption
  cd (bstFolder);
  dir_temp  = dir('*.txt'); % Find only *.txt files
  names     = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  bstFiles  = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nBstFiles = length(bstFiles);
end

cd (curPath);

% Error handling - file number mismatch:
if param.importSpkOption && ~param.reAnalyzeOption
  if (nDataFiles ~= nSpkFiles)
    error('Unequal number of files in data and spike folders - analysis will be mismatched');
  end
  
  if param.swrBstOption
    if (nDataFiles ~= nBstFiles)
      error('unequal number of files in data and burst folders - analysis will be mismatched');
    end
  end
end

% Determine individual file names:
dataFile{nDataFiles}   = [];
saveFile{nDataFiles}   = [];
spkFile{nDataFiles}    = [];
bstFile{nDataFiles}    = [];
expSpkFile{nDataFiles} = [];
expBstFile{nDataFiles} = [];
expSWRFile{nDataFiles} = [];
expAveFile{nDataFiles} = [];

for i = 1:nDataFiles
  
  % Determine individual import files/folders
  dataFile{i}  = [dataFolder slash dataFiles{i}];
  [~, dataFileName, ~] = parsePath(dataFile{i});
  
  if param.importSpkOption && ~param.reAnalyzeOption
    spkFile{i} = [spkFolder slash spkFiles{i}];
    if param.swrBstOption
      bstFile{i} = [bstFolder slash bstFiles{i}];
    end
  end

  % Determine individual output *.mat file names
  saveFile{i} = [saveFolder slash dataFileName '.mat'];
  
  % Determine individual exported spike *.txt file names (if selected)
  if ~isempty(expSpkFolder) && param.expSpkEvOption
    if (expSpkFolder ~= 0)
      expSpkFile{i} = [expSpkFolder slash dataFileName '_spkEvents.csv'];
    end
  end
  
  % Determine individual exported burst *.txt file names (if selected)
  if ~isempty(expBstFolder) && param.expBstEvOption
    if (expBstFolder ~= 0)
      expBstFile{i} = [expBstFolder slash dataFileName '_bstEvents.csv'];
    end
  end
  
  % Determine individual exported SWR event *.csv file names (if selected)
  if ~isempty(expSWRFolder) && param.expSWREvOption
    if (expSWRFolder ~= 0)
      expSWRFile{i} = [expSWRFolder slash dataFileName '_swrEvents.csv'];
    end
  end

  % Determine individual exported average *.csv file names (if selected)
  if ~isempty(expAveFolder) && param.expAveOption
    if (expAveFolder ~= 0)
      expAveFile{i} = [expAveFolder slash dataFileName '_aveStats.csv'];
    end
  end
  
end

parfor i = 1:nDataFiles
  data = load(dataFile{i});
  analyzeSpkFile(data, [], param, saveFile{i}, spkFile{i}, bstFile{i}, expSpkFile{i}, expBstFile{i}, expSWRFile{i}, expAveFile{i});
end
fprintf('complete\n');
end
