function analyzeCaBatch(param, dataFolder, saveFolder, CaFolder, timingFile, expCaFolder, expSWRFolder, expStimFolder)

%% Handle input arguments
if (nargin < 8) expStimFolder = []; end
if (nargin < 7) expSWRFolder  = []; end
if (nargin < 6) expCaFolder   = []; end
if (nargin < 5) timingFile    = []; end
if (nargin < 4) CaFolder      = []; end
if (nargin < 3) saveFolder    = []; end
if (nargin < 2) dataFolder    = []; end
if (nargin < 1) param         = struct; end

% Handle case in which empty variable is supplied:
if isempty(param) param      = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum')              param.fileNum              = 2;   end
if ~isfield(param,'interpOption')         param.interpOption         = 1;   end
if ~isfield(param,'samplingInt')          param.samplingInt          = 0.5; end
if ~isfield(param,'baseCorrectMethod')    param.baseCorrectMethod    = 2;    end
if ~isfield(param,'CaFiltLim1')           param.CaFiltLim1           = 0.03; end
if ~isfield(param,'CaFiltLim2')           param.CaFiltLim2           = 4;    end
if ~isfield(param,'CaFiltOrder')          param.CaFiltOrder          = 80;   end
if ~isfield(param,'CaFiltAlpha')          param.CaFiltAlpha          = 2.5;  end
if ~isfield(param,'smoothFactor')         param.smoothFactor         = 0.25; end
if ~isfield(param,'peakDetectCa')         param.peakDetectCa         = 1;   end
if ~isfield(param,'baseDetectMethod')     param.baseDetectMethod     = 2;   end
if ~isfield(param,'baseQuant')            param.baseQuant            = 0.8; end
if ~isfield(param,'sdMult')               param.sdMult               = 4;   end
if ~isfield(param,'skipDetectLim')        param.skipDetectLim        = 2;   end % s
if ~isfield(param,'consThreshOption')     param.consThreshOption     = 0;   end
if ~isfield(param,'swrCaOption')          param.swrCaOption          = 0;   end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 1;   end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption   = 0;   end
if ~isfield(param,'swrWindow')            param.swrWindow            = 100; end
if ~isfield(param,'expCaEvOption')        param.expCaEvOption        = 1;   end
if ~isfield(param,'expSWREvOption')       param.expSWREvOption       = 0;   end
if ~isfield(param,'spkCaOption')          param.spkCaOption          = 0;   end
if ~isfield(param,'stimCaOption')         param.stimCaOption         = 0;   end
if ~isfield(param,'stimCaLim1')           param.stimCaLim1           = 0;   end  % ms
if ~isfield(param,'stimCaLim2')           param.stimCaLim2           = 1000; end % ms
if ~isfield(param,'expStimEvOption')      param.expStimEvOption      = 0;   end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption      = 0;   end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

parentPath = [];

% Prompt first for matlab data folder of analyzed files, if options require:
if isempty(dataFolder)
  if param.reAnalyzeOption
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed dFoF *.mat files');
    if (dataFolder == 0) return; end
    [parentPath, ~, ~] = parsePath(dataFolder);
  elseif (param.swrCaOption || param.spkCaOption || param.stimCaOption)
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed LFP and/or cell *.mat files');
    if (dataFolder == 0) return; end
    [parentPath, ~, ~] = parsePath(dataFolder);
  else
    parentPath = pwd;
  end
else
  [parentPath, ~, ~] = parsePath(dataFolder);
end

% Select folders to save analyzed matlab files
if isempty(saveFolder)
  saveFolder = uigetdir(parentPath, 'Select folder to save analyzed *.mat files');
  if (saveFolder == 0)
    warning('No *.mat files will be saved - Save folder not selected');
  else
    [parentPath, ~, ~] = parsePath(saveFolder);
  end
end

% Select folder of Calcium dFoF exported from ImageJ
if isempty(CaFolder) && ~param.reAnalyzeOption
  CaFolder = uigetdir(parentPath, 'Select folder of dFoF *.csv files exported from ImageJ');
  if (CaFolder == 0) return; end
  [parentPath, ~, ~] = parsePath(CaFolder);
end

% Prompt for timing file
if isempty(timingFile) && ~param.reAnalyzeOption
  [fileName, filePath] = uigetfile('.csv', 'Select the corresponding timing.csv file', parentPath);
  timingFile = strcat(filePath,fileName);
  if ~all(timingFile) return; end
end

% Select export folder of Ca events (if selected)
if isempty(expCaFolder) && param.expCaEvOption
  expCaFolder = uigetdir(parentPath, 'Select folder to export Ca event *.csv files');
  if (expCaFolder == 0) warning('No Ca event files to be exported - folder not selected'); end
  [parentPath, ~, ~] = parsePath(expCaFolder);
end

% Select export folder of SWR events (if selected)
if isempty(expSWRFolder) && param.expSWREvOption
  expSWRFolder = uigetdir(parentPath, 'Select folder to export updated SWR event *.csv files');
  if (expSWRFolder == 0) warning('No updated SWR event files to be exported - folder not selected'); end
  [parentPath, ~, ~] = parsePath(expSWRFolder);
end

% Select export folder of stim events (if selected)
if isempty(expStimFolder) && param.expStimEvOption
  expStimFolder = uigetdir(parentPath, 'Select folder to export stimulation event *.csv files');
  if (expStimFolder == 0) warning('No stimulation event files to be exported - folder not selected'); end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract matlab data file names (if options require)
if param.swrCaOption || param.spkCaOption || param.stimCaOption || param.reAnalyzeOption
  cd (dataFolder);
  dir_temp   = dir('*.mat'); % Find only *.mat files
  names      = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  dataFiles  = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nDataFiles = length(dataFiles);
end

% Extract dFoF data file names (if options require)
if ~param.reAnalyzeOption
  cd (CaFolder);
  dir_temp   = dir('*.csv'); % Find only *.csv files
  names      = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  CaFiles    = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nCaFiles   = length(CaFiles);
end

cd (curPath);

% Error handling - file number mismatch:
if ~param.reAnalyzeOption && (param.swrCaOption || param.spkCaOption || param.stimCaOption)
  if (nDataFiles ~= nCaFiles)
    error('Unequal number of files in data and Ca folders - analysis will be mismatched');
  end
  nFiles = nDataFiles;
elseif param.reAnalyzeOption
  nFiles = nDataFiles;
else
  nFiles = nCaFiles;
end

% Determine individual file names:
dataFile{nFiles}    = [];
CaFile{nFiles}      = [];
saveFile{nFiles}    = [];
expCaFile{nFiles}   = [];
expSWRFile{nFiles}  = [];
expStimFile{nFiles} = [];

for i = 1:nFiles
  
  % Determine dataFile (if options require)
  if param.reAnalyzeOption || param.swrCaOption || param.spkCaOption || param.stimCaOption
    dataFile{i} = [dataFolder slash dataFiles{i}];
    [~, dataFileName, ~] = parsePath(dataFile{i});
  end
  
  % Determine CaFile (if options require)
  if ~param.reAnalyzeOption
    CaFile{i} = [CaFolder slash CaFiles{i}];
    if ~param.swrCaOption && ~param.spkCaOption && ~param.stimCaOption [~, dataFileName, ~] = parsePath(CaFile{i}); end
  end
  
  % Determine individual output *.mat file names
  if (saveFolder ~= 0)
    saveFile{i} = [saveFolder slash dataFileName '.mat'];
  else
    saveFile{i} = 0;
  end
  
  % Determine individual exported Ca event *.csv file names (if selected)
  if ~isempty(expCaFolder) && param.expCaEvOption
    if (expCaFolder ~= 0)
      expCaFile{i} = [expCaFolder slash dataFileName '_CaEvents.csv'];
    end
  end
  
  % Determine individual exported SWR event *.csv file names (if selected)
  if ~isempty(expSWRFolder) && param.expSWREvOption
    if (expSWRFolder ~= 0)
      expSWRFile{i} = [expSWRFolder slash dataFileName '_swrEvents.csv'];
    end
  end
  
  % Determine individual exported stim event *.csv file names (if selected)
  if ~isempty(expStimFolder) && param.expStimEvOption
    if (expStimFolder ~= 0)
      expStimFile{i} = [expStimFolder slash dataFileName '_stimEvents.csv'];
    end
  end
  
end

%% Pre-process Ca data:
% Only run if options require:
if ~param.reAnalyzeOption || (param.baseCorrectMethod > 0) || param.interpOption
  
  % Initialize local parameters to not broadcast to parfor loop:
  reAnalyzeOption = param.reAnalyzeOption;
  swrCaOption     = param.swrCaOption;
  spkCaOption     = param.spkCaOption;
  stimCaOption    = param.stimCaOption;
  
  parfor i = 1:nFiles
    if reAnalyzeOption || swrCaOption || spkCaOption || stimCaOption
      data = load(dataFile{i});
      processCaFile(data, [], param, saveFile{i}, CaFile{i}, timingFile);
    else
      processCaFile([], [], param, saveFile{i}, CaFile{i}, timingFile);
    end
  end
  
  % Safegaurd: Pause for 60 seconds to wait for all files to save:
  pause(60);
end

%% Calculate thresholds:
baseMean{nFiles}   = [];
baseSD{nFiles}     = [];
baseThresh{nFiles} = [];
peakThresh{nFiles} = [];

% Only calculate if enabled (otherwise thresholds must have been previously calculated)
if param.baseDetectMethod > 0
  
  % If consistent threshold must open each file first then calculate
  if param.consThreshOption
    tSeries{nFiles}   = [];
    
    for i = 1:nFiles
      data = load(saveFile{i});
      if stimCaOption
        tSeries{i} = trimCaFile(data.Ca.tSeries, data.Ca.samplingInt, data.LFP.samplingInt, data.stim.evStart, param);
      else
        tSeries{i} = trimCaFile(data.Ca.tSeries, data.Ca.samplingInt, [], [], param);
      end
    end
    [baseMn, baseSd, baseTh, peakTh] = calcCaThresh(tSeries, param);
    baseMean(:)   = {baseMn};
    baseSD(:)     = {baseSd};
    baseThresh(:) = {baseTh};
    peakThresh(:) = {peakTh};
    
  % Otherwise calculate each threshold individually
  else
    stimCaOption = param.stimCaOption;
    parfor i = 1:nFiles
      data = load(saveFile{i});
      if stimCaOption
        tSeries = trimCaFile(data.Ca.tSeries, data.Ca.samplingInt, data.LFP.samplingInt, data.stim.evStart, param);
      else
        tSeries = trimCaFile(data.Ca.tSeries, data.Ca.samplingInt, [], [], param);
      end
      [baseMean{i}, baseSD{i}, baseThresh{i}, peakThresh{i}] = calcCaThresh(tSeries, param);
    end
  end
end

%% Analyze Ca files:
reAnalyzeOption = param.reAnalyzeOption;
baseCorrectMethod = param.baseCorrectMethod;
interpOption = param.interpOption;
parfor i = 1:nFiles
  if ~reAnalyzeOption || (baseCorrectMethod > 0) || interpOption
    data = load(saveFile{i});
  else
    data = load(dataFile{i});
  end
  data.Ca.baseMean   = baseMean{i};
  data.Ca.baseSD     = baseSD{i};
  data.Ca.baseThresh = baseThresh{i};
  data.Ca.peakThresh = peakThresh{i};
  analyzeCaFile(data, [], param, saveFile{i}, expCaFile{i}, expSWRFile{i}, expStimFile{i});
end
fprintf('complete\n');
end