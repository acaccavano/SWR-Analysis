function analyzeCaBatch(param, dataFolder, saveFolder, CaFolder, timingFile, expCaFolder, expSWRFolder, expStimFolder)
%% analyzeCaBatch(param, dataFolder, saveFolder, CaFolder, timingFile, expCaFolder, expSWRFolder, expStimFolder)
%
%  Function to run processCaFile, calcCaThresh, and analyzeCaFile on batch of files
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   param      = structure containing all parameters including:
%     param.fileNum              = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.interpOption         = boolean flag to interpolate file (needed if comparing to LFP) (default = 1)
%     param.samplingInt          = interpolated sampling interval (default = 0.5ms)
%     param.baseCorrectMethod    = Method for baseline correction (0: none, 1: gassuian filter, 2: smoothed average (default))
%     param.CaFiltLim1           = Lower limit for gaussian filter (default = 0.03Hz)
%     param.CaFiltLim2           = Upper limit for gaussian filter (default = 4Hz)
%     param.CaFiltOrder          = Gaussian filter order (default = 80)
%     param.CaFiltAlpha          = Gaussian filter alpha (default = 2.5)
%     param.smoothFactor         = Proportion of file duration for moving linear average (default = 0.25)
%     param.peakDetectCa         = boolean option to detect calcium events (default = 1)
%     param.baseDetectMethod     = Method for baseline stats detection (0: none, 1: lower quantile, 2: iterative gaussian fitting (default))
%     param.baseQuant            = Lower quantile for baseline cutoff (default = 0.8)
%     param.pkDiffMin            = min distance between double gaussian peaks to consider them equivalent = abs(B1-B2) (default = 0.1)
%     param.pkSimLim             = Peak amplitude similarity metric = (A1^2 + A2^2)/(A1*A2) (default = 2)
%     param.kurtosisMin          = Min kurtosis limit to fit with 2 gaussians (otherwise skip 1st fit) (default = 0)
%     param.kurtosisMax          = Max kurtosis limit until exclude high points (otherwise fit can fail) (default = 5)
%     param.excludeQuant         = quantile above which to exclude if max kurtosis limit reached (default = 0.98)
%     param.plotFitHisto         = boolean option to plot histograms and fits for each cell
%     param.sdMult               = SD of baseline for threshold detection (default = 4)
%     param.sdBaseFactor         = Factor of sdMult to consider for event start/end times (default = 0.75 eg 3SD)
%     param.skipDetectLim        = Skip detection for first duration of recording for uncorrected photobleaching (default = 1s)
%     param.consThreshOption     = option to calculate same threshold for multiple files (default = 0)
%     param.swrCaOption          = option to perform coincidence detection for SWRs and Ca transients (default = 1)
%     param.useSWRDurationOption = option to use detected SWR detection for coincidence detection (default = 1)
%     param.useSWRWindowOption   = option to use standard swrWindow for coincidence detection (default = 0)
%     param.swrWindow            = +/- window around SWR peak events (default = 100 ms)
%     param.expCaEvOption        = option to export csv table of Calcium events (default = 1)
%     param.expSWREvOption       = option to export csv table of SWR events (default = 1)
%     param.spkCaOption          = option to perform coincidence detection for SWRs and Ca transients (default = 0, placeholder: code not written yet)
%     param.stimCaOption         = option to perform coincidence detection for Stim and Ca transients (default = 0)
%     param.stimCaLim1           = time after stim start to start stim window (default = 0ms)
%     param.stimCaLim2           = time after stim start to end stim window (default = 1000ms)
%     param.expStimEvOption      = option to export csv table of stim events (default = 0)
%     param.reAnalyzeOption      = option to re-analyze file (default = 0)
%   dataFolder    = full path to matlab folder to import (if not set, will prompt)
%   saveFolder    = full path to matlab folder to save (can be same, if not set, will prompt)
%   CaFolder      = full path to dFoF csv folder exported from ImageJ (if not set, will prompt)
%   timingFile    = full path to timing csv file previously setup (if not set, will prompt)
%   expCaFolder   = full path to calcium event csv folder to export (if not set, will prompt)
%   expSWRFolder  = full path to SWR event csv folder to export (if not set, will prompt)
%   expStimFolder = full path to stim event csv folder to export (if not set, will prompt)

%% Handle input arguments
if (nargin < 8); expStimFolder = []; end
if (nargin < 7); expSWRFolder  = []; end
if (nargin < 6); expCaFolder   = []; end
if (nargin < 5); timingFile    = []; end
if (nargin < 4); CaFolder      = []; end
if (nargin < 3); saveFolder    = []; end
if (nargin < 2); dataFolder    = []; end
if (nargin < 1); param         = struct; end

% Handle case in which empty variable is supplied:
if isempty(param); param      = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum');              param.fileNum              = 2;    end
if ~isfield(param,'baseCorrectMethod');    param.baseCorrectMethod    = 2;    end
if ~isfield(param,'CaFiltLim1');           param.CaFiltLim1           = 0.03; end
if ~isfield(param,'CaFiltLim2');           param.CaFiltLim2           = 4;    end
if ~isfield(param,'CaFiltOrder');          param.CaFiltOrder          = 80;   end
if ~isfield(param,'CaFiltAlpha');          param.CaFiltAlpha          = 2.5;  end
if ~isfield(param,'smoothFactor');         param.smoothFactor         = 0.25; end
if ~isfield(param,'interpOption');         param.interpOption         = 1;    end
if ~isfield(param,'samplingInt');          param.samplingInt          = 0.5;  end
if ~isfield(param,'cellTypeOption');       param.cellTypeOption       = 1;    end
if ~isfield(param,'nCellTypes');           param.nCellTypes           = 2;    end
if ~isfield(param,'cellTypeName');         param.cellTypeName{param.nCellTypes} = []; end
if ~isfield(param,'peakDetectCa');         param.peakDetectCa         = 1;    end
if ~isfield(param,'baseDetectMethod');     param.baseDetectMethod     = 2;    end
if ~isfield(param,'baseQuant');            param.baseQuant            = 0.8;  end
if ~isfield(param,'pkDiffMin');            param.pkDiffMin            = 0.1;  end 
if ~isfield(param,'pkSimLim');             param.pkSimLim             = 2;    end
if ~isfield(param,'kurtosisMin');          param.kurtosisMin          = 0;    end
if ~isfield(param,'kurtosisMax');          param.kurtosisMax          = 5;    end
if ~isfield(param,'excludeQuant');         param.excludeQuant         = 0.98; end
if ~isfield(param,'plotFitHisto');         param.plotFitHisto         = 0;    end
if ~isfield(param,'sdMult');               param.sdMult               = 4;    end
if ~isfield(param,'sdBaseFactor');         param.sdBaseFactor         = 0.75; end
if ~isfield(param,'skipDetectLim');        param.skipDetectLim        = 1;    end
if ~isfield(param,'consThreshOption');     param.consThreshOption     = 0;    end
if ~isfield(param,'expCaEvOption');        param.expCaEvOption        = 1;    end
if ~isfield(param,'swrCaOption');          param.swrCaOption          = 1;    end
if ~isfield(param,'useSWRDurationOption'); param.useSWRDurationOption = 1;    end
if ~isfield(param,'useSWRWindowOption');   param.useSWRWindowOption   = 0;    end
if ~isfield(param,'swrWindow');            param.swrWindow            = 100;  end
if ~isfield(param,'expSWREvOption');       param.expSWREvOption       = 0;    end
if ~isfield(param,'alignEndOption');       param.alignEndOption       = 0;    end
if ~isfield(param,'stimCaOption');         param.stimCaOption         = 0;    end
if ~isfield(param,'stimCaLim1');           param.stimCaLim1           = 0;    end
if ~isfield(param,'stimCaLim2');           param.stimCaLim2           = 1000; end
if ~isfield(param,'expStimEvOption');      param.expStimEvOption      = 0;    end
if ~isfield(param,'reAnalyzeOption');      param.reAnalyzeOption      = 0;    end
    
% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% Prompt first for matlab data folder of analyzed files, if options require:
if isempty(dataFolder)
  if param.reAnalyzeOption
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed dFoF *.mat files');
    if (dataFolder == 0); return; end
    [parentPath, ~, ~] = parsePath(dataFolder);
  elseif (param.swrCaOption || param.stimCaOption)
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed LFP and/or cell *.mat files');
    if (dataFolder == 0); return; end
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
  if (CaFolder == 0); return; end
  [parentPath, ~, ~] = parsePath(CaFolder);
end

% Prompt for timing file
if isempty(timingFile) && ~param.reAnalyzeOption
  [fileName, filePath] = uigetfile('.csv', 'Select the corresponding timing.csv file', parentPath);
  timingFile = strcat(filePath,fileName);
  if ~all(timingFile); return; end
end

% Select export folder of Ca events (if selected)
if isempty(expCaFolder) && param.expCaEvOption
  expCaFolder = uigetdir(parentPath, 'Select folder to export Ca event *.csv files');
  if (expCaFolder == 0); warning('No Ca event files to be exported - folder not selected'); end
  [parentPath, ~, ~] = parsePath(expCaFolder);
end

% Select export folder of SWR events (if selected)
if isempty(expSWRFolder) && param.expSWREvOption
  expSWRFolder = uigetdir(parentPath, 'Select folder to export updated SWR event *.csv files');
  if (expSWRFolder == 0); warning('No updated SWR event files to be exported - folder not selected'); end
  [parentPath, ~, ~] = parsePath(expSWRFolder);
end

% Select export folder of stim events (if selected)
if isempty(expStimFolder) && param.expStimEvOption
  expStimFolder = uigetdir(parentPath, 'Select folder to export stimulation event *.csv files');
  if (expStimFolder == 0); warning('No stimulation event files to be exported - folder not selected'); end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract matlab data file names (if options require)
if param.swrCaOption || param.stimCaOption || param.reAnalyzeOption
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
if ~param.reAnalyzeOption && (param.swrCaOption || param.stimCaOption)
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
  if param.reAnalyzeOption || param.swrCaOption || param.stimCaOption
    dataFile{i} = [dataFolder slash dataFiles{i}];
    [~, dataFileName, ~] = parsePath(dataFile{i});
  end
  
  % Determine CaFile (if options require)
  if ~param.reAnalyzeOption
    CaFile{i} = [CaFolder slash CaFiles{i}];
    if (~param.swrCaOption && ~param.stimCaOption); [~, dataFileName, ~] = parsePath(CaFile{i}); end
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
  stimCaOption    = param.stimCaOption;
  
  parfor i = 1:nFiles
    if reAnalyzeOption || swrCaOption || stimCaOption
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