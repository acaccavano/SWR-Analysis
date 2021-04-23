function [data, hand] = processCaFile(data, hand, param, saveFile, CaFile, timingFile)
%% [data, hand] = processCaFile(data, hand, param, saveFile, CaFile, timingFile)
%
%  Function to import dFoF, correct baseline, and interpolate files to prepare for event detection
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   data       = data structure containing analyzed LFP files
%   hand       = handle structure to specify where figure should be drawn
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
%     param.excludeQuant         = quantile above which to exclude if max kurtosis limit reached
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
%   saveFile    = full path to matlab file to save (if not set, will prompt)
%   CaFile      = full path to dFoF csv file exported from ImageJ (if not set, will prompt)
%   timingFile  = full path to timing csv file previously setup (if not set, will prompt)
%
%  Outputs:
%   data       = structure containing all data to be saved
%   hand       = handle structure for figure

%% Handle input arguments - if not entered
if (nargin < 6); timingFile = []; end
if (nargin < 5); CaFile     = []; end
if (nargin < 4); saveFile   = []; end
if (nargin < 3); param      = struct; end
if (nargin < 2); hand       = struct; end
if (nargin < 1); data       = struct; end

% Handle case in which empty variables are supplied:
if isempty(param); param    = struct; end
if isempty(hand);  hand     = struct; end
if isempty(data);  data     = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum');              param.fileNum              = 1;    end
if ~isfield(param,'baseCorrectMethod');    param.baseCorrectMethod    = 2;    end
if ~isfield(param,'CaFiltLim1');           param.CaFiltLim1           = 0.03; end
if ~isfield(param,'CaFiltLim2');           param.CaFiltLim2           = 4;    end
if ~isfield(param,'CaFiltOrder');          param.CaFiltOrder          = 80;   end
if ~isfield(param,'CaFiltAlpha');          param.CaFiltAlpha          = 2.5;  end
if ~isfield(param,'smoothFactor');         param.smoothFactor         = 0.25; end
if ~isfield(param,'interpOption');         param.interpOption         = 1;    end
if ~isfield(param,'samplingInt');          param.samplingInt          = 0.5;  end
if ~isfield(param,'cellTypeOption');       param.cellTypeOption       = 0;    end
if ~isfield(param,'cellType1');            param.cellType1          = 'Deep'; end
if ~isfield(param,'cellType2');            param.cellType2          = 'Supe'; end
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

% Re-initialize path variables:
parentPath   = [];
dataFileName = [];

% Import previously analyzed matlab file if options require it:
if (param.swrCaOption && ~isfield(data, 'SWR')) || (param.stimCaOption && ~isfield(data, 'stim'))
  [fileName, filePath] = uigetfile('.mat', 'Select *.mat file of analyzed LFP and/or cell channel(s)');
  dataFile = [filePath fileName];
  if ~all(dataFile)
    error('No previously analyzed *.mat file selected');
  else
    [parentPath, dataFileName, ~] = parsePath(dataFile);
  end
  data = load(dataFile);
end

% Check if necessary data structures are present
if ~isfield(data,'Ca'); data.Ca = struct; end

if param.swrCaOption || param.stimCaOption
  if ~isfield(data,'LFP'); error('Must analyze LFP before proceeding'); end
end

if param.swrCaOption
  if ~isfield(data,'SWR'); error('Must analyze LFP channel for SWR events before proceeding'); end
end

if param.stimCaOption
  if ~isfield(data,'stim'); error('Must analyze LFP and import stim events before proceeding'); end
end
    
% If not supplied, prompt for dFoF Ca file to analyze
if isempty(CaFile) && ~param.reAnalyzeOption
  
  % Option for when *.mat data file has already been imported (for SWR/Spk correlation)
  if ~isempty(parentPath)
    [fileName, filePath] = uigetfile('.csv', 'Select dFoF *.csv file exported from ImageJ', parentPath);
    CaFile = strcat(filePath,fileName);
    if ~all(CaFile) 
      return;
    else
      [parentPath, ~, ~] = parsePath(CaFile);
    end
    
  % Option for when *.mat data file has not been imported (pure Ca analysis)
  else
    [fileName, filePath] = uigetfile('.csv', 'Select dFoF *.csv file exported from ImageJ');
    CaFile = strcat(filePath,fileName);
    if ~all(CaFile)
      return;
    else
      [parentPath, dataFileName, ~] = parsePath(CaFile);
    end
  end
end

% If not supplied, prompt for timing file
if isempty(timingFile) && ~param.reAnalyzeOption
  [fileName, filePath] = uigetfile('.csv', 'Select the corresponding timing.csv file', parentPath);
  timingFile = strcat(filePath,fileName);
  if ~all(timingFile); return; end
end

% If not supplied, prompt for save file
if isempty(saveFile)
  defaultName = [parentPath dataFileName '.mat'];
  [saveName, savePath] = uiputfile('.mat','Select file to save output matlab file', defaultName);
  saveFile = strcat(savePath, saveName);
  if (saveFile == 0); warning('No *.mat file will be saved - Save file not selected'); end
end

if (saveFile ~=0)
  data.saveFile = saveFile;
  [~, dataFileName, saveExt] = parsePath(data.saveFile);
  data.saveName  = [dataFileName '.' saveExt];
end

%% Import data
if ~param.reAnalyzeOption
  CaTable             = readtable(CaFile,'Delimiter',',','ReadRowNames',1);
  timingTable         = readtable(timingFile,'ReadVariableNames',0);
  data.Ca.tSeries     = CaTable{:,:};
  data.Ca.timing      = round(1000 * timingTable{:,:}); % Convert from s to ms
  data.Ca.tSeriesR    = data.Ca.tSeries; % Static raw data variable that won't change with subsequent processing
  data.Ca.timingR     = data.Ca.timing; % Static raw data variable that won't change with subsequent processing
  data.Ca.samplingInt = data.Ca.timing(2) - data.Ca.timing(1);
  data.Ca.nChannels   = size(data.Ca.tSeries, 2);
  data.Ca.CaFile      = CaFile;
  data.Ca.timingFile  = timingFile;
  
  % Import cell type from variable names (last four characters of dFoF file column headings)
  if param.cellTypeOption
    fprintf(['sorting by arrays by cell type (file ' dataFileName ')... ']);
    
    data.Ca.cellTypeRaw{data.Ca.nChannels} = [];
    for i = 1:data.Ca.nChannels
      colName = CaTable.Properties.VariableNames{i};
      data.Ca.cellTypeRaw{i} = colName(length(colName) - 3 : length(colName));
    end
    
    % Sort cells by cell-type
    data.Ca.sortInd  = horzcat(find(strcmp(data.Ca.cellTypeRaw, param.cellType1)), find(strcmp(data.Ca.cellTypeRaw, param.cellType2)));
    data.Ca.tSeriesS = data.Ca.tSeriesR(:, data.Ca.sortInd);
    
    % Update tSeries and cellType variables
    data.Ca.tSeries  = zeros(size(data.Ca.tSeriesS, 1), size(data.Ca.tSeriesS, 2));
    data.Ca.tSeries  = data.Ca.tSeriesS;
    data.Ca.cellType = data.Ca.cellTypeRaw(data.Ca.sortInd);
    fprintf('done\n');
  end
end

data.param     = param;
data.Ca.param  = param; % Save to Ca structure, as subsequent analysis may alter data.param

%% Baseline Correction of dFoF (if selected)
if (param.baseCorrectMethod > 0)
  fprintf(['Correcting baseline dFoF (file ' dataFileName ')... ']);
  
  % Reset samplingInt (important if re-analyzing)
  data.Ca.samplingInt = data.Ca.timingR(2) - data.Ca.timingR(1);
  
  % Create static corrected variable that won't change with subsequent processing
  data.Ca.tSeriesF = zeros(size(data.Ca.tSeriesR, 1), size(data.Ca.tSeriesR, 2));
  
  for i = 1 : data.Ca.nChannels
    data.Ca.tSeriesF(:,i) = correctBaseline(data.Ca.timingR, data.Ca.tSeries(:,i), param);
  end
  
  % Update tSeries variable
  data.Ca.tSeries = zeros(size(data.Ca.tSeriesF, 1), size(data.Ca.tSeriesF, 2));
  data.Ca.tSeries = data.Ca.tSeriesF;
  fprintf('done\n');
end

%% Interpolate dFoF (if selected)
if param.interpOption
  fprintf(['interpolating dFoF (file ' dataFileName ')... ']);
  
  % Assign new timing array
  data.Ca.samplingInt = param.samplingInt;
  data.Ca.timing      = (0 : data.Ca.samplingInt : data.Ca.timingR(length(data.Ca.timingR)))';
  
  % Create static interpolated variable that won't change with subsequent processing
  data.Ca.tSeriesI = zeros(length(data.Ca.timing), size(data.Ca.tSeries, 2));

  % Interpolate data
  data.Ca.tSeriesI = interp1(data.Ca.timingR, data.Ca.tSeries, data.Ca.timing, 'linear', 'extrap');
  
  % Update tSeries variable
  data.Ca.tSeries  = zeros(size(data.Ca.tSeriesI, 1), size(data.Ca.tSeriesI, 2));
  data.Ca.tSeries  = data.Ca.tSeriesI;
  
  fprintf('done\n');
end

%% Re-order structure arrays
data       = orderfields(data);
data.param = orderfields(data.param);
data.Ca    = orderfields(data.Ca);
data.Ca.param = orderfields(data.Ca.param);

%% Save matlab file
if all(data.saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(data.saveFile,'-struct','data');
  fprintf('done\n');
end

end