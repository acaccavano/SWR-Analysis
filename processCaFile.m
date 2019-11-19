function [data, hand] = processCaFile(data, hand, param, saveFile, CaFile, timingFile)
%% [data, hand] = importCaFile(data, hand, param, saveFile, CaFile, timingFile, exportFile)

%% Handle input arguments - if not entered
if (nargin < 6) timingFile = []; end
if (nargin < 5) CaFile     = []; end
if (nargin < 4) saveFile   = []; end
if (nargin < 3) param      = struct; end
if (nargin < 2) hand       = struct; end
if (nargin < 1) data       = struct; end

% Handle case in which empty variables are supplied:
if isempty(param) param    = struct; end
if isempty(hand)  hand     = struct; end
if isempty(data)  data     = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum')              param.fileNum              = 1;   end
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
if ~isfield(param,'stimCaLim1')           param.stimCaLim1           = 0;   end % ms
if ~isfield(param,'stimCaLim2')           param.stimCaLim2           = 1000; end % ms
if ~isfield(param,'expStimEvOption')      param.expStimEvOption      = 0;   end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption      = 0;   end

% Re-initialize path variables:
parentPath   = [];
dataFileName = [];

% Import previously analyzed matlab file if options require it:
if (param.swrCaOption && ~isfield(data, 'SWR')) || (param.spkCaOption && ~isfield(data, 'C')) || (param.stimCaOption && ~isfield(data, 'stim'))
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
if ~isfield(data,'Ca') data.Ca = struct; end

if param.swrCaOption || param.spkCaOption || param.stimCaOption
  if ~isfield(data,'LFP') error('Must analyze LFP before proceeding'); end
end

if param.swrCaOption
  if ~isfield(data,'SWR') error('Must analyze LFP channel for SWR events before proceeding'); end
end

if param.spkCaOption
  if ~isfield(data,'C') error('Must analyze cell channel before proceeding'); end
end

if param.stimCaOption
  if ~isfield(data,'stim') error('Must analyze LFP and import stim events before proceeding'); end
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
  if ~all(timingFile) return; end
end

% If not supplied, prompt for save file
if isempty(saveFile)
  defaultName = [parentPath dataFileName '.mat'];
  [saveName, savePath] = uiputfile('.mat','Select file to save output matlab file', defaultName);
  saveFile = strcat(savePath, saveName);
  if (saveFile == 0) warning('No *.mat file will be saved - Save file not selected'); end
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
    data.Ca.tSeriesF(:,i) = correctBaseline(data.Ca.timingR, data.Ca.tSeriesR(:,i), param);
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
  data.Ca.tSeriesI = zeros(size(data.Ca.timing, 1), size(data.Ca.tSeries, 2));

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