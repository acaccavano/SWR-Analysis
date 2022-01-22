function analyzeLFPBatch(param, dataFolder, saveFolder, expEvFolder, expDataFolder, stimFolder, expAveFile)
%% analyzeLFPBatch(param, dataFolder, saveFolder, expEvFolder, expDataFolder, stimFolder, expAveFile)
%
%  Function to run analyzeLFPFile on batch of files
%
%   param      = structure containing all parameters including:
%     param.fileNum          = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.fileType         = 1 = pClamp (.abf), 2 = ASCII data (folder ofdata files), 3 = Matlab (.mat)
%     param.Fs               = sampling rate (ASCII recordings are usually 3000, not needed for pClamp files)
%     param.dsFactor         = downsample factor (default = 1, no downsampling)
%     param.lfpChannel       = channel to use for LFP input (default = 1, but depends on recording)
%     param.cellOption       = boolean flag to determine if second cell channel to be imported
%     param.cellChannel      = channel to use for optional cell input (default = 2, but depends on recording)
%     param.filtType         = Option for filtering (1 = built-in MATLAB bandpass (default) or 2 = custom gaussian (perfect phase-response but computationally expensive and can create artifacts)
%     param.notchOption      = option to perform comb filter to remove electrical line noise (default = 0)
%     param.notchFreq        = frequency to remove (+harmonics) (default = 60Hz)
%     param.lfpOption        = boolean flag to filter LFP signal
%     param.lfpLim1          = lower LFP band-pass lim (default = 1Hz)
%     param.lfpLim2          = upper LFP band-pass lim (default = 1000Hz)
%     param.swrOption        = boolean flag to detect SWR events
%     param.swOption         = boolean flag to filter and analyze SW signal
%     param.swLim1           = lower sharp wave band-pass lim (default = 1Hz)
%     param.swLim2           = upper sharp wave band-pass lim (default = 30Hz)
%     param.rmsPeriodSW      = root-mean square window [ms] (in Eschenko 2008 = 5ms), but had more luck with longer ~25ms
%     param.rOption          = boolean flag to filter and analyze ripple signal
%     param.rLim1            = lower ripple band-pass lim (default = 120)
%     param.rLim2            = upper ripple band-pass lim (default = 220)
%     param.rmsPeriodR       = root-mean square window [ms] (in Eschenko 2008 = 5ms)
%     param.baseDetectMethod = Method for baseline stats detection (0: none, 1: lower quantile, 2: iterative gaussian fitting (default))
%     param.baseQuant        = Lower quantile for baseline cutoff (default = 0.95)
%     param.skewedBL         = boolean option to indicate skewed BL distribution, and use both gaussians just for BL 
%     param.pkDiffMin        = min distance between double gaussian peaks to consider them equivalent = abs(B1-B2) (default = 0.01 RMS)
%     param.pkSimLim         = Peak amplitude similarity metric = (A1^2 + A2^2)/(A1*A2) (default = 2)
%     param.kurtosisMin      = Min kurtosis limit to fit with 2 gaussians (otherwise skip 1st fit) (default = 0)
%     param.kurtosisMax      = Max kurtosis limit until exclude high points (otherwise fit can fail) (default = 5)
%     param.excludeQuant     = quantile above which to exclude if max kurtosis limit reached (default = 0.98)
%     param.plotFitHisto     = boolean option to plot histograms and fits for each file
%     param.peakDetectOption = boolean flag to detect SW and ripple reaks in RMS signals
%     param.rmsMinEvDiff     = min difference between detected RMS peaks [ms] (in Eschenko 2008 = 25ms), had more luck with longer ~100ms, but may cut off doublets
%     param.rmsMinEvDur      = min duration of RMS peaks [ms]
%     param.sdMultSW         = SD of baseline for threshold detection (default = 4)
%     param.sdBaseFactorSW   = Factor of sdMult to consider for event start/end times (default = 0.5 eg 2SD)
%     param.sdMultR          = SD of baseline for threshold detection (default = 4)
%     param.sdBaseFactorR    = Factor of sdMult to consider for event start/end times (default = 0.5 eg 2SD)
%     param.swrType          = Option to determine what qualifies as SWR (1: SW & R (default), 2: SW only, 3: R only)
%     param.swrWindow        = +/- window around SWR peak events for swrData file [ms]
%     param.expSWREvOption   = boolean flag to determine whether to export csv table of SWR events
%     param.expSWRDataOption = boolean flag to determine whether to export txt file of episodic SWR events for pClamp analysis
%     param.thetaOption      = boolean flag to filter and analyze theta signal
%     param.thetaLim1        = lower theta band-pass lim (default = 4Hz)
%     param.thetaLim2        = upper theta band-pass lim (default = 8Hz)
%     param.alphaOption      = boolean flag to filter and analyze alpha signal
%     param.alphaLim1        = lower alpha band-pass lim (default = 9Hz)
%     param.alphaLim2        = upper alpha band-pass lim (default = 12Hz)
%     param.betaOption       = boolean flag to filter and analyze beta signal
%     param.betaLim1         = lower beta band-pass lim (default = 13Hz)
%     param.betaLim2         = upper beta band-pass lim (default = 24Hz)
%     param.gammaOption      = boolean flag to filter and analyze gamma signal
%     param.gammaLim1        = lower gamma band-pass lim (default = 25Hz)
%     param.gammaLim2        = upper gamma band-pass lim (default = 55Hz)
%     param.hgammaOption     = boolean flag to filter and analyze high gamma signal
%     param.hgammaLim1       = lower high gamma band-pass lim (default = 65Hz)
%     param.hgammaLim2       = upper high gamma band-pass lim (default = 85Hz)
%     param.fROption         = boolean flag to filter and analyze fast ripple signal
%     param.fRLim1           = lower fast rippple band-pass lim (default = 250Hz)
%     param.fRLim2           = lower fast rippple band-pass lim (default = 500Hz)
%     param.spectOption      = boolean flag to calculate spectrogram
%     param.spectLim1        = lower lim of spectrogram (default = 1Hz)
%     param.spectLim2        = upper lim of spectrogram (default = 500Hz)
%     param.fftOption        = boolean flag to calculate FFT
%     param.phaseOption      = boolean flag to calculate piecewise linear interpolated phase (required for many LFP cross frequency, spike-phase, and PSC-LFP correlation analyses
%     param.xFreqOption      = boolean flag to perform cross-frequency analysis
%     param.xFreqBin         = Frequency bin size for n x n PAC analysis (Default = 5 Hz)
%     param.xFreqLow         = cell: low frequency band for x-freq (Theta, Alpha, Beta, SW)
%     param.nShuffle         = # shuffles to calculate Z-value for total PAC - does not due for nxn or time PAC, very computationally expensive. (default = 200)
%     param.morlWidth        = width/number of cycles of the morlet wavelet filter, default = 7
%     param.winLength        = time binning for phase-amplitude analysis (s). Dictates min low freq (=1/winLength), so default = 0.5s results in min freq. of 2Hz
%     param.winOverlap       = Amount to overlap time bins (default = 0.2s)
%     param.importStimOption = option to import stim file from pClamp (default = 0)
%     param.reAnalyzeOption  = option to re-analyze file - will prompt for *.mat instead of raw data file
%     param.expAveOption     = boolean flag to determine whether to export csv table of average statistics
%     param.transposeOption  = boolean flag to transpose exported average stats from row to column format
%   dataFolder    = full path to folder containing raw data to be analysed (if not set, will prompt)
%   saveFolder    = full path to folder of matlab files to save (if not set, will prompt)
%   expEvFolder   = full path to folder of exported csv event table (if not set and expSWREvOption = 1, will prompt)
%   expDataFolder = full path to folder of exported txt data file (if not set and expSWRDataOption = 1, will prompt)
%   stimFolder    = full path to folder of pClamp stim events (if not set and importStimOption = 1, will prompt)
%   expAveFile    = full path to file of exported csv table of averages (if not set and expAveOption = 1, will prompt)

%% Handle optional arguments
if (nargin < 7); expAveFile  = []; end
if (nargin < 6); stimFolder    = []; end
if (nargin < 5); expDataFolder = []; end
if (nargin < 4); expEvFolder   = []; end
if (nargin < 3); saveFolder    = []; end
if (nargin < 2); dataFolder    = []; end
if (nargin < 1); param         = struct; end

% Handle case in which empty variable is supplied:
if isempty(param); param       = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum');          param.fileNum           = 2;    end
if ~isfield(param,'fileType');         param.fileType          = 1;    end
if ~isfield(param,'Fs');               param.Fs                = 3000; end  % [Hz]
if ~isfield(param,'dsFactor');         param.dsFactor          = 1;    end
if ~isfield(param,'lfpChannel');       param.lfpChannel        = 1;    end
if ~isfield(param,'cellOption');       param.cellOption        = 1;    end
if ~isfield(param,'cellChannel');      param.cellChannel       = 2;    end
if ~isfield(param,'filtType');         param.filtType          = 1;    end
if ~isfield(param,'notchOption');      param.notchOption       = 0;    end
if ~isfield(param,'notchFreq');        param.notchFreq         = 60;   end  % [Hz]
if ~isfield(param,'lfpOption');        param.lfpOption         = 1;    end
if ~isfield(param,'lfpLim1');          param.lfpLim1           = 1;    end  % [Hz]
if ~isfield(param,'lfpLim2');          param.lfpLim2           = 1000; end  % [Hz]
if ~isfield(param,'swrOption');        param.swrOption         = 1;    end
if ~isfield(param,'swOption');         param.swOption          = 1;    end
if ~isfield(param,'swLim1');           param.swLim1            = 1;    end  % [Hz]
if ~isfield(param,'swLim2');           param.swLim2            = 30;   end  % [Hz]
if ~isfield(param,'rmsPeriodSW');      param.rmsPeriodSW       = 25;   end  % [ms]
if ~isfield(param,'rOption');          param.rOption           = 1;    end
if ~isfield(param,'rLim1');            param.rLim1             = 120;  end  % [Hz]
if ~isfield(param,'rLim2');            param.rLim2             = 220;  end  % [Hz]
if ~isfield(param,'rmsPeriodR');       param.rmsPeriodR        = 5;    end  % [ms]
if ~isfield(param,'baseDetectMethod'); param.baseDetectMethod  = 2;    end
if ~isfield(param,'baseQuant');        param.baseQuant         = 0.80; end
if ~isfield(param,'skewedBL');         param.skewedBL          = 1;    end
if ~isfield(param,'pkDiffMin');        param.pkDiffMin         = 0.1;  end 
if ~isfield(param,'pkSimLim');         param.pkSimLim          = 2;    end
if ~isfield(param,'kurtosisMin');      param.kurtosisMin       = 0;    end
if ~isfield(param,'kurtosisMax');      param.kurtosisMax       = 5;    end
if ~isfield(param,'excludeQuant');     param.excludeQuant      = 0.95; end
if ~isfield(param,'plotFitHisto');     param.plotFitHisto      = 0;    end
if ~isfield(param,'peakDetectOption'); param.peakDetectOption  = 1;    end
if ~isfield(param,'rmsMinEvDiff');     param.rmsMinEvDiff      = 100;  end  % [ms]
if ~isfield(param,'rmsMinEvDur');      param.rmsMinEvDur       = 25;   end  % [ms]
if ~isfield(param,'sdMultSW');         param.sdMultSW          = 4;    end
if ~isfield(param,'sdBaseFactorSW');   param.sdBaseFactorSW    = 0.5;  end
if ~isfield(param,'sdMultR');          param.sdMultR           = 4;    end
if ~isfield(param,'sdBaseFactorR');    param.sdBaseFactorR     = 0.5;  end
if ~isfield(param,'swrType');          param.swrType           = 1;    end
if ~isfield(param,'swrWindow');        param.swrWindow         = 100;  end
if ~isfield(param,'expSWREvOption');   param.expSWREvOption    = 1;    end
if ~isfield(param,'expSWRDataOption'); param.expSWRDataOption  = 1;    end
if ~isfield(param,'thetaOption');      param.thetaOption       = 1;    end
if ~isfield(param,'thetaLim1');        param.thetaLim1         = 4;    end  % [Hz]
if ~isfield(param,'thetaLim2');        param.thetaLim2         = 8;    end  % [Hz]
if ~isfield(param,'alphaOption');      param.alphaOption       = 0;    end
if ~isfield(param,'alphaLim1');        param.alphaLim1         = 8;    end  % [Hz]
if ~isfield(param,'alphaLim2');        param.alphaLim2         = 12;   end  % [Hz]
if ~isfield(param,'betaOption');       param.betaOption        = 0;    end
if ~isfield(param,'betaLim1');         param.betaLim1          = 12;   end  % [Hz]
if ~isfield(param,'betaLim2');         param.betaLim2          = 20;   end  % [Hz]
if ~isfield(param,'gammaOption');      param.gammaOption       = 1;    end
if ~isfield(param,'gammaLim1');        param.gammaLim1         = 20;   end  % [Hz]
if ~isfield(param,'gammaLim2');        param.gammaLim2         = 50;   end  % [Hz]
if ~isfield(param,'hgammaOption');     param.hgammaOption      = 0;    end
if ~isfield(param,'hgammaLim1');       param.hgammaLim1        = 65;   end  % [Hz]
if ~isfield(param,'hgammaLim2');       param.hgammaLim2        = 85;   end  % [Hz]
if ~isfield(param,'fROption');         param.fROption          = 1;    end
if ~isfield(param,'fRLim1');           param.fRLim1            = 250;  end  % [Hz]
if ~isfield(param,'fRLim2');           param.fRLim2            = 500;  end  % [Hz]
if ~isfield(param,'spectOption');      param.spectOption       = 1;    end
if ~isfield(param,'spectLim1');        param.spectLim1         = 1;    end  % [Hz]
if ~isfield(param,'spectLim2');        param.spectLim2         = 500;  end  % [Hz]
if ~isfield(param,'fftOption');        param.fftOption         = 1;    end
if ~isfield(param,'phaseOption');      param.phaseOption       = 1;    end
if ~isfield(param,'xFreqOption');      param.xFreqOption       = 1;    end
if ~isfield(param,'xFreqBin');         param.xFreqBin          = 5;    end  % [Hz]
if ~isfield(param,'xFreqLow');         param.xFreqLow          = 'Theta'; end
if ~isfield(param,'morlWidth');        param.morlWidth         = 7;    end
if ~isfield(param,'nShuffle');         param.nShuffle          = 200;  end 
if ~isfield(param,'winLength');        param.winLength         = 0.5;  end  % [s]
if ~isfield(param,'winOverlap');       param.winOverlap        = 0.2;  end  % [s]
if ~isfield(param,'importStimOption'); param.importStimOption  = 0;    end
if ~isfield(param,'reAnalyzeOption');  param.reAnalyzeOption   = 0;    end
if ~isfield(param,'expAveOption');     param.expAveOption      = 1;    end
if ~isfield(param,'transposeOption');  param.transposeOption   = 0;    end 

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% If not supplied, prompt for folders to analyze
if isempty(dataFolder)
  if param.reAnalyzeOption
    dataFolder = uigetdir(pwd, 'Select folder containing analyzed LFP *.mat files');
  elseif (param.fileType == 1)
    dataFolder = uigetdir(pwd, 'Select folder containing set of pClamp *.abf files');
  elseif (param.fileType == 2)
    dataFolder = uigetdir(pwd, 'Select folder containing subfolders of recordings');
  elseif (param.fileType == 3)
    dataFolder = uigetdir(pwd, 'Select folder containing preprocessed single-channel LFP *.mat files');
  end
end
if (dataFolder == 0); return; end

% Parse dataFolder to determine default save name
[parentPath, ~, ~] = parsePath(dataFolder);

% Select folders to save analyzed matlab files
if isempty(saveFolder)
  saveFolder = uigetdir(parentPath, 'Select folder to save analyzed *.mat files');
  if (saveFolder == 0)
    warning('No *.mat files will be saved - Save folder not selected');
  else
    [parentPath, ~, ~] = parsePath(saveFolder);
  end
end

% Select folder to export SWR event files, if option selected
if isempty(expEvFolder) && param.expSWREvOption
  expEvFolder = uigetdir(parentPath, 'Select folder to export SWR event *.csv files');
  if (expEvFolder == 0)
    warning('No files to be exported - SWR event folder not selected');
  else
    [parentPath, ~, ~] = parsePath(expEvFolder);
  end
end

% Select folder to export SWR event-locked episodic data files, if option selected
if isempty(expDataFolder) && param.expSWRDataOption
  expDataFolder = uigetdir(parentPath, 'Select folder to export SWR event-locked episodic data *.txt files');
  if (expDataFolder == 0) 
    warning('No files to be exported - SWR data folder not selected'); 
  else
    [parentPath, ~, ~] = parsePath(expDataFolder);
  end
end

% Select folder of stimulation events, if option selected
if isempty(stimFolder) && param.importStimOption
  stimFolder = uigetdir(parentPath, 'Select folder of stimulation event *.csv files');
  if (stimFolder == 0); error('No stimulation folder selected'); end
end

% Select file to export table of average stats, if option selected
if isempty(expAveFile) && param.expAveOption
  [~, dataFolderName, ~] = parsePath(parentPath(1:end-1));
  defaultPath = [parentPath dataFolderName '_aveStats.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of average statistics', defaultPath);
  expAveFile = [exportPath exportName];
  if ~all(expAveFile); warning('No average statistics to be exported - no file selected'); end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract data file names
cd (dataFolder);
if param.reAnalyzeOption
  dir_temp = dir('*.mat'); % Find only matlab files
  names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
elseif (param.fileType == 1)
  dir_temp = dir('*.abf'); % Find only abf files
  names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
elseif (param.fileType == 2)
  dir_temp = dir(); % Find all files/folders in directory
  dir_flag = [dir_temp.isdir] & ~strcmp({dir_temp.name},'.') & ~strcmp({dir_temp.name},'..'); % assign flag only for subfolders
  dir_sub = dir_temp(dir_flag);
  file = {dir_sub.name};
elseif (param.fileType == 3)
  dir_temp = dir('*.mat'); % Find only mat files
  names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
end
nDataFiles = length(file);

% Extract stim event files names (if option selected)
if param.importStimOption
  cd (stimFolder);
  dir_temp   = dir('*.csv'); % Find only *.csv files
  names      = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  stimFiles   = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nStimFiles  = length(stimFiles);
  
  % Check for same number of data and stim files:
  if (nDataFiles ~= nStimFiles)
    error('Unequal number of files in data and stim folders - analysis will be mismatched');
  end
end

cd (curPath);

% Determine individual file names:
dataFile{nDataFiles}    = [];
saveFile{nDataFiles}    = [];
expEvFile{nDataFiles}   = [];
expDataFile{nDataFiles} = [];
stimFile{nDataFiles}    = [];

for i = 1:nDataFiles
  % Determine individual import files/folders
  dataFile{i}  = [dataFolder slash file{i}];
  [~, dataFileName, ~] = parsePath(dataFile{i});
  
  % Determine individual output *.mat file names
  if (saveFolder ~= 0)
    saveFile{i} = [saveFolder slash dataFileName '.mat'];
  else
    saveFile{i} = 0;
  end
  
  % Determine individual exported event *.csv file names (if selected)
  if ~isempty(expEvFolder) && param.expSWREvOption
    if (expEvFolder ~= 0)
      expEvFile{i} = [expEvFolder slash dataFileName '_swrEvents.csv'];
    end
  end
  
  % Determine individual exported SWR event-locked episodic data *.txt file names (if selected)
  if ~isempty(expDataFolder) && param.expSWRDataOption
    if (expDataFolder ~= 0)
      expDataFile{i} = [expDataFolder slash dataFileName '_swrData.txt'];
    end
  end
  
  % Determine individual imported stim file names (if selected)
  if param.importStimOption
    stimFile{i} = [stimFolder slash stimFiles{i}];
  end
    
end

% Initialize aveStats Cell Array and Table:
if param.expAveOption
  varNames{nDataFiles} = [];
  aveStats = table;
end

reAnalyzeOption = param.reAnalyzeOption;
parfor i = 1:nDataFiles
  
  if reAnalyzeOption
    data = load(dataFile{i});
    [~, ~, aveStats(i, :), varNames{i}] = analyzeLFPFile(data, [], param, dataFile{i}, saveFile{i}, expEvFile{i}, expDataFile{i}, stimFile{i}, expAveFile);
  else
    [~, ~, aveStats(i, :), varNames{i}] = analyzeLFPFile([], [], param, dataFile{i}, saveFile{i}, expEvFile{i}, expDataFile{i}, stimFile{i}, expAveFile);
  end

end

% Write expAveStats file:
if param.expAveOption
  fileNames = aveStats{:,1};
  
  % Convert to cell array and replace NaN values with blanks:
  tmpCell = table2cell(aveStats(:, 2:end)); % Excluding filename
  tmpCell(isnan(aveStats(:, 2:end).Variables)) = {[]};
  
  % Convert to array and transpose table:
  if param.transposeOption
    aveStats = cell2table(tmpCell');
    aveStats.Properties.RowNames       = varNames{1}(2:end);
    aveStats.Properties.VariableNames  = fileNames;
    aveStats.Properties.DimensionNames = {'Variable', 'File'};
  else
    aveStats = cell2table(tmpCell);
    aveStats.Properties.RowNames       = fileNames;
    aveStats.Properties.VariableNames  = varNames{1}(2:end);
    aveStats.Properties.DimensionNames = {'File', 'Variable'};
  end
  writetable(aveStats, expAveFile, 'Delimiter', ',', 'WriteVariableNames', true, 'WriteRowNames', true);
end

fprintf('complete\n');

end
