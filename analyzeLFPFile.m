function [data, hand, aveStats, varNames] = analyzeLFPFile(data, hand, param, dataFile, saveFile, expEvFile, expDataFile, stimFile, expAveFile)
%% [data, hand, aveStats, varNames] = analyzeLFPFile(data, hand, param, dataFile, saveFile, expEvFile, expDataFile, stimFile, expAveFile)
%
%  Function to detect sharp wave ripple (SWR) events, theta, beta, and gamma analysis, time-frequency 
%  spectrogram analysis, and/or stimulation event pre-processing of single LFP recording. SWRs are detected
%  by performing both a low frequency (sharp wave) and high frequency (ripple) band-pass filter of 
%  an LFP recording. The root-mean-square (RMS) of these filtered channels are taken, and events
%  are counted if they exceed a given standard deviation of the RMS channels. SWR events are counted if 
%  both a sharp wave and ripple occur simultaneously.
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   data       = structure - specify if appending previously analyzed files
%   hand       = handle structure to specify where figure should be drawn
%   param      = structure containing all parameters including:
%     param.fileNum          = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.fileType         = 1 = pClamp (.abf), 2 = ASCII data (folder of data files), 3 = Matlab (.mat)
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
%   dataFile    = full path to file/folder containing data to be analysed (if not set, will prompt)
%   saveFile    = full path to matlab file to save (if not set, will prompt)
%   expEvFile   = full path to exported csv event table (if not set and expSWREvOption = 1, will prompt
%   expDataFile = full path to exported txt data file (if not set and expSWRDataOption = 1, will prompt
%   stimFile    = full path to pClamp stim event file (if not set and importStimOption = 1, will prompt
%   expAveFile  = full path to exported csv average table (if not set and expAveOption = 1, will prompt
%
%  Outputs:
%   data     = structure containing all data to be saved
%   hand     = handle structure for figure
%   aveStats = 1D array of average statistics (size variable, depending on options selected)
%   varNames = 1D array of average statistics (size variable, depending on options selected)

%% Handle input arguments - if not entered
if (nargin < 9); expAveFile  = []; end
if (nargin < 8); stimFile    = []; end
if (nargin < 7); expDataFile = []; end
if (nargin < 6); expEvFile   = []; end
if (nargin < 5); saveFile    = []; end
if (nargin < 4); dataFile    = []; end
if (nargin < 3); param       = struct; end
if (nargin < 2); hand        = struct; end
if (nargin < 1); data        = struct; end

% Handle case in which empty variables are supplied:
if isempty(param); param     = struct; end
if isempty(hand);  hand      = struct; end
if isempty(data);  data      = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum');          param.fileNum           = 1;    end
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

% Initialize LFP structure if it doesn't already exist
if ~isfield(data,'LFP'); data.LFP = struct; end

% If cell option is selected, initialize cell structure if it doesn't already exist
if param.cellOption
  if ~isfield(data,'C'); data.C = struct; end
end

% If not supplied, prompt for files/folder to analyze
if isfield(data.LFP, 'dataFile')
  dataFile = data.LFP.dataFile;
elseif isempty(dataFile)
  if (param.fileType == 1)
    [fileName, filePath] = uigetfile('.abf', 'Select *.abf file to analyze');
    dataFile = strcat(filePath, fileName);
  elseif (param.fileType == 2)
    dataFile = uigetdir();
  elseif (param.fileType == 3)
    [fileName, filePath] = uigetfile('.mat', 'Select *.mat file to analyze');
    dataFile = strcat(filePath, fileName);
  end
  if ~all(dataFile); return; end
end

% Parse dataFile to determine default save name
[parentPath, dataFileName, ~] = parsePath(dataFile);

% If not supplied, prompt for save file
if isempty(saveFile)
  defaultPath = [parentPath dataFileName '.mat'];
  [data.saveName, savePath] = uiputfile('.mat','Select file to save output matlab file', defaultPath);
  saveFile = [savePath data.saveName];
  if ~all(saveFile)
    warning('No file to be saved - no file selected');
  else
    [parentPath, dataFileName, ~] = parsePath(saveFile);
  end
end

% Select export SWR event file, if option selected
if isempty(expEvFile) && param.expSWREvOption
  defaultPath = [parentPath dataFileName '_swrEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of SWR events', defaultPath);
  expEvFile = [exportPath exportName];
  if ~all(expEvFile)
    warning('No SWR events to be exported - no file selected');
  else
    [parentPath, ~, ~] = parsePath(expEvFile);
  end
end

% Select export SWR data file, if option selected
if isempty(expDataFile) && param.expSWRDataOption
  defaultPath = [parentPath dataFileName '_swrData.txt'];
  [exportName, exportPath] = uiputfile('.txt','Select *.txt file to export episodic SWR data', defaultPath);
  expDataFile = [exportPath exportName];
  if ~all(expDataFile) 
    warning('No SWR events to be exported - no file selected'); 
  else
    [parentPath, ~, ~] = parsePath(expDataFile);
  end
end

% Select stimulation event file, if option selected
if param.importStimOption
  if isempty(stimFile)
    [stimName, stimPath] = uigetfile('.csv', 'Select stimulation event *.csv file from pClamp', parentPath);
    stimFile = [stimPath stimName];
    if ~all(stimFile); error('No stimulation file selected'); end
  end
  [~, stimFileName, ~] = parsePath(stimFile);
end

% Select export average statistics file, if option selected
if isempty(expAveFile) && param.expAveOption
  defaultPath = [parentPath dataFileName '_aveStats.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of average statistics', defaultPath);
  expAveFile = [exportPath exportName];
  if ~all(expAveFile); warning('No average statistics to be exported - no file selected'); end
end

%% Import data
if ~isfield(data.LFP, 'dataFile')
  fprintf(['importing file ' dataFileName '... ']);
  if (param.fileType == 1) % abf file
    [dataIn, samplingInt, ~] = abfload(dataFile);

    data.LFP.samplingInt = samplingInt;
    data.LFP.tSeriesRaw  = dataIn(:, param.lfpChannel);
    data.LFP.samplingInt = data.LFP.samplingInt / 1000; % convert from um to ms
    
    if param.cellOption
      data.C.samplingInt = samplingInt;
      data.C.tSeries = dataIn(:, param.cellChannel);
      data.C.samplingInt = data.C.samplingInt / 1000; % convert from um to ms
    end
    
  elseif (param.fileType == 2) % ASCII files
    
    % Ensure current directory is in path so helper functions work
    curPath = pwd;
    path(path, curPath);
    
    % Extract file names
    cd (dataFile);
    dir_temp = dir;
    names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
    files = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
    
    % Sort file names
    if iscell(files)
      for i = 1:length(files)
        temp = files{i};
        if ~contains(temp,'_')
          string = temp(1:strfind(temp,'.')-1);
        else
          if ~contains(temp,'.')
            string = temp(max(strfind(temp,'_'))+1:end);
          else
            string = temp(max(strfind(temp,'_'))+1:strfind(temp,'.')-1);
          end
        end
        % extract the number part
        str_ascii = str2num(int2str(string)); %#ok<*ST2NM>
        index_nums = find(str_ascii <=57 & str_ascii >=48);
        filenum(i) = str2num(string(index_nums)); %#ok<FNDSB,AGROW>
      end
      [filenum, index]=sort(filenum,'ascend'); %#ok<ASGLU>
      files = files(index);
    else
      temp_file = files;
      clear files;
      files = cell(1);
      files{1}=temp_file;
    end
    
    % Initialize cell arrays:
    data.nFiles = length(files);
    dataIn{data.nFiles} = []; % raw data cell array
    
    %% Read Data file by file
    for i = 1:data.nFiles
      buf = fopen(files{i},'r'); % open data file was files{i}
      firstline = fgetl(buf);
      while (isempty(firstline))
        firstline = fgetl(buf);
      end
      n = length(sscanf(firstline,'%g'));
      data_firstline = sscanf(firstline,'%g',[n inf]);
      
      dataTemp = fscanf(buf,'%g',[n inf]);
      dataTemp = [data_firstline dataTemp]; %#ok<AGROW>
      dataTemp = dataTemp';
      
      fclose(buf);  % close data file after read
      
      % Import LFP channel
      dataIn{i} = dataTemp(:, param.lfpChannel);
      
      % Import cell channel if selected
      if param.cellOption
        dataIn{i} = horzcat(dataIn{i}, dataTemp(:, param.cellChannel));
      end
    end
    
    % Concatenate data
    dataIn = vertcat(dataIn{:});
    data.LFP.tSeriesRaw  = dataIn(:,1);
    data.LFP.samplingInt = 1000 / param.Fs; % (ms)
    if param.cellOption
      data.C.tSeries = dataIn(:,2);
      data.C.samplingInt = 1000 / param.Fs; % (ms)
    end
    
  elseif (param.fileType == 3) % mat file
    inStruct = load(dataFile);
    varNames = fieldnames(inStruct);
    data.LFP.tSeriesRaw  = inStruct.(varNames{1})/1000;
    data.LFP.samplingInt = 1000 / param.Fs; % (ms)
    
  end
  fprintf('done\n');
  
  %% Assign parameters to structure array
  data.LFP.dataFile  = dataFile;
  
  % Downsample data if selected
  if (param.dsFactor > 1)
    fprintf(['downsampling by factor of ' num2str(param.dsFactor) ' (file ' dataFileName ')... ']);
    
    % Display warning if downsampling results in non-integer Fs:
    FsTemp = 1000 / (data.LFP.samplingInt * param.dsFactor);
    if (round(FsTemp,3) ~= round(FsTemp))
      warning('Downsampling to non-integer sampling rate, may cause timing mismatch')
    end
    
    data.LFP.samplingInt  = data.LFP.samplingInt * param.dsFactor;
    data.LFP.tSeriesRaw   = downsampleMean(data.LFP.tSeriesRaw, param.dsFactor);

    if param.cellOption
      data.C.samplingInt  = data.C.samplingInt * param.dsFactor;
      data.C.tSeries      = downsampleMean(data.C.tSeries, param.dsFactor);
    end
    
    fprintf('done\n');
  end

  % Determine timing:
  data.LFP.nSamples = length(data.LFP.tSeriesRaw);
  data.LFP.timing   = (0: data.LFP.samplingInt : (data.LFP.nSamples-1) * data.LFP.samplingInt)';
  
  if param.cellOption
    data.C.nSamples = length(data.C.tSeries);
    data.C.timing   = (0: data.C.samplingInt : (data.C.nSamples-1) * data.C.samplingInt)';
  end

end

data.saveFile  = saveFile;
[~, saveName, saveExt] = parsePath(saveFile);
data.saveName  = [saveName '.' saveExt];

% Assign Fs based on samplingInt imported (may differ from input Fs)
param.Fs = 1000 / data.LFP.samplingInt;

data.param     = param;
data.LFP.param = param; % Save to LFP structure, as subsequent analysis may alter data.param

%% Filter data
if param.notchOption
  fprintf(['Notch filter %4.1fHz noise (file ' dataFileName ')... '], param.notchFreq);
  nHarm =   5;  % Number of harmonics
  Q     = 100;  % Quality Factor to determine bandwidth
  
  for i = 1:nHarm
    W0 = (i * param.notchFreq / (round(param.Fs)/2));
    BW = W0/Q;
    [b, a] = iirnotch(W0, BW);
    data.LFP.tSeries = filtfilt(b, a, data.LFP.tSeriesRaw);
  end
  fprintf('done\n');
else
  data.LFP.tSeries = data.LFP.tSeriesRaw;
end

% If filtering with MATLAB bandpass, append signal ends with flipped and reversed data to eliminate edge effects (already handled with custom gaussian filter)
if (param.filtType == 1)
  tSeriesApp1 = 2*data.LFP.tSeries(1) - flipud(data.LFP.tSeries(2 : round(param.Fs) + 1));
  tSeriesApp2 = 2*data.LFP.tSeries(end) - flipud(data.LFP.tSeries(end - round(param.Fs) : end - 1));
  tSeriesApp  = [tSeriesApp1; data.LFP.tSeries; tSeriesApp2];
end

if param.lfpOption
  % Apply filter to LFP signal for DC drift and HF noise
  fprintf(['band-pass filtering LFP between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.lfpLim1, param.lfpLim2);
  if (param.filtType == 1)
    data.LFP.tSeries = bandpass(tSeriesApp, [param.lfpLim1 param.lfpLim2], round(param.Fs));
    data.LFP.tSeries = data.LFP.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.LFP.tSeries = gaussianFilt(data.LFP.tSeries, param.lfpLim1, param.lfpLim2, data.LFP.samplingInt, 20);
  end
  data.LFP.tPower  = bandpower(data.LFP.tSeries);
  data.LFP.lim1    = param.lfpLim1;
  data.LFP.lim2    = param.lfpLim2;
end

if param.swOption
  % Apply filter to extract SW signal
  fprintf(['band-pass filtering sharp wave between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.swLim1, param.swLim2);
  if ~isfield(data,'SW'); data.SW = struct; end
  if (param.filtType == 1)
    data.SW.tSeries = bandpass(tSeriesApp, [param.swLim1 param.swLim2], round(param.Fs));
    data.SW.tSeries = data.SW.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.SW.tSeries = gaussianFilt(data.LFP.tSeries, param.swLim1, param.swLim2, data.LFP.samplingInt, 20);
  end
  data.SW.tPower  = bandpower(data.SW.tSeries);
  data.SW.lim1    = param.swLim1;
  data.SW.lim2    = param.swLim2;
  fprintf('done\n');
end

if param.rOption
  % Apply filter to extract ripple signal
  fprintf(['band-pass filtering ripple between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.rLim1, param.rLim2);
  if ~isfield(data,'R'); data.R = struct; end
  if (param.filtType == 1)
    data.R.tSeries = bandpass(tSeriesApp, [param.rLim1 param.rLim2], round(param.Fs));
    data.R.tSeries = data.R.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.R.tSeries = gaussianFilt(data.LFP.tSeries, param.rLim1, param.rLim2, data.LFP.samplingInt, 1);
  end
  data.R.tPower  = bandpower(data.R.tSeries);
  data.R.lim1    = param.rLim1;
  data.R.lim2    = param.rLim2;
  fprintf('done\n');
end

if param.thetaOption
  % Apply Gaussian filter to extract theta signal
  fprintf(['band-pass filtering theta between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.thetaLim1, param.thetaLim2);
  if ~isfield(data,'theta'); data.theta = struct; end
  if (param.filtType == 1)
    data.theta.tSeries = bandpass(tSeriesApp, [param.thetaLim1 param.thetaLim2], round(param.Fs));
    data.theta.tSeries = data.theta.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.theta.tSeries = gaussianFilt(data.LFP.tSeries, param.thetaLim1, param.thetaLim2, data.LFP.samplingInt, 2);
  end
  data.theta.tPower  = bandpower(data.theta.tSeries);
  data.theta.lim1    = param.thetaLim1;
  data.theta.lim2    = param.thetaLim2;
  fprintf('done\n');
end

if param.alphaOption
  % Apply Gaussian filter to extract alpha signal
  fprintf(['band-pass filtering alpha between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.alphaLim1, param.alphaLim2);
  if ~isfield(data,'alpha'); data.alpha = struct; end
  if (param.filtType == 1)
    data.alpha.tSeries = bandpass(tSeriesApp, [param.alphaLim1 param.alphaLim2], round(param.Fs));
    data.alpha.tSeries = data.alpha.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.alpha.tSeries = gaussianFilt(data.LFP.tSeries, param.alphaLim1, param.alphaLim2, data.LFP.samplingInt, 2);
  end
  data.alpha.tPower  = bandpower(data.alpha.tSeries);
  data.alpha.lim1    = param.alphaLim1;
  data.alpha.lim2    = param.alphaLim2;
  fprintf('done\n');
end

if param.betaOption
  % Apply Gaussian filter to extract beta signal
  fprintf(['band-pass filtering beta between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.betaLim1, param.betaLim2);
  if ~isfield(data,'beta'); data.beta = struct; end
  if (param.filtType == 1)
    data.beta.tSeries = bandpass(tSeriesApp, [param.betaLim1 param.betaLim2], round(param.Fs));
    data.beta.tSeries = data.beta.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.beta.tSeries = gaussianFilt(data.LFP.tSeries, param.betaLim1, param.betaLim2, data.LFP.samplingInt, 1);
  end
  data.beta.tPower  = bandpower(data.beta.tSeries);
  data.beta.lim1    = param.betaLim1;
  data.beta.lim2    = param.betaLim2;
  fprintf('done\n');
end

if param.gammaOption
  % Apply Gaussian filter to extract gamma signal
  fprintf(['band-pass filtering gamma between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.gammaLim1, param.gammaLim2);
  if ~isfield(data,'gamma'); data.gamma = struct; end
  if (param.filtType == 1)
    data.gamma.tSeries = bandpass(tSeriesApp, [param.gammaLim1 param.gammaLim2], round(param.Fs));
    data.gamma.tSeries = data.gamma.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.gamma.tSeries = gaussianFilt(data.LFP.tSeries, param.gammaLim1, param.gammaLim2, data.LFP.samplingInt, 1);
  end
  data.gamma.tPower  = bandpower(data.gamma.tSeries);
  data.gamma.lim1    = param.gammaLim1;
  data.gamma.lim2    = param.gammaLim2;
  fprintf('done\n');
end

if param.hgammaOption
  % Apply Gaussian filter to extract high gamma signal
  fprintf(['band-pass filtering high gamma between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.hgammaLim1, param.hgammaLim2);
  if ~isfield(data,'hgamma'); data.hgamma = struct; end
  if (param.filtType == 1)
    data.hgamma.tSeries = bandpass(tSeriesApp, [param.hgammaLim1 param.hgammaLim2], round(param.Fs));
    data.hgamma.tSeries = data.hgamma.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.hgamma.tSeries = gaussianFilt(data.LFP.tSeries, param.hgammaLim1, param.hgammaLim2, data.LFP.samplingInt, 1);
  end
  data.hgamma.tPower  = bandpower(data.hgamma.tSeries);
  data.hgamma.lim1    = param.hgammaLim1;
  data.hgamma.lim2    = param.hgammaLim2;
  fprintf('done\n');
end

if param.fROption
  % Apply Gaussian filter to extract fast ripple signal
  fprintf(['band-pass filtering fast ripple between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.fRLim1, param.fRLim2);
  if ~isfield(data,'fR'); data.fR = struct; end
  if (param.filtType == 1)
    data.fR.tSeries = bandpass(tSeriesApp, [param.fRLim1 param.fRLim2], round(param.Fs));
    data.fR.tSeries = data.fR.tSeries(round(param.Fs) + 1 : end - round(param.Fs));
  elseif (param.filtType == 2)
    data.fR.tSeries = gaussianFilt(data.LFP.tSeries, param.fRLim1, param.fRLim2, data.LFP.samplingInt, 1);
  end
  data.fR.tPower  = bandpower(data.fR.tSeries);
  data.fR.lim1    = param.fRLim1;
  data.fR.lim2    = param.fRLim2;
  fprintf('done\n');
end


%% SWR event detection if both SW and ripple option enabled
if param.swrOption
  if ~isfield(data,'SWR'); data.SWR = struct; end
  
  % RMS Signal calculation of SW signal
  if isfield(data,'SW')
    fprintf(['calculating root mean square (RMS) of SW in a %4.1fms sliding window (file ' dataFileName ')... '], param.rmsPeriodSW);
    rmsFactor = round(param.rmsPeriodSW / data.LFP.samplingInt);
    nSamplesRMS = floor(data.LFP.nSamples / rmsFactor);
    data.SW.RMS = zeros(nSamplesRMS, 1);
    
    for i = 1:nSamplesRMS
      loAve          = max((i - 2) * rmsFactor, 1);
      hiAve          = min(i * rmsFactor, data.LFP.nSamples);
      data.SW.RMS(i) = rms(data.SW.tSeries(loAve:hiAve));
    end
    data.param.rmsPeriodSW = rmsFactor * data.LFP.samplingInt; % Corrected rmsPeriod
    timingRMS        = (0: data.param.rmsPeriodSW: (nSamplesRMS-1) * data.param.rmsPeriodSW)';
    data.SW.RMS      = interp1(timingRMS, data.SW.RMS, data.LFP.timing, 'spline');
    fprintf('done\n');
  end
  
  % RMS Signal calculation of ripple signal
  if isfield(data,'R')
    fprintf(['calculating root mean square (RMS) of ripple in a %4.1fms sliding window (file ' dataFileName ')... '], param.rmsPeriodR);
    rmsFactor = round(param.rmsPeriodR / data.LFP.samplingInt);
    nSamplesRMS = floor(data.LFP.nSamples / rmsFactor);
    data.R.RMS  = zeros(nSamplesRMS, 1);
    
    for i = 1:nSamplesRMS
      loAve          = max((i - 2) * rmsFactor, 1);
      hiAve          = min(i * rmsFactor, data.LFP.nSamples);
      data.R.RMS(i)  = rms(data.R.tSeries(loAve:hiAve));
    end
    data.param.rmsPeriodR = rmsFactor * data.LFP.samplingInt; % Corrected rmsPeriod
    timingRMS        = (0: data.param.rmsPeriodR: (nSamplesRMS-1) * data.param.rmsPeriodR)';
    data.R.RMS       = interp1(timingRMS, data.R.RMS, data.LFP.timing, 'spline');
    fprintf('done\n');
  end

  
  %% Event Detection
  if param.peakDetectOption
    %% Find sharp wave peak based on standard deviation of RMS of SW signal
    if isfield(data,'SW')
      fprintf(['detecting SW events %4.0f standard deviations above baseline (file ' dataFileName ')... '], param.sdMultSW);
      
      % Re-initialize data structures
      data.SW.evStatus  = [];
      data.SW.evStart   = [];
      data.SW.evPeak    = [];
      data.SW.evEnd     = [];
      data.SW.IEI       = [];
      data.SW.power     = [];
      data.SW.duration  = [];
      data.SW.frequency = 0;
      
      % Estimate baseline based on options selected, either taking param.baseQuant quantile of signal
      % (unreliable for active/quiet recordings, or through an iterative gaussian fitting process (more robust)
      [mn, sd, hand] = calcBaseline(data.SW.RMS, hand, param);
      data.SW.peakThresh = mn + sd * param.sdMultSW;
      data.SW.baseThresh = mn + param.sdBaseFactorSW * sd * param.sdMultSW;
      [data.SW.evStatus, data.SW.evStart, data.SW.evPeak, data.SW.evEnd] = peakFindUnique(data.SW.RMS, data.LFP.timing, data.SW.peakThresh, data.SW.baseThresh, 1, param.rmsMinEvDiff, param.rmsMinEvDur);
      
      if ~isnan(data.SW.evStart)
        for i = 1:length(data.SW.evStart)
          data.SW.power(i)    = bandpower(data.SW.tSeries(data.SW.evStart(i) : data.SW.evEnd(i)));
          data.SW.duration(i) = (data.LFP.timing(data.SW.evEnd(i)) - data.LFP.timing(data.SW.evStart(i)));
          if (i > 1); data.SW.IEI = horzcat(data.SW.IEI, (data.LFP.timing(data.SW.evPeak(i)) - data.LFP.timing(data.SW.evPeak(i-1))) / 1000); end
        end
        data.SW.frequency = length(data.SW.evStart) / ((data.LFP.timing(length(data.LFP.timing)) - data.LFP.timing(1)) / 1000);
        data.SW.power     = data.SW.power';
        data.SW.duration  = data.SW.duration';
        data.SW.IEI       = data.SW.IEI';
      end
    end
    
    %% Find ripple peak based on standard deviation of RMS of ripple signal
    if isfield(data,'R')
      fprintf(['detecting ripple events %4.0f standard deviations above baseline (file ' dataFileName ')... '], param.sdMultR);
      
      % Re-initialize data structures
      data.R.evStatus  = [];
      data.R.evStart   = [];
      data.R.evPeak    = [];
      data.R.evEnd     = [];
      data.R.IEI       = [];
      data.R.power     = [];
      data.R.duration  = [];
      data.R.frequency = 0;
      
      % Estimate baseline based on options selected, either taking param.baseQuant quantile of signal
      % (unreliable for active/quiet recordings, or through an iterative gaussian fitting process
      [mn, sd, hand] = calcBaseline(data.R.RMS, hand, param);
      data.R.peakThresh = mn + sd * param.sdMultR;
      data.R.baseThresh = mn + param.sdBaseFactorR * sd * param.sdMultR;
      [data.R.evStatus, data.R.evStart, data.R.evPeak, data.R.evEnd] = peakFindUnique(data.R.RMS, data.LFP.timing, data.R.peakThresh, data.R.baseThresh, 1, param.rmsMinEvDiff, param.rmsMinEvDur);
      
      if ~isnan(data.R.evStart)
        for i = 1:length(data.R.evStart)
          data.R.power(i)    = bandpower(data.R.tSeries(data.R.evStart(i) : data.R.evEnd(i)));
          data.R.duration(i) = (data.LFP.timing(data.R.evEnd(i)) - data.LFP.timing(data.R.evStart(i)));
          if (i > 1); data.R.IEI = horzcat(data.R.IEI, (data.LFP.timing(data.R.evPeak(i)) - data.LFP.timing(data.R.evPeak(i-1))) / 1000); end
        end
        data.R.frequency = length(data.R.evStart) / ((data.LFP.timing(length(data.LFP.timing)) - data.LFP.timing(1)) / 1000);
        data.R.power     = data.R.power';
        data.R.duration  = data.R.duration';
        data.R.IEI       = data.R.IEI';
      end
    end
    
    %% SWR Event Calculation
    
    % (Re)Initialize data arrays
    data.SWR.evStatus  = [];
    data.SWR.evStart   = [];
    data.SWR.evPeak    = [];
    data.SWR.evEnd     = [];
    data.SWR.IEI       = [];
    data.SWR.duration  = [];
    data.SWR.amp       = [];
    data.SWR.power     = [];
    data.SWR.area      = [];
    data.SWR.event     = [];  
    data.SWR.frequency = 0;
    
    % SW arrays:
    if isfield(data,'SW')
      if ~isfield(data.SW,'SWR'); data.SW.SWR = struct; end
      data.SW.SWR.event = [];
      data.SW.SWR.power = [];
      data.SW.SWR.area  = [];
    end
    
    % Ripple arrays:
    if isfield(data,'R')
      if ~isfield(data.R,'SWR'); data.R.SWR = struct; end
      data.R.SWR.event  = [];
      data.R.SWR.power  = [];
    end

    % Gamma arrays:
    if isfield(data,'gamma')
      if ~isfield(data.gamma,'SWR'); data.gamma.SWR = struct; end
      data.gamma.SWR.event = [];
      data.gamma.SWR.power = [];
    end

    % High gamma arrays:
    if isfield(data,'hgamma')
      if ~isfield(data.hgamma,'SWR'); data.hgamma.SWR = struct; end
      data.hgamma.SWR.event = [];
      data.hgamma.SWR.power = [];
    end
    
    % Fast ripple arrays:
    if isfield(data,'fR')
      if ~isfield(data.fR,'SWR'); data.fR.SWR = struct; end
      data.fR.SWR.event = [];
      data.fR.SWR.power = [];
    end
    
    if param.swrType == 1 % Overlap of SW and ripple
      [data.SWR.evStatus, data.SWR.evStart, data.SWR.evEnd, data.SWR.evIndex] = eventOverlap(data.SW.evStatus, data.SW.evStart, data.SW.evEnd, data.R.evStatus, data.R.evStart, data.R.evEnd, data.LFP.timing, 0);
      
    elseif param.swrType == 2 % SW only
      data.SWR.evStatus = data.SW.evStatus;
      data.SWR.evStart  = data.SW.evStart;
      data.SWR.evEnd    = data.SW.evEnd;
      data.SWR.evIndex  = 1:length(data.SWR.evStart);
      data.SWR.evIndex  = data.SWR.evIndex';
      data.SWR.evIndex  = horzcat(data.SWR.evIndex, zeros(length(data.SWR.evStart), 1));
      
    elseif param.swrType == 3 % Ripple only
      data.SWR.evStatus = data.R.evStatus;
      data.SWR.evStart  = data.R.evStart;
      data.SWR.evEnd    = data.R.evEnd;
      data.SWR.evIndex  = 1:length(data.SWR.evStart);
      data.SWR.evIndex  = data.SWR.evIndex';
      data.SWR.evIndex  = horzcat(zeros(length(data.SWR.evStart), 1), data.SWR.evIndex);
      
    end
    
    if ~isnan(data.SWR.evStart)
      
      % Initialize event locked data window cell arrays
      data.SWR.event{length(data.SWR.evStart)}    = [];
      if isfield(data,'SW');     data.SW.SWR.event{length(data.SWR.evStart)}     = []; end
      if isfield(data,'R');      data.R.SWR.event{length(data.SWR.evStart)}      = []; end
      if isfield(data,'gamma');  data.gamma.SWR.event{length(data.SWR.evStart)}  = []; end
      if isfield(data,'hgamma'); data.hgamma.SWR.event{length(data.SWR.evStart)} = []; end
      if isfield(data,'fR');     data.fR.SWR.event{length(data.SWR.evStart)}     = []; end
      
      % Determine baseline for amplitude determination
      if param.swrType == 1 || param.swrType == 2 % use SW signal
        baseAmp = data.SW.tSeries;
        baseAmp(baseAmp > quantile(baseAmp, param.baseQuant)) = [];
        baseAmp = mean(baseAmp);
      elseif param.swrType == 3 % Use LFP signal, SW filter may not have been performed
        baseAmp = data.LFP.tSeries;
        baseAmp(baseAmp > quantile(baseAmp, param.baseQuant)) = [];
        baseAmp = mean(baseAmp);
      end
      
      for i = 1:length(data.SWR.evStart)
        data.SWR.power(i)    = bandpower(data.LFP.tSeries(data.SWR.evStart(i) : data.SWR.evEnd(i)));
        data.SWR.duration(i) = (data.LFP.timing(data.SWR.evEnd(i)) - data.LFP.timing(data.SWR.evStart(i)));
        
        % Peak determination:
        if param.swrType == 1 || param.swrType == 2 % use SW-RMS peak
          data.SWR.evPeak(i) = data.SW.evPeak(find((data.SW.evStart >= data.SWR.evStart(i)) .* (data.SW.evEnd <= data.SWR.evEnd(i)),1));
          data.SWR.amp(i)    = data.SW.tSeries(data.SWR.evPeak(i)) - baseAmp;
        elseif param.swrType == 3 % use R-RMS peak
          data.SWR.evPeak(i) = data.R.evPeak(i);
          data.SWR.amp(i)    = data.LFP.tSeries(data.SWR.evPeak(i)) - baseAmp; % Use LFP signal, SW filter may not have been performed
        end

        if (i > 1); data.SWR.IEI = horzcat(data.SWR.IEI, (data.LFP.timing(data.SWR.evPeak(i)) - data.LFP.timing(data.SWR.evPeak(i-1))) / 1000); end
        
        % Calculate SWR-locked event data
        loWin = max(round(data.SWR.evPeak(i) - param.swrWindow / data.LFP.samplingInt), 1);
        hiWin = min(round(data.SWR.evPeak(i) + param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
        loBaseWin = max(round(data.SWR.evPeak(i) - 0.5 * param.swrWindow / data.LFP.samplingInt), 1);
        hiBaseWin = min(round(data.SWR.evPeak(i) + 0.5 * param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
        
        % SWR data:
        data.SWR.event{i}  = data.LFP.tSeries(loWin : hiWin);
        data.SWR.evTiming  = -param.swrWindow : data.LFP.samplingInt : param.swrWindow;
        data.SWR.area(i)   = data.LFP.samplingInt * sum(sum(data.LFP.tSeries(loBaseWin : hiBaseWin)));
        
        % SW data:
        if isfield(data,'SW')
          data.SW.SWR.event{i} = data.SW.tSeries(loWin : hiWin);
          data.SW.SWR.power(i) = bandpower(data.SW.tSeries(loBaseWin : hiBaseWin));
          data.SW.SWR.area(i)  = data.LFP.samplingInt * sum(sum(data.SW.tSeries(loBaseWin : hiBaseWin)));
        end
        
        % Ripple data:
        if isfield(data,'R')
          data.R.SWR.event{i}  = data.R.tSeries(loWin : hiWin);
          data.R.SWR.power(i)  = bandpower(data.R.tSeries(loBaseWin : hiBaseWin));
        end
        
        % Gamma data:
        if isfield(data,'gamma')
          data.gamma.SWR.event{i} = data.gamma.tSeries(loWin : hiWin);
          data.gamma.SWR.power(i) = bandpower(data.gamma.tSeries(loBaseWin : hiBaseWin));
        end

        % High gamma data:
        if isfield(data,'hgamma')
          data.hgamma.SWR.event{i} = data.hgamma.tSeries(loWin : hiWin);
          data.hgamma.SWR.power(i) = bandpower(data.hgamma.tSeries(loBaseWin : hiBaseWin));
        end
        
        % Fast ripple data:
        if isfield(data,'fR')
          data.fR.SWR.event{i} = data.fR.tSeries(loWin : hiWin);
          data.fR.SWR.power(i) = bandpower(data.fR.tSeries(loBaseWin : hiBaseWin));
        end
        
      end
      
      data.SWR.frequency = length(data.SWR.evStart) / ((data.LFP.timing(length(data.LFP.timing)) - data.LFP.timing(1)) / 1000);
      
      % Standardize array transposition
      data.SWR.power    = data.SWR.power';
      data.SWR.duration = data.SWR.duration';
      data.SWR.evPeak   = data.SWR.evPeak';
      data.SWR.amp      = data.SWR.amp';
      data.SWR.IEI      = data.SWR.IEI';
      data.SWR.event    = data.SWR.event';
      data.SWR.evTiming = data.SWR.evTiming';
      data.SWR.area     = data.SWR.area';
      
      if isfield(data,'SW')
        data.SW.SWR.event = data.SW.SWR.event';
        data.SW.SWR.power = data.SW.SWR.power';
        data.SW.SWR.area  = data.SW.SWR.area';
      end
      
      if isfield(data,'R')
        data.R.SWR.event  = data.R.SWR.event';
        data.R.SWR.power  = data.R.SWR.power';
      end
      
      if isfield(data,'gamma')
        data.gamma.SWR.event = data.gamma.SWR.event';
        data.gamma.SWR.power = data.gamma.SWR.power';
      end

      if isfield(data,'hgamma')
        data.hgamma.SWR.event = data.hgamma.SWR.event';
        data.hgamma.SWR.power = data.hgamma.SWR.power';
      end
      
      if isfield(data,'fR')
        data.fR.SWR.event = data.fR.SWR.event';
        data.fR.SWR.power = data.fR.SWR.power';
      end
      
    end
    fprintf('done\n');
  end
end

%% Spectral Analysis

% Calculate Spectrograms:
if param.spectOption
  fprintf(['calculating spectrogram of total LFP signal (file ' dataFileName ')... ']);
  fRange = param.spectLim1 : param.spectLim2;
  [data.LFP, ~] = calcSpect(data.LFP, [], fRange, round(param.Fs), 30, 0);
  fprintf('done\n');
      
  % If SWR events analyzed, detect spectrogram for event-locked data
  if (param.swrOption)
    fprintf(['calculating spectrograms of SWR-locked events (file ' dataFileName ')... ']);
    [data.SWR, ~] = calcSpect(data.SWR, [], fRange, round(param.Fs), 3, 0);
    fprintf('done\n');
  end
end

% Calculate FFTs:
if param.fftOption
  fprintf(['calculating FFT of total filtered signals (file ' dataFileName ')... ']);
  data.LFP = calcTotFFT(data.LFP, data.param);

  % Compute total FFTs for any other selected bandwidths
  if isfield(data,'theta');  data.theta  = calcTotFFT(data.theta, data.param);  end
  if isfield(data,'alpha');  data.alpha  = calcTotFFT(data.alpha, data.param);  end
  if isfield(data,'beta');   data.beta   = calcTotFFT(data.beta, data.param);   end
  if isfield(data,'SW');     data.SW     = calcTotFFT(data.SW, data.param);     end
  if isfield(data,'gamma');  data.gamma  = calcTotFFT(data.gamma, data.param);  end
  if isfield(data,'hgamma'); data.hgamma = calcTotFFT(data.hgamma, data.param); end
  if isfield(data,'R');      data.R      = calcTotFFT(data.R, data.param);      end
  if isfield(data,'fR');     data.fR     = calcTotFFT(data.fR, data.param);     end
  
  fprintf('done\n');
      
  % If SWR events analyzed, calculate FFT for event-locked data
  if (param.swrOption)
    fprintf(['calculating FFTs of SWR-locked events (file ' dataFileName ')... ']);
    data.SWR = calcEvFFT(data.SWR, data.param, data.param.spectLim1, data.param.spectLim2);
    
    % Compute SWR event-locked FFT and for gamma, high gamma, ripple, and fast ripple (if selected)
    if isfield(data,'gamma');  data.gamma.SWR  = calcEvFFT(data.gamma.SWR, data.param, data.gamma.lim1, data.gamma.lim2);    end
    if isfield(data,'hgamma'); data.hgamma.SWR = calcEvFFT(data.hgamma.SWR, data.param, data.hgamma.lim1, data.hgamma.lim2); end
    if isfield(data,'R');      data.R.SWR      = calcEvFFT(data.R.SWR, data.param, data.R.lim1, data.R.lim2);                end
    if isfield(data,'fR');     data.fR.SWR     = calcEvFFT(data.fR.SWR, data.param, data.fR.lim1, data.fR.lim2);             end
    
    fprintf('done\n');
    
  end
end

% Calculate Phase:
if param.phaseOption
  % Compute interpolated piecewise-linear phase of of total filtered signals
  fprintf(['calculating interpolated piecewise-linear phase of total filtered signals (file ' dataFileName ')... ']);
  if isfield(data,'theta');  data.theta  = calcTotPhase(data.theta, data.LFP, data.param);  end
  if isfield(data,'alpha');  data.alpha  = calcTotPhase(data.alpha, data.LFP, data.param);  end
  if isfield(data,'beta');   data.beta   = calcTotPhase(data.beta, data.LFP, data.param);   end
  if isfield(data,'SW');     data.SW     = calcTotPhase(data.SW, data.LFP, data.param);     end
  if isfield(data,'gamma');  data.gamma  = calcTotPhase(data.gamma, data.LFP, data.param);  end
  if isfield(data,'hgamma'); data.hgamma = calcTotPhase(data.hgamma, data.LFP, data.param); end
  if isfield(data,'R');      data.R      = calcTotPhase(data.R, data.LFP, data.param);      end
  if isfield(data,'fR');     data.fR     = calcTotPhase(data.fR, data.LFP, data.param);     end
  fprintf('done\n');
      
  % If SWR events analyzed, interpolated piecewise-linear phase of of SWR-locked signals
  if (param.swrOption)
    fprintf(['calculating interpolated piecewise-linear phase of SWR-locked signals (file ' dataFileName ')... ']);
    if isfield(data,'gamma');  data.gamma.SWR  = calcEvPhase(data.gamma.SWR, data.SWR, data.param, data.gamma.lim1, data.gamma.lim2);    end
    if isfield(data,'hgamma'); data.hgamma.SWR = calcEvPhase(data.hgamma.SWR, data.SWR, data.param, data.hgamma.lim1, data.hgamma.lim2); end
    if isfield(data,'R');      data.R.SWR      = calcEvPhase(data.R.SWR, data.SWR, data.param, data.R.lim1, data.R.lim2);                end
    if isfield(data,'fR');     data.fR.SWR     = calcEvPhase(data.fR.SWR, data.SWR, data.param, data.fR.lim1, data.fR.lim2);             end
    fprintf('done\n');
  end
end


%% Cross-Frequency Phase-Amplitude Coupling (PAC)
%  Adapted from Canolty et al 2006 and Onslow et al 2011
if param.xFreqOption
  % Assign a temp LFP tSeries (it will be trimmed later)
  tSeries  = data.LFP.tSeries;
  nSample  = length(data.LFP.tSeries);

  %% Total PAC Analysis for n x n matrix:
  if isfield(param, 'spectLim1') && isfield(param, 'spectLim2')
    
    fprintf(['calculating n x n phase-amplitude coupling (file ' dataFileName ')... ']);
    
    % Initialize data structures:
    if ~isfield(data.LFP,'xFreq'); data.LFP.xFreq = struct; end
    
    % Determine frequency vectors:
    nDig1 = numel(num2str(param.spectLim1));
    nDig2 = numel(num2str(param.spectLim2));
    data.LFP.xFreq.freqRange = round5sd(param.spectLim1, nDig1) : param.xFreqBin : round5sd(param.spectLim2, nDig2); % Array in units of param.xFreqBin
    data.LFP.xFreq.morlFreq  = data.LFP.xFreq.freqRange(1 : end-1) + (0.5 * param.xFreqBin); % Final frequency values are mid-points
    nFreq = length(data.LFP.xFreq.morlFreq);
    
    % Initialize phase, amplitude, and PAC arrays
    data.LFP.xFreq.phsPAC = zeros(nSample, nFreq);
    data.LFP.xFreq.ampPAC = zeros(nSample, nFreq);
    data.LFP.xFreq.pacMI  = zeros(nFreq, nFreq);
    
    % Determine phase and amplitude via Morlet wavelet:
    for i = 1 : nFreq
      data.LFP.xFreq.phsPAC(:,i) = morletPhase(data.LFP.xFreq.morlFreq(i), tSeries, round(param.Fs), param.morlWidth);
      data.LFP.xFreq.ampPAC(:,i) = morletAmp(data.LFP.xFreq.morlFreq(i), tSeries, round(param.Fs), param.morlWidth);
    end
    
    % Calculate PAC modulation index (MI):
    data.LFP.xFreq = calcPACMI(data.LFP.xFreq, round(param.Fs), 0); % nShuffle > 1 gets very computationally expensive!
    
    % Order structure:
    data.LFP.xFreq = orderStruct(data.LFP.xFreq);
    
    fprintf('done\n');
  end
  
  %% Total PAC Analysis for each higher frequency band:
  % Calculate modulating phase (lower frequency) via Morlet wavelet:
  if isfield(data, param.xFreqLow)
    fprintf(['calculating Z-corrected phase-amplitude coupling for frequency bands of interest (file ' dataFileName ')... ']);
    
    morlFreqP = data.(param.xFreqLow).lim1 + floor((data.(param.xFreqLow).lim2 - data.(param.xFreqLow).lim1)/2);
    phsPAC    = morletPhase(morlFreqP, tSeries, round(param.Fs), param.morlWidth);
    
    % Determine available higher phase-modulated frequency(ies) available:
    nHi = 0;
    if isfield(data,'gamma')
      nHi = nHi + 1;
      xFreqHi{nHi} = 'gamma';
    end
    if isfield(data,'hgamma')
      nHi = nHi + 1;
      xFreqHi{nHi} = 'hgamma';
    end
    if isfield(data,'R')
      nHi = nHi + 1;
      xFreqHi{nHi} = 'R';
    end
    if isfield(data,'fR')
      nHi = nHi + 1;
      xFreqHi{nHi} = 'fR';
    end
    
    if nHi > 0
      for i = 1:nHi
        if ~isfield(data.(xFreqHi{i}),'xFreq'); data.(xFreqHi{i}).xFreq = struct; end
        
        % Calculate amplitude of higher frequency(ies) via Morlet wavelet:
        morlFreqA = data.(xFreqHi{i}).lim1 + floor((data.(xFreqHi{i}).lim2 - data.(xFreqHi{i}).lim1)/2);
        ampPAC    = morletAmp(morlFreqA, tSeries, round(param.Fs), param.morlWidth);
        
        % Update data structure:
        data.(xFreqHi{i}).xFreq.xFreqLow  = param.xFreqLow;
        data.(xFreqHi{i}).xFreq.xFreqHi   = xFreqHi{i};
        data.(xFreqHi{i}).xFreq.morlFreqP = morlFreqP;
        data.(xFreqHi{i}).xFreq.morlFreqA = morlFreqA;
        data.(xFreqHi{i}).xFreq.phsPAC    = phsPAC;
        data.(xFreqHi{i}).xFreq.ampPAC    = ampPAC;

        % Calculate total PAC measure:
        data.(xFreqHi{i}).xFreq = calcPACMI(data.(xFreqHi{i}).xFreq, round(param.Fs), param.nShuffle);
      end
      
      %% Time PAC Analysis
      % Truncate signals to get integer number of time windows
      nSampWn   = ceil(param.winLength * round(param.Fs));
      nSampOl   = ceil(param.winOverlap * round(param.Fs));
      remSample = mod(nSample, nSampWn);
      phsPAC    = phsPAC(1 : nSample - remSample);
      ampPAC    = ampPAC(1 : nSample - remSample, :);
      
      % Update nSample
      nSample = length(phsPAC);
      idx     = bsxfun(@plus, (1:nSampWn)', 1+(0:(fix((nSample - nSampOl)/(nSampWn - nSampOl)) - 1))*(nSampWn - nSampOl)) - 1;
      nWin    = size(idx,2);
      
      % Determine average time of windows:
      timingWin = zeros(nWin, 1);
      for k = 1:size(idx, 2)
        timingWin(k) = mean(data.LFP.timing(idx(:, k)));
      end
      
      % Calculate windowed time-series PAC:
      pacMIWin       = zeros(nWin, 1);
      pacMIWin_Len   = zeros(nWin, 1);
      pacMIWin_Phase = zeros(nWin, 1);
      
      for i = 1:nHi
        for k = 1:size(idx, 2)
          z = ampPAC(idx(:,k), i) .* exp(1i * phsPAC(idx(:,k))); % Create composite signal
          pacMIWin(k)       = mean(z);  % Compute the mean length of composite signal
          pacMIWin_Len(k)   = abs(pacMIWin(k));
          pacMIWin_Phase(k) = angle(pacMIWin(k));
        end
        
        % Update data structure:
        data.(xFreqHi{i}).xFreq.pacMIWin        = pacMIWin;
        data.(xFreqHi{i}).xFreq.pacMIWin_Len    = pacMIWin_Len;
        data.(xFreqHi{i}).xFreq.pacMIWin_Phase  = pacMIWin_Phase;
        data.(xFreqHi{i}).xFreq.timingWin       = timingWin;
        
        % Linear interpolation of pacMIWin_Len to data samplingInt for later correlations
        data.(xFreqHi{i}).xFreq.timingI = data.LFP.timing(1:nSample);
        data.(xFreqHi{i}).xFreq.pacMIWin_LenI = interp1(timingWin, pacMIWin_Len, data.(xFreqHi{i}).xFreq.timingI);
        
        % Order structure:
        data.(xFreqHi{i}).xFreq = orderStruct(data.(xFreqHi{i}).xFreq);
      end
    end
  end
  fprintf('done\n');
end


%% Import and process stim file (if selected)
if param.importStimOption
  
  % Import stim event file into matlab table
  warning('off')
  fprintf(['importing stim event file ' stimFileName '... ']);
  stimTable = readtable(stimFile);
  fprintf('done\n');
  warning('on')
  
  data.stim = struct;
  
  % Re-initialize spike data structures
  data.stim.evStatus = [];
  data.stim.evStart  = [];
  data.stim.evPeak   = [];
  data.stim.evEnd    = [];
  
  fprintf(['processing stim events (file ' dataFileName ')... ']);
  for i = 1:size(stimTable,2)
    if ~isempty(strfind(stimTable.Properties.VariableNames{i},"stimStart"))
      stimStartTime = stimTable{:,i};
      for ev = 1:length(stimStartTime)
        if ~isempty(find(data.LFP.timing >= stimStartTime(ev), 1))
          data.stim.evStart(ev) = find(data.LFP.timing >= stimStartTime(ev), 1);
        end
      end
      data.stim.evStart = data.stim.evStart';
    end
    
    if ~isempty(strfind(stimTable.Properties.VariableNames{i},"stimEnd"))
      stimEndTime = stimTable{:,i};
      for ev = 1:length(stimEndTime)
        if ~isempty(find(data.LFP.timing <= stimEndTime(ev), 1, 'last'))
          data.stim.evEnd(ev) = find(data.LFP.timing <= stimEndTime(ev), 1, 'last');
        end
      end
      data.stim.evEnd = data.stim.evEnd';
    end
    
    if ~isempty(strfind(stimTable.Properties.VariableNames{i},"stimPeak"))
      stimPeakTime = stimTable{:,i};
      for ev = 1:length(stimPeakTime)
        if ~isempty(find(data.LFP.timing >= stimPeakTime(ev), 1, 'last'))
          data.stim.evPeak(ev) = find(data.LFP.timing >= stimPeakTime(ev), 1);
        end
      end
      data.stim.evPeak = data.stim.evPeak';
    end
    
    if ~isempty(strfind(stimTable.Properties.VariableNames{i},"lfpAmp"))
      data.stim.lfpAmp = stimTable{:,i};
    end
    
    if ~isempty(strfind(stimTable.Properties.VariableNames{i},"lfpBL"))
      data.stim.lfpBL = stimTable{:,i};
    end
    
    if ~isempty(strfind(stimTable.Properties.VariableNames{i},"lfpSlope"))
      data.stim.lfpSlope = stimTable{:,i};
    end
    
  end
  fprintf('done\n');
  
  % compute stim status array
  data.stim.evStatus = zeros(length(data.LFP.timing),1);
  for ev = 1:length(data.stim.evStart)
    for i = data.stim.evStart(ev) : data.stim.evEnd(ev)
      data.stim.evStatus(i) = 1;
    end
  end
end

%% Re-order structure arrays
data       = orderStruct(data);
data.param = orderStruct(data.param);
data.LFP   = orderStruct(data.LFP);
data.LFP.param = orderStruct(data.LFP.param);

if isfield(data, 'SW')
  data.SW = orderStruct(data.SW);
  if isfield(data.SW, 'SWR'); data.SW.SWR = orderStruct(data.SW.SWR); end
end

if isfield(data, 'R')
  data.R = orderStruct(data.R);
  if isfield(data.R, 'SWR'); data.R.SWR = orderStruct(data.R.SWR); end
end

if isfield(data, 'SWR');    data.SWR    = orderStruct(data.SWR); end
if isfield(data, 'theta');  data.theta  = orderStruct(data.theta); end
if isfield(data, 'alpha');  data.alpha  = orderStruct(data.alpha); end
if isfield(data, 'beta');   data.beta   = orderStruct(data.beta); end

if isfield(data, 'gamma')
  data.gamma = orderStruct(data.gamma);
  if isfield(data.gamma, 'SWR'); data.gamma.SWR = orderStruct(data.gamma.SWR); end
end

if isfield(data, 'hgamma')
  data.hgamma = orderStruct(data.hgamma);
  if isfield(data.hgamma, 'SWR'); data.hgamma.SWR = orderStruct(data.hgamma.SWR); end
end

if isfield(data, 'fR')
  data.fR = orderStruct(data.fR);
  if isfield(data.fR, 'SWR'); data.fR.SWR = orderStruct(data.fR.SWR); end
end

if isfield(data, 'C')
  data.C = orderStruct(data.C);
  if isfield(data.C, 'SWR'); data.C.SWR = orderStruct(data.C.SWR); end
end

if isfield(data, 'stim')
  data.stim = orderStruct(data.stim);
end

%% Save and export results
% Export SWR event table
if all(expEvFile) && param.expSWREvOption
  fprintf(['exporting SWR events (file ' dataFileName ')... ']);
  exportSWREvents(data, saveFile, expEvFile);
  data.SWR.expEvFile = expEvFile;
  fprintf('done\n');
end

% Export SWR event-locked episodic data files
if all(expDataFile) && param.expSWRDataOption
  fprintf(['exporting SWR event-locked data (file ' dataFileName ')... ']);
  exportSWRData(data, param, expDataFile);
  data.SWR.expDataFile = expDataFile;
  fprintf('done\n');
end

% Export average statistics
if all(expAveFile) && param.expAveOption
  fprintf(['calculating average statistics (file ' dataFileName ')... ']);
  [aveStats, varNames] = calcAveStats(data, param);
  data.LFP.expAveFile = expAveFile;
  if (param.fileNum == 1); writetable(aveStats, expAveFile, 'Delimiter', ',', 'WriteVariableNames', true, 'WriteRowNames', true); end
  fprintf('done\n');
else
  aveStats{1} = [];
  varNames{1} = [];
end

% Save matlab file
if all(saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(saveFile,'-struct','data');
  fprintf('done\n');
end

if (param.fileType == 2) && ~param.reAnalyzeOption
  cd (curPath);
end

end
