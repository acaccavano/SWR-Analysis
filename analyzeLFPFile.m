function [data, hand] = analyzeLFPFile(data, hand, param, dataFile, saveFile, expEvFile, expDataFile)
%% [data, h] = analyzeLFPFile(data, hand, param, dataFile, saveFile, expEvFile, expDataFile)
%
%  Function to detect sharp wave ripple (SWR) events, gamma and theta analysis,
%  and spectrogram analysis of single LFP recording. SWRs are detected by
%  performing both a low frequency (sharp wave) and high frequency (ripple) band-pass
%  filter of an LFP recording. The root-mean-square (RMS) of these filtered
%  channels are taken, and events are counted if they exceed a given
%  standard deviation of the RMS channels. SWR events are counted if both
%  a sharp wave and ripple occur simultaneously.
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   data       = structure - specify if appending previously analyzed files
%   hand       = handle structure to specify where figure should be drawn
%   param      = structure containing all parameters including:
%     param.fileNum          = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.fileType         = 1 = pClamp (.abf), 2 = Wu data (folder of data files)
%     param.Fs               = sampling rate (JYs recordings are usually 3000, not needed for pClamp files)
%     param.dsFactor         = downsample factor (default = 1, no downsampling)
%     param.lfpChannel       = channel to use for LFP input (default = 1, but depends on recording)
%     param.cellOption       = boolean flag to determine if second cell channel to be imported
%     param.cellChannel      = channel to use for optional cell input (default = 2, but depends on recording)
%     param.lfpOption        = boolean flag to filter LFP signal
%     param.lfpLim1          = lower LFP band-pass lim (default = 1Hz)
%     param.lfpLim2          = upper LFP band-pass lim (default = 1000Hz)
%     param.swrOption        = boolean flag to detect SWR events
%     param.swOption         = boolean flag to filter and analyze SW signal
%     param.swLim1           = lower sharp wave band-pass lim (default = 1Hz)
%     param.swLim2           = upper sharp wave band-pass lim (default = 30Hz)
%     param.rOption          = boolean flag to filter and analyze ripple signal
%     param.rLim1            = lower ripple band-pass lim (default = 120)
%     param.rLim2            = upper ripple band-pass lim (default = 220)
%     param.rmsOption        = boolean flag to calculate RMS of SW and ripple
%     param.rmsPeriod        = root-mean square window [ms] (in Eschenko 2006, rmsPeriod = 5ms)
%     param.peakDetectOption = boolean flag to detect SW and ripple reaks in RMS signals
%     param.sdMult           = standard deviation threshold to detect peaks (default = 4)
%     param.baseQuant        = quantile with which to determine baseline signal (default = 0.95)
%     param.swrWindow        = +/- window around SWR peak events for swrData file [ms]
%     param.expSWREvOption   = boolean flag to determine whether to export txt table of SWR events
%     param.expSWRDataOption = boolean flag to determine whether to export txt file of episodic SWR events for pClamp analysis
%     param.thetaOption      = boolean flag to filter and analyze theta signal
%     param.thetaLim1        = lower theta band-pass lim (default = 4Hz)
%     param.thetaLim2        = upper theta band-pass lim (default = 12Hz)
%     param.betaOption       = boolean flag to filter and analyze beta signal
%     param.betaLim1         = lower beta band-pass lim (default = 13Hz)
%     param.betaLim2         = upper beta band-pass lim (default = 24Hz)
%     param.gammaOption      = boolean flag to filter and analyze gamma signal
%     param.gammaLim1        = lower gamma band-pass lim (default = 25Hz)
%     param.gammaLim2        = upper gamma band-pass lim (default = 55Hz)
%     param.hgammaOption     = boolean flag to filter and analyze high gamma signal
%     param.hgammaLim1       = lower high gamma band-pass lim (default = 65Hz)
%     param.hgammaLim2       = upper high gamma band-pass lim (default = 85Hz)
%     param.spectOption      = boolean flag to perform spectrogram analysis
%     param.spectLim1        = lower lim of spectrogram (default = 1Hz)
%     param.spectLim2        = upper lim of spectrogram (default = 250Hz)
%   dataFile    = full path to file/folder containing data to be analysed (if not set, will prompt)
%   saveFile    = full path to matlab file to save (if not set, will prompt)
%   expEvFile   = full path to exported txt event table (if not set and expSWREvOption = 1, will prompt
%   expDataFile = full path to exported txt data file (if not set and expSWRDataOption = 1, will prompt
%
%  Outputs:
%   data       = structure containing all data to be saved
%   hand       = handle structure for figure

%% Handle input arguments - if not entered
if (nargin < 7) expDataFile = []; end
if (nargin < 6) expEvFile   = []; end
if (nargin < 5) saveFile    = []; end
if (nargin < 4) dataFile    = []; end
if (nargin < 3) param       = struct; end
if (nargin < 2) hand        = struct; end
if (nargin < 1) data        = struct; end

% Handle case in which empty variables are supplied:
if isempty(param) param     = struct; end
if isempty(hand)  hand      = struct; end
if isempty(data)  data      = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum')          param.fileNum           = 1; end
if ~isfield(param,'fileType')         param.fileType          = 2;    end
if ~isfield(param,'Fs')               param.Fs                = 3000; end
if ~isfield(param,'dsFactor')         param.dsFactor          = 1;    end
if ~isfield(param,'lfpChannel')       param.lfpChannel        = 1;    end
if ~isfield(param,'cellOption')       param.cellOption        = 1;    end
if ~isfield(param,'cellChannel')      param.cellChannel       = 2;    end
if ~isfield(param,'notchOption')      param.notchOption       = 0;    end
if ~isfield(param,'lfpOption')        param.lfpOption         = 1;    end
if ~isfield(param,'lfpLim1')          param.lfpLim1           = 1;    end
if ~isfield(param,'lfpLim2')          param.lfpLim2           = 1000; end
if ~isfield(param,'swrOption')        param.swrOption         = 1;    end
if ~isfield(param,'swOption')         param.swOption          = 1;    end
if ~isfield(param,'swLim1')           param.swLim1            = 1;    end
if ~isfield(param,'swLim2')           param.swLim2            = 30;   end
if ~isfield(param,'rOption')          param.rOption           = 1;    end
if ~isfield(param,'rLim1')            param.rLim1             = 120;  end
if ~isfield(param,'rLim2')            param.rLim2             = 220;  end
if ~isfield(param,'rmsOption')        param.rmsOption         = 1;    end
if ~isfield(param,'rmsPeriod')        param.rmsPeriod         = 5;    end
if ~isfield(param,'peakDetectOption') param.peakDetectOption  = 1;    end
if ~isfield(param,'sdMult')           param.sdMult            = 3;    end
if ~isfield(param,'baseQuant')        param.baseQuant         = 0.95; end
if ~isfield(param,'swrWindow')        param.swrWindow         = 100;  end
if ~isfield(param,'expSWREvOption')   param.expSWREvOption    = 1;    end
if ~isfield(param,'expSWRDataOption') param.expSWRDataOption  = 1;    end
if ~isfield(param,'thetaOption')      param.thetaOption       = 0;    end
if ~isfield(param,'thetaLim1')        param.thetaLim1         = 4;    end
if ~isfield(param,'thetaLim2')        param.thetaLim2         = 12;   end
if ~isfield(param,'betaOption')       param.betaOption        = 0;    end
if ~isfield(param,'betaLim1')         param.betaLim1          = 13;   end
if ~isfield(param,'betaLim2')         param.betaLim2          = 24;   end
if ~isfield(param,'gammaOption')      param.gammaOption       = 1;    end
if ~isfield(param,'gammaLim1')        param.gammaLim1         = 20;   end
if ~isfield(param,'gammaLim2')        param.gammaLim2         = 50;   end
if ~isfield(param,'hgammaOption')     param.hgammaOption      = 0;    end
if ~isfield(param,'hgammaLim1')       param.hgammaLim1        = 65;   end
if ~isfield(param,'hgammaLim2')       param.hgammaLim2        = 85;   end
if ~isfield(param,'fROption')         param.fROption          = 1;    end
if ~isfield(param,'fRLim1')           param.fRLim1            = 250;  end
if ~isfield(param,'fRLim2')           param.fRLim2            = 600;  end
if ~isfield(param,'spectOption')      param.spectOption       = 1;    end
if ~isfield(param,'spectLim1')        param.spectLim1         = 1;    end
if ~isfield(param,'spectLim2')        param.spectLim2         = 600;  end
if ~isfield(param,'reAnalyzeOption')  param.reAnalyzeOption   = 0;    end

% Initialize LFP structure if it doesn't already exist
if ~isfield(data,'LFP') data.LFP = struct; end

% If cell option is selected, initialize cell structure if it doesn't already exist
if param.cellOption
  if ~isfield(data,'C') data.C = struct; end
end

% If not supplied, prompt for files/folder to analyze
if isfield(data.LFP, 'dataFile')
  dataFile = data.LFP.dataFile;
elseif isempty(dataFile)
  if (param.fileType == 1)
    [fileName, filePath] = uigetfile('.abf', 'Select *.abf file to analyze');
    dataFile = strcat(filePath,fileName);
  elseif (param.fileType == 2)
    dataFile = uigetdir();
  end
  if ~all(dataFile) return; end
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
  if ~all(expDataFile) warning('No SWR events to be exported - no file selected'); end
end

%% Import data
if ~isfield(data.LFP, 'dataFile')
  fprintf(['importing file ' dataFileName '... ']);
  if (param.fileType == 1)
    [dataIn, samplingInt] = abfload(dataFile);
    data.LFP.samplingInt = samplingInt;
    data.LFP.tSeries = dataIn(:, param.lfpChannel);
    data.LFP.samplingInt = data.LFP.samplingInt / 1000; % convert from um to ms
    
    if param.cellOption
      data.C.samplingInt = samplingInt;
      data.C.tSeries = dataIn(:, param.cellChannel);
      data.C.samplingInt = data.C.samplingInt / 1000; % convert from um to ms
    end
    
  elseif (param.fileType == 2)
    % Ensure current directory is in path so helper functions work
    curPath = pwd;
    path(path, curPath);
    
    % Extract file names
    cd (dataFile);
    dir_temp = dir('2*');
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
        str_ascii= str2num(int2str(string));
        index_nums = find(str_ascii <=57 & str_ascii >=48);
        filenum(i) = str2num(string(index_nums));
      end
      [filenum, index]=sort(filenum,'ascend');
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
      dataTemp = [data_firstline dataTemp];
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
    data.LFP.tSeries = dataIn(:,1);
    data.LFP.samplingInt = 1000 / param.Fs; % (ms)
    if param.cellOption
      data.C.tSeries = dataIn(:,2);
      data.C.samplingInt = 1000 / param.Fs; % (ms)
    end
  end
  fprintf('done\n');
  
  %% Assign parameters to structure array
  data.LFP.dataFile  = dataFile;
  
  % Downsample data if selected
  if (param.dsFactor >= 2)
    fprintf(['downsampling by factor of ' num2str(param.dsFactor) ' (file ' dataFileName ')... ']);
    data.LFP.samplingInt = data.LFP.samplingInt * param.dsFactor;
    data.LFP.tSeries = downsampleMean(data.LFP.tSeries, param.dsFactor);
    
    if param.cellOption
      data.C.samplingInt = data.C.samplingInt * param.dsFactor;
      data.C.tSeries = downsampleMean(data.C.tSeries, param.dsFactor);
    end
    
    fprintf('done\n');
  end
  
  data.LFP.nSamples = length(data.LFP.tSeries);
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
if param.lfpOption
  % Apply Gaussian filter to LFP signal for DC drift and HF noise
  fprintf(['band-pass filtering LFP between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.lfpLim1, param.lfpLim2);
  data.LFP.tSeries = gaussianFilt(data.LFP.tSeries, param.lfpLim1, param.lfpLim2, data.LFP.samplingInt, 5);
  data.LFP.tPower  = bandpower(data.LFP.tSeries);
  data.LFP.lim1    = param.lfpLim1;
  data.LFP.lim2    = param.lfpLim2;
  
  if param.notchOption
    fprintf(['Notch filter 60Hz noise (file ' dataFileName ')... ']);
    
    Ord   = 50;  % Order
    BW    = 10;  % Bandwidth
    Apass = 1;   % Bandwidth Attenuation
    
    [b, a] = iircomb(Ord, BW/(data.param.Fs/2), Apass);
    data.LFP.tSeries = filtfilt(b, a, data.LFP.tSeries);
    
  end
  
  fprintf('done\n');
end

if param.swOption
  % Apply Gaussian filter to extract SW signal
  fprintf(['band-pass filtering sharp wave between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.swLim1, param.swLim2);
  if ~isfield(data,'SW') data.SW = struct; end
  data.SW.tSeries = gaussianFilt(data.LFP.tSeries, param.swLim1, param.swLim2, data.LFP.samplingInt, 3);
  data.SW.tPower  = bandpower(data.SW.tSeries);
  data.SW.lim1    = param.swLim1;
  data.SW.lim2    = param.swLim2;
  fprintf('done\n');
end

if param.rOption
  % Apply Gaussian filter to extract ripple signal
  fprintf(['band-pass filtering ripple between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.rLim1, param.rLim2);
  if ~isfield(data,'R') data.R = struct; end
  data.R.tSeries = gaussianFilt(data.LFP.tSeries, param.rLim1, param.rLim2, data.LFP.samplingInt, 1);
  data.R.tPower  = bandpower(data.R.tSeries);
  data.R.lim1    = param.rLim1;
  data.R.lim2    = param.rLim2;
  fprintf('done\n');
end

if param.thetaOption
  % Apply Gaussian filter to extract theta signal
  fprintf(['band-pass filtering theta between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.thetaLim1, param.thetaLim2);
  if ~isfield(data,'theta') data.theta = struct; end
  data.theta.tSeries = gaussianFilt(data.LFP.tSeries, param.thetaLim1, param.thetaLim2, data.LFP.samplingInt, 2);
  data.theta.tPower  = bandpower(data.theta.tSeries);
  data.theta.lim1    = param.thetaLim1;
  data.theta.lim2    = param.thetaLim2;
  fprintf('done\n');
end

if param.betaOption
  % Apply Gaussian filter to extract beta signal
  fprintf(['band-pass filtering beta between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.betaLim1, param.betaLim2);
  if ~isfield(data,'beta') data.beta = struct; end
  data.beta.tSeries = gaussianFilt(data.LFP.tSeries, param.betaLim1, param.betaLim2, data.LFP.samplingInt, 1);
  data.beta.tPower  = bandpower(data.beta.tSeries);
  data.beta.lim1    = param.betaLim1;
  data.beta.lim2    = param.betaLim2;
  fprintf('done\n');
end

if param.gammaOption
  % Apply Gaussian filter to extract gamma signal
  fprintf(['band-pass filtering gamma between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.gammaLim1, param.gammaLim2);
  if ~isfield(data,'gamma') data.gamma = struct; end
  data.gamma.tSeries = gaussianFilt(data.LFP.tSeries, param.gammaLim1, param.gammaLim2, data.LFP.samplingInt, 1);
  data.gamma.tPower  = bandpower(data.gamma.tSeries);
  data.gamma.lim1    = param.gammaLim1;
  data.gamma.lim2    = param.gammaLim2;
  fprintf('done\n');
end

if param.hgammaOption
  % Apply Gaussian filter to extract high gamma signal
  fprintf(['band-pass filtering high gamma between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.hgammaLim1, param.hgammaLim2);
  if ~isfield(data,'hgamma') data.hgamma = struct; end
  data.hgamma.tSeries = gaussianFilt(data.LFP.tSeries, param.hgammaLim1, param.hgammaLim2, data.LFP.samplingInt, 1);
  data.hgamma.tPower  = bandpower(data.hgamma.tSeries);
  data.hgamma.lim1    = param.hgammaLim1;
  data.hgamma.lim2    = param.hgammaLim2;
  fprintf('done\n');
end

if param.fROption
  % Apply Gaussian filter to extract fast ripple signal
  fprintf(['band-pass filtering fast ripple between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.fRLim1, param.fRLim2);
  if ~isfield(data,'fR') data.fR = struct; end
  data.fR.tSeries = gaussianFilt(data.LFP.tSeries, param.fRLim1, param.fRLim2, data.LFP.samplingInt, 1);
  data.fR.tPower  = bandpower(data.fR.tSeries);
  data.fR.lim1    = param.fRLim1;
  data.fR.lim2    = param.fRLim2;
  fprintf('done\n');
end


%% SWR event detection if both SW and ripple option enabled
if param.swrOption
  if ~isfield(data,'SWR') data.SWR = struct; end
  
  % RMS Signal calculations for SWR detection
  if param.rmsOption
    fprintf(['calculating root mean square (RMS) in a %4.1fms sliding window (file ' dataFileName ')... '], param.rmsPeriod);
    rmsFactor = round(param.rmsPeriod / data.LFP.samplingInt);
    nSamplesRMS = floor(data.LFP.nSamples / rmsFactor);
    data.SW.RMS = zeros(nSamplesRMS, 1);
    data.R.RMS  = zeros(nSamplesRMS, 1);
    
    for i = 1:nSamplesRMS
      loAve          = max((i - 2) * rmsFactor, 1);
      hiAve          = min(i * rmsFactor, data.LFP.nSamples);
      data.SW.RMS(i) = rms(data.SW.tSeries(loAve:hiAve));
      data.R.RMS(i)  = rms(data.R.tSeries(loAve:hiAve));
    end
    data.param.rmsPeriod = rmsFactor * data.LFP.samplingInt; % Corrected rmsPeriod
    timingRMS        = (0: data.param.rmsPeriod: (nSamplesRMS-1) * data.param.rmsPeriod)';
    data.SW.RMS      = interp1(timingRMS, data.SW.RMS, data.LFP.timing, 'spline');
    data.R.RMS       = interp1(timingRMS, data.R.RMS, data.LFP.timing, 'spline');
    fprintf('done\n');
  end
  
  %% Event Detection
  if param.peakDetectOption
    fprintf(['detecting events %4.0f standard deviations above baseline (%4.2f quantile) (file ' dataFileName ')... '], param.sdMult, param.baseQuant);
    
    %% Find sharp wave peak based on standard deviation of RMS of SW signal
    % Re-initialize data structures
    data.SW.evStatus = [];
    data.SW.evStart  = [];
    data.SW.evPeak   = [];
    data.SW.evEnd    = [];
    data.SW.IEI      = [];
    data.SW.power    = [];
    data.SW.duration = [];
    
    % Calculate baseline by taking bottom baseQuant quantile of signal to minimize false negatives for active fields
    baseline = data.SW.RMS;
    baseline(baseline > quantile(baseline, param.baseQuant)) = [];
    sd = std(baseline);
    mn = mean(baseline);
    data.SW.peakThresh = mn + sd * param.sdMult;
    data.SW.baseThresh = mn + 0.5 * sd * param.sdMult;
    [data.SW.evStatus, data.SW.evStart, data.SW.evPeak, data.SW.evEnd] = peakFindUnique(data.SW.RMS, data.LFP.timing, data.SW.peakThresh, data.SW.baseThresh, 1);
    
    if ~isnan(data.SW.evStart)
      for i = 1:length(data.SW.evStart)
        data.SW.power(i)    = bandpower(data.SW.tSeries(data.SW.evStart(i) : data.SW.evEnd(i)));
        data.SW.duration(i) = (data.LFP.timing(data.SW.evEnd(i)) - data.LFP.timing(data.SW.evStart(i)));
        if (i > 1) data.SW.IEI = horzcat(data.SW.IEI, (data.LFP.timing(data.SW.evPeak(i)) - data.LFP.timing(data.SW.evPeak(i-1))) / 1000); end
      end
      data.SW.frequency = length(data.SW.evStart) / ((data.LFP.timing(length(data.LFP.timing)) - data.LFP.timing(1)) / 1000);
      data.SW.power     = data.SW.power';
      data.SW.duration  = data.SW.duration';
      data.SW.IEI       = data.SW.IEI';
    end
    
    %% Find ripple peak based on standard deviation of RMS of ripple signal
    % Re-initialize data structures
    data.R.evStatus = [];
    data.R.evStart  = [];
    data.R.evPeak   = [];
    data.R.evEnd    = [];
    data.R.IEI      = [];
    data.R.power    = [];
    data.R.duration = [];
    
    % Calculate baseline by taking bottom baseQuant quantile of signal to minimize false negatives for active fields
    baseline = data.R.RMS;
    baseline(baseline > quantile(baseline, param.baseQuant)) = [];
    sd = std(baseline);
    mn = mean(baseline);
    data.R.peakThresh = mn + sd * param.sdMult;
    data.R.baseThresh = mn + 0.5 * sd * param.sdMult;
    [data.R.evStatus, data.R.evStart, data.R.evPeak, data.R.evEnd] = peakFindUnique(data.R.RMS, data.LFP.timing, data.R.peakThresh, data.R.baseThresh, 1);
    
    if ~isnan(data.R.evStart)
      for i = 1:length(data.R.evStart)
        data.R.power(i)    = bandpower(data.R.tSeries(data.R.evStart(i) : data.R.evEnd(i)));
        data.R.duration(i) = (data.LFP.timing(data.R.evEnd(i)) - data.LFP.timing(data.R.evStart(i)));
        if (i > 1) data.R.IEI = horzcat(data.R.IEI, (data.LFP.timing(data.R.evPeak(i)) - data.LFP.timing(data.R.evPeak(i-1))) / 1000); end
      end
      data.R.frequency = length(data.R.evStart) / ((data.LFP.timing(length(data.LFP.timing)) - data.LFP.timing(1)) / 1000);
      data.R.power     = data.R.power';
      data.R.duration  = data.R.duration';
      data.R.IEI       = data.R.IEI';
    end
    
    %% SWR Event Calculation
    
    % (Re)Initialize data arrays
    data.SWR.evStatus = [];
    data.SWR.evStart  = [];
    data.SWR.evPeak   = [];
    data.SWR.evEnd    = [];
    data.SWR.IEI      = [];
    data.SWR.duration = [];
    data.SWR.amp      = [];
    data.SWR.power    = [];
    data.SWR.area     = [];
    data.SWR.event    = [];
    
    % SW arrays:
    if ~isfield(data.SW,'SWR') data.SW.SWR = struct; end
    data.SW.SWR.event = [];
    data.SW.SWR.power = [];
    data.SW.SWR.area  = [];
    
    % Ripple arrays:
    if ~isfield(data.R,'SWR') data.R.SWR = struct; end
    data.R.SWR.event  = [];
    data.R.SWR.power  = [];
    
    % Gamma arrays:
    if isfield(data,'gamma')
      if ~isfield(data.gamma,'SWR') data.gamma.SWR = struct; end
      data.gamma.SWR.event = [];
      data.gamma.SWR.power = [];
    end
    
    % Fast ripple arrays:
    if isfield(data,'fR')
      if ~isfield(data.fR,'SWR') data.fR.SWR = struct; end
      data.fR.SWR.event = [];
      data.fR.SWR.power = [];
    end
    
    [data.SWR.evStatus, data.SWR.evStart, data.SWR.evEnd, data.SWR.evIndex] = eventOverlap(data.SW.evStatus, data.SW.evStart, data.SW.evEnd, data.R.evStatus, data.R.evStart, data.R.evEnd, data.LFP.timing, 0);
    if ~isnan(data.SWR.evStart)
      
      % Initialize event locked data window cell arrays
      data.SWR.event{length(data.SWR.evStart)}    = [];
      data.SW.SWR.event{length(data.SWR.evStart)} = [];
      data.R.SWR.event{length(data.SWR.evStart)}  = [];
      if isfield(data,'gamma') data.gamma.SWR.event{length(data.SWR.evStart)} = []; end
      if isfield(data,'fR')    data.fR.SWR.event{length(data.SWR.evStart)}    = []; end
      
      baseAmp = data.SW.tSeries;
      baseAmp(baseAmp > quantile(baseAmp, param.baseQuant)) = [];
      baseAmp = mean(baseAmp);
      
      for i = 1:length(data.SWR.evStart)
        data.SWR.power(i)    = bandpower(data.LFP.tSeries(data.SWR.evStart(i) : data.SWR.evEnd(i)));
        data.SWR.duration(i) = (data.LFP.timing(data.SWR.evEnd(i)) - data.LFP.timing(data.SWR.evStart(i)));
        data.SWR.evPeak(i)   = data.SW.evPeak(find((data.SW.evStart >= data.SWR.evStart(i)) .* (data.SW.evEnd <= data.SWR.evEnd(i)),1));
        data.SWR.amp(i)      = data.SW.tSeries(data.SWR.evPeak(i)) - baseAmp;
        if (i > 1) data.SWR.IEI = horzcat(data.SWR.IEI, (data.LFP.timing(data.SWR.evPeak(i)) - data.LFP.timing(data.SWR.evPeak(i-1))) / 1000); end
        
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
        data.SW.SWR.event{i} = data.SW.tSeries(loWin : hiWin);
        data.SW.SWR.power(i) = bandpower(data.SW.tSeries(loBaseWin : hiBaseWin));
        data.SW.SWR.area(i)  = data.LFP.samplingInt * sum(sum(data.SW.tSeries(loBaseWin : hiBaseWin)));
        
        % Ripple data:
        data.R.SWR.event{i}  = data.R.tSeries(loWin : hiWin);
        data.R.SWR.power(i)  = bandpower(data.R.tSeries(loBaseWin : hiBaseWin));
        
        % Gamma data:
        if isfield(data,'gamma')
          data.gamma.SWR.event{i} = data.gamma.tSeries(loWin : hiWin);
          data.gamma.SWR.power(i) = bandpower(data.gamma.tSeries(loBaseWin : hiBaseWin));
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
      
      data.SW.SWR.event = data.SW.SWR.event';
      data.SW.SWR.power = data.SW.SWR.power';
      data.SW.SWR.area  = data.SW.SWR.area';
      
      data.R.SWR.event  = data.R.SWR.event';
      data.R.SWR.power  = data.R.SWR.power';
      
      if isfield(data,'gamma')
        data.gamma.SWR.event = data.gamma.SWR.event';
        data.gamma.SWR.power = data.gamma.SWR.power';
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
if param.spectOption
  fprintf(['spectral analysis of total LFP signal (file ' dataFileName ')... ']);
  fRange = param.spectLim1 : param.spectLim2;
  [data.LFP, ~] = calcSpect(data.LFP, [], fRange, data.param.Fs, 30, 0);
  fprintf('done\n');
  
  % If SWR events analyzed, detect spectrogram for event-locked data
  if (param.swrOption)
    fprintf(['spectral analysis of SWR-locked events (file ' dataFileName ')... ']);
    [data.SWR, ~] = calcSpect(data.SWR, [], fRange, data.param.Fs, 3, 0);
    data.SWR = calcEvFFT(data.SWR, data.param, param.spectLim1, param.spectLim2);
    fprintf('done\n');
    
    % Compute FFTs and phase over time for SWR-locked ripple & gamma (if selected)
    data.R.SWR = calcEvFFT(data.R.SWR, data.param, data.R.lim1, data.R.lim2);
    data.R.SWR = calcEvPhase(data.R.SWR, data.SWR, data.R.lim1, data.R.lim2);

    if isfield(data,'gamma')
      data.gamma.SWR = calcEvFFT(data.gamma.SWR, data.param, data.gamma.lim1, data.gamma.lim2);
      data.gamma.SWR = calcEvPhase(data.gamma.SWR, data.SWR, data.gamma.lim1, data.gamma.lim2);
    end
    
    if isfield(data,'fR')
      data.fR.SWR = calcEvFFT(data.fR.SWR, data.param, data.fR.lim1, data.fR.lim2);
      data.fR.SWR = calcEvPhase(data.fR.SWR, data.SWR, data.fR.lim1, data.fR.lim2);
    end
    
  end
end

% Re-order structure arrays
data       = orderStruct(data);
data.param = orderStruct(data.param);
data.LFP   = orderStruct(data.LFP);
data.LFP.param = orderStruct(data.LFP.param);

if (param.swOption)
  data.SW = orderStruct(data.SW);
  if isfield(data.SW, 'SWR') data.SW.SWR = orderStruct(data.SW.SWR); end
end

if (param.rOption)
  data.R = orderStruct(data.R);
  if isfield(data.R, 'SWR') data.R.SWR = orderStruct(data.R.SWR); end
end

if (param.swrOption)    data.SWR    = orderStruct(data.SWR); end
if (param.thetaOption)  data.theta  = orderStruct(data.theta); end
if (param.betaOption)   data.beta   = orderStruct(data.beta); end

if (param.gammaOption)
  data.gamma = orderStruct(data.gamma);
  if isfield(data.gamma, 'SWR') data.gamma.SWR = orderStruct(data.gamma.SWR); end
end

if (param.hgammaOption) data.hgamma = orderStruct(data.hgamma); end

if (param.fROption)
  data.fR = orderStruct(data.fR);
  if isfield(data.fR, 'SWR') data.fR.SWR = orderStruct(data.fR.SWR); end
end

if (param.cellOption)
  data.C = orderStruct(data.C);
  if isfield(data.C, 'SWR') data.C.SWR = orderStruct(data.C.SWR); end
end

%% Save file
if all(saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(saveFile,'-struct','data');
  fprintf('done\n');
end

%% Export SWR event table
if all(expEvFile) && param.expSWREvOption && param.swrOption
  fprintf(['exporting SWR events (file ' dataFileName ')... ']);
  exportSWREvents(data, saveFile, expEvFile);
  data.SWR.expEvFile = expEvFile;
  fprintf('done\n');
end

%% Export SWR event-locked episodic data files
if all(expDataFile) && param.expSWRDataOption && param.swrOption
  fprintf(['exporting SWR event-locked data (file ' dataFileName ')... ']);
  exportSWRData(data, param, expDataFile);
  data.SWR.expDataFile = expDataFile;
  fprintf('done\n');
end

if (param.fileType == 2)
  cd (curPath);
end

end
