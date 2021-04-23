function [data, hand] = analyzeSpkFile(data, hand, param, saveFile, spkFile, bstFile, expSpkFile, expBstFile, expSWRFile, expAveFile)
%% [data, hand] = analyzeSpkFile(data, hand, param, saveFile, spkFile, bstFile, expSpkFile, expBstFile, expSWRFile)
%
%  Function to detect coincidence of spikes/bursts and SWR events, construct 
%  peri-SWR spike histograms and spike-phase analysis
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   data       = data structure containing analyzed LFP files
%   hand       = handle structure to specify where figure should be drawn
%   param      = structure containing all parameters including:
%     param.fileNum              = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.importSpkOption      = boolean flag to import spike file (needed unless reanalyzing) (default = 1)
%     param.swrSpkOption         = boolean flag to calculate coincidence of SWRs and spikes (default = 1)
%     param.swrBstOption         = boolean flag to calculate coincidence of SWRs and bursts (default = 1)
%     param.lfpSpkOption         = boolean flag to calculate spike-phase over duration of recording (instead of swrSpkOption - only 1 selectable, default = 0)
%     param.useSWRDurationOption = boolean flag to use detected SWR detection for coincidence detection (default = 1)
%     param.useSWRWindowOption   = boolean flag to use standard swrWindow for coincidence detection (default = 0)
%     param.swrWindow            = +/- window around SWR peak events (default = 100 ms)
%     param.spkPhaseSWROption    = boolean flag to calculate spike-phase coupling within SWR events (default = 1)
%     param.expSWREvOption       = boolean flag to determine whether to export csv table of SWR events (default = 1)
%     param.spkPhaseLFPOption    = boolean flag to calculate spike-phase coupling over duration of file (default = 0)
%     param.spkAmpLFPOption      = boolean flag to calculate spike-amplitude coupling over duration of file (default = 0)
%     param.sdMultPhase          = for spike-phase coupling, the SD multiple the peak-trough amplitude must exceed in the oscillation surrounding the spike to be be considered
%     param.calcPhaseStats       = Calculates statistics, WARNING: requires circ_stat toolbox (default = 1)
%     param.expSpkEvOption       = boolean flag to determine whether to export csv table of Spk events (default = 1)
%     param.expBstEvOption       = boolean flag to determine whether to export csv table of Bst events (default = 1)
%     param.expAveOption         = boolean flag to determine whether to export csv table of average statistics
%     param.reAnalyzeOption      = option to re-analyze file (default = 0)
%     param.nBins                = For spike histogram, not currently selectable from UI (default = 100)
%   saveFile    = full path to matlab file to save (if not set, will prompt)
%   spkFile     = full path to pClamp spike event file to import (if not set, will prompt)
%   bstFile     = full path to pClamp burst event file to import (if not set, will prompt)
%   expSpkFile  = full path to spike event csv file to export (if not set, will prompt)
%   expBstFile  = full path to burst event csv file to export (if not set, will prompt)
%   expSWRFile  = full path to SWR event csv file to export (if not set, will prompt)
%   expAveFile  = full path to file of exported csv table of average (if not set and expAveOption = 1, will prompt)
%
%  Outputs:
%   data       = structure containing all data to be saved
%   hand       = handle structure for figure

%% Handle input arguments - if not entered
if (nargin < 10); expAveFile = []; end
if (nargin < 9); expSWRFile = []; end
if (nargin < 8); expBstFile = []; end
if (nargin < 7); expSpkFile = []; end
if (nargin < 6); bstFile    = []; end
if (nargin < 5); spkFile    = []; end
if (nargin < 4); saveFile   = []; end
if (nargin < 3); param      = struct; end
if (nargin < 2); hand       = struct; end
if (nargin < 1); data       = struct; end

% Handle case in which empty variables are supplied:
if isempty(param); param    = struct; end
if isempty(hand);  hand     = struct; end
if isempty(data);  data     = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum');              param.fileNum              = 1;   end
if ~isfield(param,'importSpkOption');      param.importSpkOption      = 1;   end
if ~isfield(param,'swrSpkOption');         param.swrSpkOption         = 1;   end
if ~isfield(param,'swrBstOption');         param.swrBstOption         = 1;   end
if ~isfield(param,'lfpSpkOption');         param.lfpSpkOption         = 0;   end
if ~isfield(param,'useSWRDurationOption'); param.useSWRDurationOption = 1;   end
if ~isfield(param,'useSWRWindowOption');   param.useSWRWindowOption   = 0;   end
if ~isfield(param,'swrWindow');            param.swrWindow            = 100; end
if ~isfield(param,'spkPhaseSWROption');    param.spkPhaseSWROption    = 1;   end
if ~isfield(param,'expSWREvOption');       param.expSWREvOption       = 1;   end
if ~isfield(param,'spkPhaseLFPOption');    param.spkPhaseLFPOption    = 0;   end
if ~isfield(param,'spkAmpLFPOption');      param.spkAmpLFPOption      = 0;   end
if ~isfield(param,'sdMultPhase');          param.sdMultPhase          = 2;   end
if ~isfield(param,'calcPhaseStats');       param.calcPhaseStats       = 1;   end
if ~isfield(param,'expSpkEvOption');       param.expSpkEvOption       = 1;   end
if ~isfield(param,'expBstEvOption');       param.expBstEvOption       = 1;   end
if ~isfield(param,'expAveOption');         param.expAveOption         = 1;   end
if ~isfield(param,'reAnalyzeOption');      param.reAnalyzeOption      = 0;   end
if ~isfield(param,'nBins');                param.nBins                = 100; end % Not yet selectable in UI

if ~isfield(data, 'LFP')
  [fileName, filePath] = uigetfile('.mat', 'Select *.mat file of analyzed LFP + imported cell channel');
  dataFile = [filePath fileName];
  if ~all(dataFile); error('No file data file selected'); end
  data = load(dataFile);
end

% Check if necessary data structures are present
if ~isfield(data,'LFP'); error('Must import LFP channel to data structure before proceeding'); end
if ~isfield(data,'C');   error('Must import cell channel to data structure before proceeding'); end
if ~isfield(data.C,'spike'); data.C.spike = struct; end

if param.swrSpkOption
  if ~isfield(data,'SWR'); error('Must analyze SWRs before proceeding'); end
  if ~isfield(data.SWR,'spike'); data.SWR.spike = struct; end
  if ~isfield(data.C,'SWR'); data.C.SWR = struct; end
  
  if param.swrBstOption
    if ~isfield(data.C,'burst'); data.C.burst = struct; end
    if ~isfield(data.SWR,'burst'); data.SWR.burst = struct; end
  end

  if ~isfield(data.C.SWR,'spike'); data.C.SWR.spike = struct; end
end

% Parse dataFile to determine default save name
[parentPath, dataFileName, ~] = parsePath(data.saveFile);

% Prompt for change in saveName
if isempty(saveFile)
  [saveName, savePath] = uiputfile('.mat','Select same or alternate *.mat file to save to', data.saveFile);
  saveFile = [savePath saveName];
  if ~all(saveFile)
    warning('No file to be saved - no file selected');
  else
    [parentPath, ~, ~] = parsePath(saveFile);
  end
end

% Prompt for event file(s) exported from pClamp if not supplied
if param.importSpkOption && ~param.reAnalyzeOption
  if isempty(spkFile)
    [fileName, filePath] = uigetfile('.txt', 'Select spike event *.txt file exported from pClamp', parentPath);
    spkFile = [filePath fileName];
    if ~all(spkFile); error('No spike file selected'); end
    data.C.spike.evFile = spkFile;
  end
  [parentPath, spkFileName, ~] = parsePath(spkFile);
  
  if param.swrBstOption
    if isempty(bstFile)
      [fileName, filePath] = uigetfile('.txt', 'Select burst event *.txt file exported from pClamp', parentPath);
      bstFile = [filePath fileName];
      if ~all(bstFile); error('No burst file selected'); end
      data.C.burst.evFile = bstFile;
    end
    [parentPath, bstFileName, ~] = parsePath(bstFile);
  end
end

% Select export spike event file, if option selected
if param.expSpkEvOption
  if isempty(expSpkFile)
    defaultPath = [parentPath dataFileName '_spkEvents.csv'];
    [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of spike events', defaultPath);
    expSpkFile = [exportPath exportName];
    if ~all(expSpkFile)
      warning('No spike events to be exported - no file selected');
    else
      [parentPath, expSpkFileName, ~] = parsePath(expSpkFile);
    end
  else
    [parentPath, expSpkFileName, ~] = parsePath(expSpkFile);
  end
end

% Select export burst event file, if option selected
if param.expBstEvOption
  if isempty(expBstFile)
    defaultPath = [parentPath dataFileName '_bstEvents.csv'];
    [expBstName, exportPath] = uiputfile('.csv','Select *.csv file to export table of burst events', defaultPath);
    expBstFile = [exportPath expBstName];
    if ~all(expBstFile)
      warning('No burst events to be exported - no file selected');
    else
      [parentPath, expBstFileName, ~] = parsePath(expBstFile);
    end
  else
    [parentPath, expBstFileName, ~] = parsePath(expBstFile);
  end
end

% Select re-export file for SWR events (if selected)
if param.expSWREvOption
  if isempty(expSWRFile)
    defaultName = [parentPath dataFileName '_swrEvents.csv'];
    [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export updated table of SWR events', defaultName);
    expSWRFile = [exportPath exportName];
    if ~all(expSWRFile)
      warning('No SWR events to be exported - no file selected');
    else
      [parentPath, expSWRFileName, ~] = parsePath(expSWRFile);
    end
  else
    [parentPath, expSWRFileName, ~] = parsePath(expSWRFile);
  end
end

% Select export average statistics file, if option selected
if param.expAveOption
  if isempty(expAveFile)
    defaultName = [parentPath dataFileName '_aveStats.csv'];
    [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of average statistics', defaultName);
    expAveFile = [exportPath exportName];
    if ~all(expAveFile)
      warning('No average statistics to be exported - no file selected');
    else
      [~, expAveFileName, ~] = parsePath(expAveFile);
    end
  else
    [~, expAveFileName, ~] = parsePath(expAveFile);
  end
end

% Save information to data structures
[~, saveName, saveExt] = parsePath(saveFile);
data.saveFile = saveFile;
data.saveName = [saveName '.' saveExt];
data.param    = param;
data.C.param  = param; % Save to Cell structure, as subsequent analysis may alter data.param

%% Import and process event file(s)
if param.importSpkOption && ~param.reAnalyzeOption
  
  % Import event files into matlab tables
  warning('off')
  fprintf(['importing spike event file ' spkFileName '... ']);
  spikeTable = readtable(spkFile,'Delimiter','\t','ReadRowNames', 1, 'TreatAsEmpty', 'Not found');
  fprintf('done\n');
  if param.swrBstOption
    fprintf(['importing burst event file ' bstFileName '... ']);
    burstTable = readtable(bstFile,'Delimiter','\t','ReadRowNames', 1, 'TreatAsEmpty', 'Not found');
    fprintf('done\n');
  end
  warning('on')
  
  % Re-initialize spike data structures
  data.C.spike.evStatus = [];
  data.C.spike.evStart  = [];
  data.C.spike.evPeak   = [];
  data.C.spike.evEnd    = [];
  
  fprintf(['processing spike events (file ' dataFileName ')... ']);
  for i = 1:size(spikeTable,2)
    if ~isempty(strfind(spikeTable.Properties.VariableNames{i},"EventStartTime"))
      spikeStartTime = spikeTable{:,i};
      for ev = 1:length(spikeStartTime)
        if ~isempty(find(data.C.timing >= spikeStartTime(ev), 1))
          data.C.spike.evStart(ev) = find(data.C.timing >= spikeStartTime(ev), 1);
        end
      end
      data.C.spike.evStart = data.C.spike.evStart';
    end
    
    if ~isempty(strfind(spikeTable.Properties.VariableNames{i},"EventEndTime"))
      spikeEndTime = spikeTable{:,i};
      for ev = 1:length(spikeEndTime)
        if ~isempty(find(data.C.timing <= spikeEndTime(ev), 1, 'last'))
          data.C.spike.evEnd(ev) = find(data.C.timing <= spikeEndTime(ev), 1, 'last');
        end
      end
      data.C.spike.evEnd = data.C.spike.evEnd';
    end
    
    if ~isempty(strfind(spikeTable.Properties.VariableNames{i},"TimeOfPeak"))
      spikePeakTime = spikeTable{:,i};
      for ev = 1:length(spikePeakTime)
        if ~isempty(find(data.C.timing >= spikePeakTime(ev), 1, 'last'))
          data.C.spike.evPeak(ev) = find(data.C.timing >= spikePeakTime(ev), 1);
        end
      end
      data.C.spike.evPeak = data.C.spike.evPeak';
    end
  end
  data.C.spike.nEvents = length(data.C.spike.evStart);
  data.C.spike.frequency = data.C.spike.nEvents / ((data.C.timing(end) - data.C.timing(1)) / 1000);
  fprintf('done\n');
  
  % compute spike status array
  data.C.spike.evStatus = zeros(length(data.C.timing),1);
  for ev = 1:length(data.C.spike.evStart)
    for i = data.C.spike.evStart(ev) : data.C.spike.evEnd(ev)
      data.C.spike.evStatus(i) = 1;
    end
  end
  

  %% Import burst file (if selected)
  if param.swrBstOption
    
    % Re-initialize data structures
    data.C.burst.nSpike   = [];
    data.C.burst.evStart  = [];
    data.C.burst.evEnd    = [];
    data.C.burst.intraBI  = [];
    data.C.burst.evStatus = [];
    
    fprintf(['processing burst events (file ' dataFileName ')... ']);
    for i = 1:size(burstTable,2)
      if ~isempty(strfind(burstTable.Properties.VariableNames{i},"EventsInBurst"))
        data.C.burst.nSpike = burstTable{:,i};
      end
      
      if ~isempty(strfind(burstTable.Properties.VariableNames{i},"StartTime"))
        burstStartTime = burstTable{:,i};
        for ev = 1:length(burstStartTime)
          if ~isempty(find(data.C.timing >= burstStartTime(ev), 1))
            data.C.burst.evStart(ev) = find(data.C.timing >= burstStartTime(ev), 1);
          end
        end
        data.C.burst.evStart = data.C.burst.evStart';
      end
      
      if ~isempty(strfind(burstTable.Properties.VariableNames{i},"EndTime"))
        burstEndTime = burstTable{:,i};
        for ev = 1:length(burstEndTime)
          if ~isempty(find(data.C.timing <= burstEndTime(ev), 1, 'last'))
            data.C.burst.evEnd(ev) = find(data.C.timing <= burstEndTime(ev), 1, 'last');
          end
        end
        data.C.burst.evEnd = data.C.burst.evEnd';
      end
      
      if ~isempty(strfind(burstTable.Properties.VariableNames{i},"MeanIntraburstInterval"))
        data.C.burst.intraBI = burstTable{:,i};
      end
    end
    data.C.burst.nEvents   = length(data.C.burst.evStart);
    data.C.burst.frequency = data.C.burst.nEvents / ((data.C.timing(end) - data.C.timing(1)) / 1000);
    fprintf('done\n');
    
    % compute burst status array
    data.C.burst.evStatus = zeros(length(data.C.timing),1);
    for ev = 1:length(data.C.burst.evStart)
      for i = data.C.burst.evStart(ev) : data.C.burst.evEnd(ev)
        data.C.burst.evStatus(i) = 1;
      end
    end
  end
end

%% Calculate standard SWR status (if selected)
if param.useSWRWindowOption && param.swrSpkOption
  if ~isfield(data.SWR,'evStatusStand')
    data.SWR.evStatusStand = zeros(length(data.SWR.evStatus),1);
    if ~isempty(data.SWR.evStart)
      for swr = 1:length(data.SWR.evStart)
        loBaseWin = max(round(data.SWR.evPeak(swr) - 0.5 * param.swrWindow / data.LFP.samplingInt), 1);
        hiBaseWin = min(round(data.SWR.evPeak(swr) + 0.5 * param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
        data.SWR.evStatusStand(loBaseWin:hiBaseWin) = 1;
      end
    end
  end
end

%% Spike-SWR Analysis
if param.swrSpkOption
  
  %% Align SWR-Spike Arrays:
  fprintf(['aligning time arrays for SWR-spike coincidence analysis (file ' dataFileName ')... ']);
  if param.useSWRDurationOption
    [data.SWR.spike.evStatusA, data.SWR.spike.evStartA, data.SWR.spike.evEndA, data.C.spike.evStatusA, data.C.spike.evStartA, data.C.spike.evEndA, data.SWR.spike.timingA] = ...
      timeAlign(data.SWR.evStatus, data.C.spike.evStatus, data.LFP.timing, data.C.timing);
    % Re-calculate peak - truncating if necessary:
    data.SWR.spike.evPeakA = data.SWR.evPeak(1:length(data.SWR.spike.evStartA));
    
  elseif param.useSWRWindowOption
    [data.SWR.spike.evStatusA, data.SWR.spike.evStartA, data.SWR.spike.evEndA, data.C.spike.evStatusA, data.C.spike.evStartA, data.C.spike.evEndA, data.SWR.spike.timingA] = ...
      timeAlign(data.SWR.evStatusStand, data.C.spike.evStatus, data.LFP.timing, data.C.timing);
    % Re-calculate SWR peak from mid-point of start and end:
    data.SWR.spike.evPeakA = round(0.5*(data.SWR.spike.evStartA + data.SWR.spike.evEndA));
    
  end
  data.C.spike.evPeakA = data.C.spike.evPeak(1:length(data.C.spike.evStartA));
  
  % Align SWR-Burst Arrays:
  if param.swrBstOption
    if param.useSWRDurationOption
      [data.SWR.burst.evStatusA, data.SWR.burst.evStartA, data.SWR.burst.evEndA, data.C.burst.evStatusA, data.C.burst.evStartA, data.C.burst.evEndA, data.SWR.burst.timingA] = ...
        timeAlign(data.SWR.evStatus, data.C.burst.evStatus, data.LFP.timing, data.C.timing);
      % Re-calculate SWR peak - truncating if necessary:
      data.SWR.burst.evPeakA = data.SWR.evPeak(1:length(data.SWR.burst.evStartA));
      
    elseif param.useSWRWindowOption
      [data.SWR.burst.evStatusA, data.SWR.burst.evStartA, data.SWR.burst.evEndA, data.C.burst.evStatusA, data.C.burst.evStartA, data.C.burst.evEndA, data.SWR.burst.timingA] = ...
        timeAlign(data.SWR.evStatusStand, data.C.burst.evStatus, data.LFP.timing, data.C.timing);
      % Re-calculate SWR peak from mid-point of start and end:
      data.SWR.burst.evPeakA = round(0.5*(data.SWR.burst.evStartA + data.SWR.burst.evEndA));
      
    end
    % In case any burst events truncated:
    data.C.burst.nSpike  = data.C.burst.nSpike(1:length(data.C.burst.evStartA));
    data.C.burst.intraBI = data.C.burst.intraBI(1:length(data.C.burst.evStartA));
    
  end
  fprintf('done\n');
  
  %% Parse events into SWR blocks
  fprintf(['parsing spikes into SWR events (file ' dataFileName ')... ']);
  data.C.SWR.spike.evStartA  = [];
  data.C.SWR.spike.evPeakA   = [];
  data.C.SWR.spike.evEndA    = [];
  data.C.SWR.spike.evStatusA = [];
  
  if ~isempty(data.SWR.evStart)
    
    % Initialize event locked data window cell arrays
    data.C.SWR.spike.evStartA{length(data.SWR.spike.evStartA)}  = [];
    data.C.SWR.spike.evPeakA{length(data.SWR.spike.evStartA)}   = [];
    data.C.SWR.spike.evEndA{length(data.SWR.spike.evStartA)}    = [];
    data.C.SWR.spike.evStatusA{length(data.SWR.spike.evStartA)} = [];
    
    % For each SWR event, parse spike events
    for swr = 1:length(data.SWR.spike.evStartA)

      % Re-calculate start and end windows for each SWR event:
      loWin = max(round(data.SWR.spike.evPeakA(swr) - param.swrWindow / data.LFP.samplingInt), 1);
      hiWin = min(round(data.SWR.spike.evPeakA(swr) + param.swrWindow / data.LFP.samplingInt), length(data.C.spike.evStatusA));
      
      % Parse spike events into SWR cell array
      data.C.SWR.spike.evStartA{swr}  = data.C.spike.evStartA(data.C.spike.evStartA>=loWin & data.C.spike.evEndA<=hiWin) - (loWin - 1);
      data.C.SWR.spike.evPeakA{swr}   = data.C.spike.evPeakA(data.C.spike.evStartA>=loWin & data.C.spike.evEndA<=hiWin) - (loWin - 1);
      data.C.SWR.spike.evEndA{swr}    = data.C.spike.evEndA(data.C.spike.evStartA>=loWin & data.C.spike.evEndA<=hiWin) - (loWin - 1);
      data.C.SWR.spike.evStatusA{swr} = data.C.spike.evStatusA(loWin : hiWin);

    end

    data.C.SWR.spike.evStartA  = data.C.SWR.spike.evStartA';
    data.C.SWR.spike.evPeakA   = data.C.SWR.spike.evPeakA';
    data.C.SWR.spike.evEndA    = data.C.SWR.spike.evEndA';
    data.C.SWR.spike.evStatusA = data.C.SWR.spike.evStatusA';
    
    % Calculate summed event spike status and normalized histogram over all SWR events 
    nSamples = max(cellfun(@length, data.C.SWR.spike.evStatusA));
    Bins     = linspace(data.SWR.evTiming(1), data.SWR.evTiming(length(data.SWR.evTiming)), param.nBins)';
    binWidth = Bins(2) - Bins(1);
    
    data.C.SWR.spike.evStatusSum = zeros(nSamples, 1);
    data.C.SWR.spike.evHist      = zeros(param.nBins, 1);
    data.C.SWR.spike.nSpksHist   = 0;
    
    for swr = 1:length(data.C.SWR.spike.evStatusA)
      
      if length(data.C.SWR.spike.evStatusA{swr}) == nSamples
        data.C.SWR.spike.evStatusSum = data.C.SWR.spike.evStatusSum + data.C.SWR.spike.evStatusA{swr};
      end
      
      nSpkSWR = length(data.C.SWR.spike.evPeakA{swr});
      data.C.SWR.spike.nSpksHist = data.C.SWR.spike.nSpksHist + nSpkSWR;
      
      for spk = 1:nSpkSWR
        peakTime = data.SWR.evTiming(data.C.SWR.spike.evPeakA{swr}(spk));
        data.C.SWR.spike.evHist = data.C.SWR.spike.evHist + double(peakTime >= Bins & peakTime < Bins + binWidth);
      end
    end
    % Normalize histogram to all SWRs to get prob. of spiking per bin
    data.C.SWR.spike.evHist = data.C.SWR.spike.evHist / length(data.C.SWR.spike.evPeakA); 
  end
  fprintf('done\n');

  %% Calculate overlap of SWRs and spikes
  fprintf(['detecting spikes coincident with SWRs (file ' dataFileName ')... ']);
  [data.C.spike.evStatusC, data.C.spike.evStartC, data.C.spike.evEndC, data.C.spike.evIndex] = eventOverlap(data.SWR.spike.evStatusA, data.SWR.spike.evStartA, data.SWR.spike.evEndA, ...
    data.C.spike.evStatusA, data.C.spike.evStartA, data.C.spike.evEndA, data.SWR.spike.timingA, 2);
  fprintf('done\n');
  
  % Find SWRs that coincide with at least one spike
  fprintf(['detecting SWRs coincident with spikes (file ' dataFileName ')... ']);
  [data.SWR.spike.evStatusC, data.SWR.spike.evStartC, data.SWR.spike.evEndC, data.SWR.spike.evIndex] = eventOverlap(data.SWR.spike.evStatusA, data.SWR.spike.evStartA, data.SWR.spike.evEndA, ...
    data.C.spike.evStatusA, data.C.spike.evStartA, data.C.spike.evEndA, data.SWR.spike.timingA, 1);
  fprintf('done\n');
  
  data.SWR.spike.nEventsA = size(data.SWR.spike.evStartA,1);
  data.SWR.spike.nEventsC = size(data.SWR.spike.evStartC,1);
  data.C.spike.nEventsA   = size(data.C.spike.evStartA,1);
  data.C.spike.nEventsC   = size(data.C.spike.evStartC,1);
  
  % Calculate simplified event matrix
  data.C.spike.swrMatrix  = zeros(data.C.spike.nEventsA, 1);
  data.SWR.spike.evMatrix = zeros(data.SWR.spike.nEventsA, 1);
  
  evStatusC = data.C.spike.evStatusA .* data.SWR.spike.evStatusA;
  
  for spk = 1:data.C.spike.nEventsA
    if (sum(evStatusC(data.C.spike.evStartA(spk) : data.C.spike.evEndA(spk))) > 0)
      data.C.spike.swrMatrix(spk) = 1;
    end
  end
  
  for swr = 1:data.SWR.spike.nEventsA
    if (sum(evStatusC(data.SWR.spike.evStartA(swr) : data.SWR.spike.evEndA(swr))) > 0)
      data.SWR.spike.evMatrix(swr) = 1;
    end
  end
  
  %% Calculate overlap of SWRs and bursts
  if param.swrBstOption
    
    % Calculate overlap for SWRs and spikes
    fprintf(['detecting bursts coincident with SWRs (file ' dataFileName ')... ']);
    [data.C.burst.evStatusC, data.C.burst.evStartC, data.C.burst.evEndC] = eventOverlap(data.SWR.burst.evStatusA, data.SWR.burst.evStartA, data.SWR.burst.evEndA, ...
      data.C.burst.evStatusA, data.C.burst.evStartA, data.C.burst.evEndA, data.SWR.burst.timingA, 2);
    fprintf('done\n');
    
    % Find SWRs that coincide with at least one spike
    fprintf(['detecting SWRs coincident with bursts (file ' dataFileName ')... ']);
    [data.SWR.burst.evStatusC, data.SWR.burst.evStartC, data.SWR.burst.evEndC] = eventOverlap(data.SWR.burst.evStatusA, data.SWR.burst.evStartA, data.SWR.burst.evEndA, ...
      data.C.burst.evStatusA, data.C.burst.evStartA, data.C.burst.evEndA, data.SWR.burst.timingA, 1);
    fprintf('done\n');
    
    data.SWR.burst.nEventsA = size(data.SWR.burst.evStartA,1);
    data.SWR.burst.nEventsC = size(data.SWR.burst.evStartC,1);
    data.C.burst.nEventsA   = size(data.C.burst.evStartA,1);
    data.C.burst.nEventsC   = size(data.C.burst.evStartC,1);
    
    % Calculate event matrix
    data.C.burst.swrMatrix  = zeros(data.C.burst.nEventsA, 1);
    data.SWR.burst.evMatrix = zeros(data.SWR.burst.nEventsA, 1);
    
    evStatusC = data.C.burst.evStatusA .* data.SWR.burst.evStatusA;
    
    for bst = 1:data.C.burst.nEventsA
      if (sum(evStatusC(data.C.burst.evStartA(bst) : data.C.burst.evEndA(bst))) > 0)
        data.C.burst.swrMatrix(bst) = 1;
      end
    end
    
    for swr = 1:data.SWR.burst.nEventsA
      if (sum(evStatusC(data.SWR.burst.evStartA(swr) : data.SWR.burst.evEndA(swr))) > 0)
        data.SWR.burst.evMatrix(swr) = 1;
      end
    end
  end
  
  %% Calculate spike-phase within SWRs
  if param.spkPhaseSWROption && ~isempty(data.SWR.spike.evStartA) && isfield(data.R.SWR.phase, 'evPhase') 
    fprintf(['detecting oscillation phase of spikes during SWRs (file ' dataFileName ')... ']);

    % Gamma:
    if isfield(data, 'gamma')
      if isfield(data.gamma, 'SWR')
        if isfield(data.gamma.SWR, 'phase')
          data.C.SWR.spike.gamma = struct;
          data.C.SWR.spike.gamma.phase = [];
          data.C.SWR.spike.gamma.phase{length(data.C.SWR.spike.evStartA)} = [];
          data.C.SWR.spike.gamma.phase = data.C.SWR.spike.gamma.phase';
          data.C.spike.gamma.phase  = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.gamma.phaseX = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.gamma.phaseY = NaN*ones(length(data.C.spike.evStartA), 1);
          minAmpG = std(data.gamma.tSeries) * param.sdMultPhase;
        end
      end
    end
    
    % High Gamma:
    if isfield(data, 'hgamma')
      if isfield(data.hgamma, 'SWR')
        if isfield(data.hgamma.SWR, 'phase')
          data.C.SWR.spike.hgamma = struct;
          data.C.SWR.spike.hgamma.phase = [];
          data.C.SWR.spike.hgamma.phase{length(data.C.SWR.spike.evStartA)} = [];
          data.C.SWR.spike.hgamma.phase = data.C.SWR.spike.hgamma.phase';
          data.C.spike.hgamma.phase  = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.hgamma.phaseX = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.hgamma.phaseY = NaN*ones(length(data.C.spike.evStartA), 1);
          minAmpHG = std(data.hgamma.tSeries) * param.sdMultPhase;
        end
      end
    end
    
    % Ripple:
    if isfield(data, 'R')
      if isfield(data.R, 'SWR')
        if isfield(data.R.SWR, 'phase')
          data.C.SWR.spike.R = struct;
          data.C.SWR.spike.R.phase = [];
          data.C.SWR.spike.R.phase{length(data.C.SWR.spike.evStartA)} = [];
          data.C.SWR.spike.R.phase = data.C.SWR.spike.R.phase';
          data.C.spike.R.phase  = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.R.phaseX = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.R.phaseY = NaN*ones(length(data.C.spike.evStartA), 1);
          minAmpR = std(data.R.tSeries) * param.sdMultPhase;
        end
      end
    end
    
    % Fast Ripple:
    if isfield(data, 'fR')
      if isfield(data.fR, 'SWR')
        if isfield(data.fR.SWR, 'phase')
          data.C.SWR.spike.fR = struct;
          data.C.SWR.spike.fR.phase = [];
          data.C.SWR.spike.fR.phase{length(data.C.SWR.spike.evStartA)} = [];
          data.C.SWR.spike.fR.phase = data.C.SWR.spike.fR.phase';
          data.C.spike.fR.phase  = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.fR.phaseX = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.fR.phaseY = NaN*ones(length(data.C.spike.evStartA), 1);
          minAmpFR = std(data.fR.tSeries) * param.sdMultPhase;
        end
      end
    end

    nSamples = max(cellfun(@length, data.SWR.event));
    
    for spk1 = 1:data.C.spike.nEventsA
      if (data.C.spike.swrMatrix(spk1) == 1) % Coincident spike
        spkInd = data.C.spike.evIndex(:, 2);
        swr    = data.C.spike.evIndex(spkInd == spk1, 1);
        loWin  = max(round(data.SWR.spike.evPeakA(swr) - param.swrWindow / data.LFP.samplingInt), 1);
        spk2   = find(data.C.SWR.spike.evStartA{swr} == data.C.spike.evStartA(spk1) - (loWin - 1));
        
        % Only consider non-truncated SWRs with spikes
        if(length(data.SWR.event{swr}) == nSamples) && ~isempty(spk2)
          
          % Gamma:
          if isfield(data.C.SWR.spike, 'gamma')
            ampG = calcEvPhaseAmp(data.gamma.SWR.phase, swr, data.C.SWR.spike.evPeakA{swr}(spk2));
            if ampG > minAmpG
              data.C.SWR.spike.gamma.phase{swr}(spk2) = data.gamma.SWR.phase.evPhase{swr}(data.C.SWR.spike.evPeakA{swr}(spk2));
              data.C.spike.gamma.phase(spk1) = data.C.SWR.spike.gamma.phase{swr}(spk2);
              data.C.spike.gamma.phaseX(spk1) = cos(data.C.spike.gamma.phase(spk1));
              data.C.spike.gamma.phaseY(spk1) = sin(data.C.spike.gamma.phase(spk1));
            end
          end
          
          % High Gamma:
          if isfield(data.C.SWR.spike, 'hgamma')
            ampHG = calcEvPhaseAmp(data.hgamma.SWR.phase, swr, data.C.SWR.spike.evPeakA{swr}(spk2));
            if ampHG > minAmpHG
              data.C.SWR.spike.hgamma.phase{swr}(spk2) = data.hgamma.SWR.phase.evPhase{swr}(data.C.SWR.spike.evPeakA{swr}(spk2));
              data.C.spike.hgamma.phase(spk1)  = data.C.SWR.spike.hgamma.phase{swr}(spk2);
              data.C.spike.hgamma.phaseX(spk1) = cos(data.C.spike.hgamma.phase(spk1));
              data.C.spike.hgamma.phaseY(spk1) = sin(data.C.spike.hgamma.phase(spk1));
            end
          end
          
          % Ripple:
          if isfield(data.C.SWR.spike, 'R')
            ampR = calcEvPhaseAmp(data.R.SWR.phase, swr, data.C.SWR.spike.evPeakA{swr}(spk2));
            if ampR > minAmpR
              data.C.SWR.spike.R.phase{swr}(spk2) = data.R.SWR.phase.evPhase{swr}(data.C.SWR.spike.evPeakA{swr}(spk2));
              data.C.spike.R.phase(spk1)  = data.C.SWR.spike.R.phase{swr}(spk2);
              data.C.spike.R.phaseX(spk1) = cos(data.C.spike.R.phase(spk1));
              data.C.spike.R.phaseY(spk1) = sin(data.C.spike.R.phase(spk1));
            end
          end 
          
          % Fast Ripple:
          if isfield(data.C.SWR.spike, 'fR')
            ampFR = calcEvPhaseAmp(data.fR.SWR.phase, swr, data.C.SWR.spike.evPeakA{swr}(spk2));
            if ampFR > minAmpFR
              data.C.SWR.spike.gamma.phase{swr}(spk2) = data.fR.SWR.phase.evPhase{swr}(data.C.SWR.spike.evPeakA{swr}(spk2));
              data.C.spike.fR.phase(spk1)  = data.C.SWR.spike.fR.phase{swr}(spk2);
              data.C.spike.fR.phaseX(spk1) = cos(data.C.spike.fR.phase(spk1));
              data.C.spike.fR.phaseY(spk1) = sin(data.C.spike.fR.phase(spk1));
            end
          end
        end
      end
    end
    fprintf('done\n');
  end
end


%% Spike-LFP Analysis
if param.lfpSpkOption
  
  %% Calculate spike-phase coupling of entire file (if option selected)
  if param.spkPhaseLFPOption
    fprintf(['detecting oscillation phase of spikes over entire file (' dataFileName ')... ']);
    
    % Theta:
    if isfield(data, 'theta')
      if isfield(data.theta, 'phase')
        data.C.spike.theta = struct;
        data.C.spike.theta.phase  = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.theta.phaseX = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.theta.phaseY = NaN*ones(length(data.C.spike.evStart), 1);
        minAmpT = std(data.theta.tSeries) * param.sdMultPhase;
      end
    end
    
    % Alpha:
    if isfield(data, 'alpha')
      if isfield(data.alpha, 'phase')
        data.C.spike.alpha = struct;
        data.C.spike.alpha.phase  = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.alpha.phaseX = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.alpha.phaseY = NaN*ones(length(data.C.spike.evStart), 1);
        minAmpA = std(data.alpha.tSeries) * param.sdMultPhase;
      end
    end
    
    % Beta:
    if isfield(data, 'beta')
      if isfield(data.beta, 'phase')
        data.C.spike.beta = struct;
        data.C.spike.beta.phase  = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.beta.phaseX = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.beta.phaseY = NaN*ones(length(data.C.spike.evStart), 1);
        minAmpB = std(data.beta.tSeries) * param.sdMultPhase;
      end
    end
    
    % Gamma:
    if isfield(data, 'gamma')
      if isfield(data.gamma, 'phase')
        data.C.spike.gamma = struct;
        data.C.spike.gamma.phase  = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.gamma.phaseX = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.gamma.phaseY = NaN*ones(length(data.C.spike.evStart), 1);
        minAmpG = std(data.gamma.tSeries) * param.sdMultPhase;
      end
    end
    
    % High Gamma:
    if isfield(data, 'hgamma')
      if isfield(data.hgamma, 'phase')
        data.C.spike.hgamma = struct;
        data.C.spike.hgamma.phase  = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.hgamma.phaseX = NaN*ones(length(data.C.spike.evStart), 1);
        data.C.spike.hgamma.phaseY = NaN*ones(length(data.C.spike.evStart), 1);
        minAmpHG = std(data.hgamma.tSeries) * param.sdMultPhase;
      end
    end
    
    %% Spike by spike, calculate phase if oscillation above threshold
    for spk = 1:data.C.spike.nEvents
      
      % Theta:
      if isfield(data.C.spike, 'theta')
        ampT = calcTotPhaseAmp(data.theta.phase, data.C.spike.evPeak(spk));
        if ampT > minAmpT
          data.C.spike.theta.phase(spk)  = data.theta.phase.tPhase(data.C.spike.evPeak(spk));
          data.C.spike.theta.phaseX(spk) = cos(data.C.spike.theta.phase(spk));
          data.C.spike.theta.phaseY(spk) = sin(data.C.spike.theta.phase(spk));
        end
      end
      
      % Alpha:
      if isfield(data.C.spike, 'alpha')
        ampA = calcTotPhaseAmp(data.alpha.phase, data.C.spike.evPeak(spk));
        if ampA > minAmpA
          data.C.spike.alpha.phase(spk)  = data.alpha.phase.tPhase(data.C.spike.evPeak(spk));
          data.C.spike.alpha.phaseX(spk) = cos(data.C.spike.alpha.phase(spk));
          data.C.spike.alpha.phaseY(spk) = sin(data.C.spike.alpha.phase(spk));
        end
      end
      
      % Beta:
      if isfield(data.C.spike, 'beta')
        ampB = calcTotPhaseAmp(data.beta.phase, data.C.spike.evPeak(spk));
        if ampB > minAmpB
          data.C.spike.beta.phase(spk)  = data.beta.phase.tPhase(data.C.spike.evPeak(spk));
          data.C.spike.beta.phaseX(spk) = cos(data.C.spike.beta.phase(spk));
          data.C.spike.beta.phaseY(spk) = sin(data.C.spike.beta.phase(spk));
        end
      end
      
      % Gamma:
      if isfield(data.C.spike, 'gamma')
        ampG = calcTotPhaseAmp(data.gamma.phase, data.C.spike.evPeak(spk));
        if ampG > minAmpG
          data.C.spike.gamma.phase(spk)  = data.gamma.phase.tPhase(data.C.spike.evPeak(spk));
          data.C.spike.gamma.phaseX(spk) = cos(data.C.spike.gamma.phase(spk));
          data.C.spike.gamma.phaseY(spk) = sin(data.C.spike.gamma.phase(spk));
        end
      end
      
      % High Gamma:
      if isfield(data.C.spike, 'hgamma')
        ampHG = calcTotPhaseAmp(data.hgamma.phase, data.C.spike.evPeak(spk));
        if ~isnan(ampHG) && ampHG > minAmpHG
          data.C.spike.hgamma.phase(spk)  = data.hgamma.phase.tPhase(data.C.spike.evPeak(spk));
          data.C.spike.hgamma.phaseX(spk) = cos(data.C.spike.hgamma.phase(spk));
          data.C.spike.hgamma.phaseY(spk) = sin(data.C.spike.hgamma.phase(spk));
        end
      end
    end
    fprintf('done\n');
  end
  
  %% Calculate spike-amplitude coupling of entire file (if option selected)
  if param.spkAmpLFPOption
    
  end
end
    
%% Calculate additional phase stats if option selected (requires circular statistics package)
if param.calcPhaseStats
  fprintf(['calculating circular phase statistics (file: ' dataFileName ')... ']);
  if isfield(data.C.spike, 'theta');  data.C.spike.theta  = calcPhaseStats(data.C.spike.theta);  end
  if isfield(data.C.spike, 'alpha');  data.C.spike.alpha  = calcPhaseStats(data.C.spike.alpha);  end
  if isfield(data.C.spike, 'beta');   data.C.spike.beta   = calcPhaseStats(data.C.spike.beta);   end
  if isfield(data.C.spike, 'gamma');  data.C.spike.gamma  = calcPhaseStats(data.C.spike.gamma);  end
  if isfield(data.C.spike, 'hgamma'); data.C.spike.hgamma = calcPhaseStats(data.C.spike.hgamma); end
  if isfield(data.C.spike, 'R');      data.C.spike.R      = calcPhaseStats(data.C.spike.R);      end
  if isfield(data.C.spike, 'fR');     data.C.spike.fR     = calcPhaseStats(data.C.spike.fR);     end
  fprintf('done\n');
end

% Re-order structure arrays:
data.C         = orderStruct(data.C);
data.C.spike   = orderStruct(data.C.spike);
data.C.param   = orderStruct(data.C.param);

if param.swrSpkOption
  data.SWR         = orderStruct(data.SWR);
  data.SWR.spike   = orderStruct(data.SWR.spike);
  data.C.SWR.spike = orderStruct(data.C.SWR.spike);
end

if param.swrBstOption
  data.C.burst   = orderStruct(data.C.burst);
  data.SWR.burst = orderStruct(data.SWR.burst);
end


%% Save and export results
% Export Spike event table
if all(expSpkFile) && param.expSpkEvOption
  fprintf(['exporting spike event file ' expSpkFileName '... ']);
  exportSpkEvents(data, saveFile, expSpkFile);
  data.C.expSpkFile = expSpkFile;
  fprintf('done\n');
end

% Export Burst event table
if all(expBstFile) && param.expBstEvOption
  fprintf(['exporting burst event file ' expBstFileName '... ']);
  exportBstEvents(data, saveFile, expBstFile);
  data.C.expBstFile = expBstFile;
  fprintf('done\n');
end

% Export SWR event file
if (all(expSWRFile) && param.expSWREvOption)
  fprintf(['exporting updated SWR events (file ' expSWRFileName ')... ']);
  exportSWREvents(data, saveFile, expSWRFile)
  data.SWR.expEvFile = expSWRFile;
  fprintf('done\n');
end

% Export average table
if all(expAveFile) && param.expAveOption
  fprintf(['exporting averages statistics (file ' expAveFileName ')... ']);
  exportAveStats(data, saveFile, expAveFile);
  data.C.expAveFile = expAveFile;
  fprintf('done\n');
end

% Save matlab file
if all(saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(saveFile,'-struct','data');
  fprintf('done\n');
end

end
