function [data, hand] = analyzeSpkFile(data, hand, param, saveFile, spkFile, bstFile, expSpkFile, expBstFile, expSWRFile)

%% Handle input arguments - if not entered
if (nargin < 9) expSWRFile = []; end
if (nargin < 8) expBstFile = []; end
if (nargin < 7) expSpkFile = []; end
if (nargin < 6) bstFile    = []; end
if (nargin < 5) spkFile    = []; end
if (nargin < 4) saveFile   = []; end
if (nargin < 3) param      = struct; end
if (nargin < 2) hand       = struct; end
if (nargin < 1) data       = struct; end

% Handle case in which empty variables are supplied:
if isempty(param) param    = struct; end
if isempty(hand)  hand     = struct; end
if isempty(data)  data     = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileNum')              param.fileNum = 1; end
if ~isfield(param,'importSpkOption')      param.importSpkOption = 1; end
if ~isfield(param,'swrSpkOption')         param.swrSpkOption = 1; end
if ~isfield(param,'swrBstOption')         param.swrBstOption = 1; end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 1; end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption = 0; end
if ~isfield(param,'swrWindow')            param.swrWindow = 100; end
if ~isfield(param,'parseSpkOption')       param.parseSpkOption = 1; end
if ~isfield(param,'calcEvMatrixOption')   param.calcEvMatrixOption = 1; end
if ~isfield(param,'expSpkEvOption')       param.expSpkEvOption = 1; end
if ~isfield(param,'expBstEvOption')       param.expBstEvOption = 1; end
if ~isfield(param,'expSWREvOption')       param.expSWREvOption = 1; end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption = 0; end

if ~isfield(data, 'LFP')
  [fileName, filePath] = uigetfile('.mat', 'Select *.mat file of analyzed LFP + imported cell channel');
  dataFile = [filePath fileName];
  if ~all(dataFile) error('No file data file selected'); end
  data = load(dataFile);
end

% Check if necessary data structures are present
if ~isfield(data,'LFP') error('Must import LFP channel to data structure before proceeding'); end
if ~isfield(data,'C')   error('Must import cell channel to data structure before proceeding'); end
if ~isfield(data.C,'spike') data.C.spike = struct; end

if param.swrSpkOption
  if ~isfield(data,'SWR') error('Must analyze SWRs before proceeding'); end
  if ~isfield(data.SWR,'spike') data.SWR.spike = struct; end
  if ~isfield(data.C,'SWR') data.C.SWR = struct; end
  
  if param.swrBstOption
    if ~isfield(data.C,'burst') data.C.burst = struct; end
    if ~isfield(data.SWR,'burst') data.SWR.burst = struct; end
  end
  
  if param.parseSpkOption
    if ~isfield(data.C.SWR,'spike') data.C.SWR.spike = struct; end
  end
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
    if ~all(spkFile) error('No spike file selected'); end
    data.C.spike.evFile = spkFile;
  end
  [parentPath, spkFileName, ~] = parsePath(spkFile);
  
  if param.swrBstOption
    if isempty(bstFile)
      [fileName, filePath] = uigetfile('.txt', 'Select burst event *.txt file exported from pClamp', parentPath);
      bstFile = [filePath fileName];
      if ~all(bstFile) error('No burst file selected'); end
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
      [~, expSWRFileName, ~] = parsePath(expSWRFile);
    end
  else
    [~, expSWRFileName, ~] = parsePath(expSWRFile);
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

%% Align files
if param.swrSpkOption
  
  % Align SWR-Spike Arrays:
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
end

%% Parse events into SWR blocks (currently mandatory)
if param.parseSpkOption
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
    
    % Calculate summed event spike status over all SWR events:
    nSamples = max(cellfun(@length, data.C.SWR.spike.evStatusA));
    data.C.SWR.spike.evStatusSum = zeros(nSamples, 1);
    for swr = 1:length(data.C.SWR.spike.evStatusA)
      if length(data.C.SWR.spike.evStatusA{swr}) == nSamples
        data.C.SWR.spike.evStatusSum = data.C.SWR.spike.evStatusSum + data.C.SWR.spike.evStatusA{swr};
      end
    end
    
  end
end

%% Calculate overlap of SWRs and spikes
if param.swrSpkOption
  
  % Calculate overlap for SWRs and spikes
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
  
  % Calculate simplified event matrix (currently mandatory)
  if param.calcEvMatrixOption
    data.C.spike.swrMatrix  = zeros(data.C.spike.nEventsA, 1);
    data.SWR.spike.evMatrix = zeros(data.SWR.spike.nEventsA, 1);
    
    coincStatus = data.C.spike.evStatusA .* data.SWR.spike.evStatusA;
    
    for spk = 1:data.C.spike.nEventsA
      if (sum(coincStatus(data.C.spike.evStartA(spk) : data.C.spike.evEndA(spk))) > 0)
        data.C.spike.swrMatrix(spk) = 1;
      end
    end
    
    for swr = 1:data.SWR.spike.nEventsA
      if (sum(coincStatus(data.SWR.spike.evStartA(swr) : data.SWR.spike.evEndA(swr))) > 0)
        data.SWR.spike.evMatrix(swr) = 1;
      end
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
    
    % Calculate event matrix (currently mandatory)
    if param.calcEvMatrixOption
      data.C.burst.swrMatrix  = zeros(data.C.burst.nEventsA, 1);
      data.SWR.burst.evMatrix = zeros(data.SWR.burst.nEventsA, 1);
      
      coincStatus = data.C.burst.evStatusA .* data.SWR.burst.evStatusA;
      
      for bst = 1:data.C.burst.nEventsA
        if (sum(coincStatus(data.C.burst.evStartA(bst) : data.C.burst.evEndA(bst))) > 0)
          data.C.burst.swrMatrix(bst) = 1;
        end
      end
      
      for swr = 1:data.SWR.burst.nEventsA
        if (sum(coincStatus(data.SWR.burst.evStartA(swr) : data.SWR.burst.evEndA(swr))) > 0)
          data.SWR.burst.evMatrix(swr) = 1;
        end
      end
    end
  end
  
  %% Calculate ripple phase of spikes
  if ~isempty(data.SWR.spike.evStartA) && isfield(data.R.SWR.phase, 'evPhase')
    fprintf(['detecting oscillation phase of spikes (file ' dataFileName ')... ']);
    
    % Initialize phase variables
    data.C.SWR.spike.ripplePhase = [];
    data.C.SWR.spike.ripplePhase{length(data.C.SWR.spike.evStartA)} = [];
    data.C.SWR.spike.ripplePhase = data.C.SWR.spike.ripplePhase';
    
    data.C.spike.ripplePhase  = NaN*ones(length(data.C.spike.evStartA), 1);
    data.C.spike.ripplePhaseX = NaN*ones(length(data.C.spike.evPeak), 1);
    data.C.spike.ripplePhaseY = NaN*ones(length(data.C.spike.evPeak), 1);
    
    % Calculate gamma phase if available (have to iterate through structure to not throw error)
    if isfield(data, 'gamma')
      if isfield(data.gamma, 'SWR')
        if isfield(data.gamma.SWR, 'phase')
          data.C.SWR.spike.gammaPhase = [];
          data.C.SWR.spike.gammaPhase{length(data.C.SWR.spike.evStartA)} = [];
          data.C.SWR.spike.gammaPhase = data.C.SWR.spike.gammaPhase';
          
          data.C.spike.gammaPhase  = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.gammaPhaseX = NaN*ones(length(data.C.spike.evStartA), 1);
          data.C.spike.gammaPhaseY = NaN*ones(length(data.C.spike.evStartA), 1);
        end
      end
    end
    
    nSamples = max(cellfun(@length, data.SWR.event));
    
    % Thresholds for phase using previously saved param.sdMult from analyzeLFP
    minAmpR = std(data.R.tSeries) * data.LFP.param.sdMult;
    if isfield(data.C.SWR.spike, 'gammaPhase') minAmpG = std(data.gamma.tSeries) * data.LFP.param.sdMult; end
    
    for spk1 = 1:data.C.spike.nEventsA
      if (data.C.spike.swrMatrix(spk1) == 1) % Coincident spike
        spkInd = data.C.spike.evIndex(:, 2);
        swr    = data.C.spike.evIndex(spkInd == spk1, 1);
        loWin  = max(round(data.SWR.spike.evPeakA(swr) - param.swrWindow / data.LFP.samplingInt), 1);
        spk2   = find(data.C.SWR.spike.evStartA{swr} == data.C.spike.evStartA(spk1) - (loWin - 1));
        
        % Only consider non-truncated SWRs with spikes
        if(length(data.SWR.event{swr}) == nSamples) && ~isempty(spk2)
          swr
          peakTime = data.SWR.evTiming(data.C.SWR.spike.evPeakA{swr}(spk2));
          
          % Only select spikes for which trough-to-peak amplitude is greater than minAmpR
          precMin  = max(data.R.SWR.phase.minLoc{swr}(data.R.SWR.phase.minLoc{swr} < peakTime));
          procMin  = min(data.R.SWR.phase.minLoc{swr}(data.R.SWR.phase.minLoc{swr} >= peakTime));
          precMax  = max(data.R.SWR.phase.maxLoc{swr}(data.R.SWR.phase.maxLoc{swr} < peakTime));
          procMax  = min(data.R.SWR.phase.maxLoc{swr}(data.R.SWR.phase.maxLoc{swr} >= peakTime));
          
          if ~isempty(precMin) && ~isempty(precMax) && ~isempty(procMin) && ~isempty(procMax)
            if precMin > precMax
              ampR = data.R.SWR.phase.maxVal{swr}(data.R.SWR.phase.maxLoc{swr} == procMax) - data.R.SWR.phase.minVal{swr}(data.R.SWR.phase.minLoc{swr} == precMin);
            else
              ampR = data.R.SWR.phase.maxVal{swr}(data.R.SWR.phase.maxLoc{swr} == precMax) - data.R.SWR.phase.minVal{swr}(data.R.SWR.phase.minLoc{swr} == procMin);
            end
            
            if ampR > minAmpR
              data.C.SWR.spike.ripplePhase{swr}(spk2) = data.R.SWR.phase.evPhase{swr}(data.C.SWR.spike.evPeakA{swr}(spk2));
              data.C.spike.ripplePhase(spk1)  = data.C.SWR.spike.ripplePhase{swr}(spk2);
              data.C.spike.ripplePhaseX(spk1) = cos(data.C.spike.ripplePhase(spk1));
              data.C.spike.ripplePhaseY(spk1) = sin(data.C.spike.ripplePhase(spk1));
            end
          end
          
          % Calculate gamma phase if available
          if isfield(data.C.SWR.spike, 'gammaPhase')
            
            % Only select spikes for which trough-to-peak amplitude is greater than minAmpG
            precMin  = max(data.gamma.SWR.phase.minLoc{swr}(data.gamma.SWR.phase.minLoc{swr} < peakTime));
            procMin  = min(data.gamma.SWR.phase.minLoc{swr}(data.gamma.SWR.phase.minLoc{swr} >= peakTime));
            precMax  = max(data.gamma.SWR.phase.maxLoc{swr}(data.gamma.SWR.phase.maxLoc{swr} < peakTime));
            procMax  = min(data.gamma.SWR.phase.maxLoc{swr}(data.gamma.SWR.phase.maxLoc{swr} >= peakTime));
            
            if ~isempty(precMin) && ~isempty(precMax) && ~isempty(procMin) && ~isempty(procMax)
              if precMin > precMax
                ampG = data.gamma.SWR.phase.maxVal{swr}(data.gamma.SWR.phase.maxLoc{swr} == procMax) - data.gamma.SWR.phase.minVal{swr}(data.gamma.SWR.phase.minLoc{swr} == precMin);
              else
                ampG = data.gamma.SWR.phase.maxVal{swr}(data.gamma.SWR.phase.maxLoc{swr} == precMax) - data.gamma.SWR.phase.minVal{swr}(data.gamma.SWR.phase.minLoc{swr} == procMin);
              end
              
              if ampG > minAmpG
                data.C.SWR.spike.gammaPhase{swr}(spk2) = data.gamma.SWR.phase.evPhase{swr}(data.C.SWR.spike.evPeakA{swr}(spk2));
                data.C.spike.gammaPhase(spk1) = data.C.SWR.spike.gammaPhase{swr}(spk2);
                data.C.spike.gammaPhaseX(spk1) = cos(data.C.spike.gammaPhase(spk1));
                data.C.spike.gammaPhaseY(spk1) = sin(data.C.spike.gammaPhase(spk1));
              end
            end
          end
        end
      end
    end
    fprintf('done\n');
  end
end

% Re-order structure arrays:
data.C         = orderStruct(data.C);
data.C.spike   = orderStruct(data.C.spike);
data.C.param   = orderStruct(data.C.param);
data.SWR       = orderStruct(data.SWR);
data.SWR.spike = orderStruct(data.SWR.spike);

if param.parseSpkOption
  data.C.SWR.spike = orderStruct(data.C.SWR.spike);
end

if param.swrBstOption
  data.C.burst   = orderStruct(data.C.burst);
  data.SWR.burst = orderStruct(data.SWR.burst);
end


%% Save matlab file
if all(saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(saveFile,'-struct','data');
  fprintf('done\n');
end

%% Export Spike event table
if all(expSpkFile) && param.expSpkEvOption
  fprintf(['exporting spike event file ' expSpkFileName '... ']);
  exportSpkEvents(data, saveFile, expSpkFile);
  fprintf('done\n');
end

%% Export Burst event table
if all(expBstFile) && param.expBstEvOption
  fprintf(['exporting burst event file ' expBstFileName '... ']);
  exportBstEvents(data, saveFile, expBstFile);
  fprintf('done\n');
end

%% Export SWR event file
if (all(expSWRFile) && param.expSWREvOption)
  fprintf(['exporting updated SWR events (file ' expSWRFileName ')... ']);
  exportSWREvents(data, saveFile, expSWRFile)
  fprintf('done\n');
end

end
