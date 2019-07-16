function [data, hand] = analyzeCaFile(data, hand, param, saveFile, expCaFile, expSWRFile)
%% [data, hand] = analyzeCaFile(data, hand, param, saveFile, expCaFile, expSWRFile)

%% Handle input arguments - if not entered
if (nargin < 6) expSWRFile = []; end
if (nargin < 5) expCaFile  = []; end
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
if ~isfield(param,'sdBaseFactor')         param.sdBaseFactor         = 1;   end
if ~isfield(param,'skipDetectLim')        param.skipDetectLim        = 2;   end
if ~isfield(param,'consThreshOption')     param.consThreshOption     = 0;   end
if ~isfield(param,'calcCaWindowOption')   param.calcCaWindowOption   = 0;   end
if ~isfield(param,'limPeakDetectOption')  param.limPeakDetectOption  = 0;   end
if ~isfield(param,'swrCaOption')          param.swrCaOption          = 0;   end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 1;   end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption   = 0;   end
if ~isfield(param,'swrWindow')            param.swrWindow            = 100; end
if ~isfield(param,'expCaEvOption')        param.expCaEvOption        = 1;   end
if ~isfield(param,'expSWREvOption')       param.expSWREvOption       = 0;   end
if ~isfield(param,'spkCaOption')          param.spkCaOption          = 0;   end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption      = 0;   end

% Check if necessary data structures are present
if ~isfield(data,'Ca')
  [fileName, filePath] = uigetfile('.mat', 'Select *.mat file with imported Ca data');
  dataFile = [filePath fileName];
  if ~all(dataFile) error('No previously analyzed *.mat file selected'); end
  data = load(dataFile);
end

if ~isfield(data,'Ca') error('missing Calcium data, run importCaFile.m first'); end

if param.swrCaOption
  if ~isfield(data,'SWR') error('Must analyze LFP channel for SWR events before proceeding'); end
end

if param.spkCaOption
  if ~isfield(data,'C') error('Must analyze cell channel before proceeding'); end
end

% If not supplied, prompt for save file
if isempty(saveFile)
  if isfield(data, 'saveFile')
    saveFile = data.saveFile;
  else
    error('Missing file information, problem with upstream scripts');
  end
end
[parentPath, dataFileName, saveExt] = parsePath(saveFile);
data.saveFile = saveFile;
data.saveName = [dataFileName '.' saveExt];

% Select export file for Calcium events (if selected)
if isempty(expCaFile) && param.expCaEvOption
  defaultName = [parentPath dataFileName '_CaEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of Calcium events', defaultName);
  expCaFile = [exportPath exportName];
end
[parentPath, ~, ~] = parsePath(saveFile);

% Select re-export file for SWR events (if selected)
if isempty(expSWRFile) && param.expSWREvOption
  defaultName = [parentPath dataFileName '_swrEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export updated table of SWR events', defaultName);
  expSWRFile = [exportPath exportName];
end


%% Calcium Event detection
if param.peakDetectCa
  fprintf('detecting events %4.0f SD above baseline... ', param.sdMult);
  
  % Throw error if thresholds not yet calculated:
  if ~isfield(data.Ca, 'baseMean') || ~isfield(data.Ca, 'baseThresh') || ~isfield(data.Ca, 'peakThresh') 
    error('Missing threholds, first run calcThresh.m');
  end

  % Set detection range based on param.skipDetectLim
  CaRange = find(data.Ca.timing >= 1000 * param.skipDetectLim);
  
  % Initialize cell arrays:
  data.Ca.evStatus  = [];
  data.Ca.evStart   = [];
  data.Ca.evPeak    = [];
  data.Ca.evEnd     = [];
  data.Ca.frequency = [];
  data.Ca.IEI       = [];
  data.Ca.duration  = [];
  data.Ca.amp       = [];
  
  data.Ca.evStatus = zeros(length(data.Ca.timing(CaRange)), data.Ca.nChannels);
  data.Ca.evStart{1, data.Ca.nChannels}  = [];
  data.Ca.evPeak{1, data.Ca.nChannels}   = [];
  data.Ca.evEnd{1, data.Ca.nChannels}    = [];
  data.Ca.IEI{1, data.Ca.nChannels}      = [];
  data.Ca.duration{1, data.Ca.nChannels} = [];
  data.Ca.amp{1, data.Ca.nChannels}      = [];
  data.Ca.frequency = zeros(1, data.Ca.nChannels);
  data.Ca.ampAve    = NaN * zeros(1, data.Ca.nChannels);
  data.Ca.durAve    = NaN * zeros(1, data.Ca.nChannels);
  data.Ca.IEIAve    = NaN * zeros(1, data.Ca.nChannels);
  data.Ca.nEvents   = NaN * zeros(1, data.Ca.nChannels);

  % Detect peaks based on previously identified thresholds
  for ch = 1:data.Ca.nChannels
    timingWin  = data.Ca.timing(CaRange);
    tSeriesWin = data.Ca.tSeries(CaRange,ch);
    
    warning ('off','all');
    [data.Ca.evStatus(:,ch), data.Ca.evStart{ch}, data.Ca.evPeak{ch}, data.Ca.evEnd{ch}] = ...
      peakFindUnique(tSeriesWin, timingWin, data.Ca.peakThresh(ch), data.Ca.baseThresh(ch), 1);
    warning ('on','all');
    
    if ~isnan(data.Ca.evStart{ch})
      for i = 1:length(data.Ca.evStart{ch})
        data.Ca.duration{ch}(i) = (timingWin(data.Ca.evEnd{ch}(i)) - timingWin(data.Ca.evStart{ch}(i)));
        data.Ca.amp{ch}(i)      = tSeriesWin(data.Ca.evPeak{ch}(i)) - data.Ca.baseMean(ch);
        if (i > 1) data.Ca.IEI{ch} = horzcat(data.Ca.IEI{ch}, (timingWin(data.Ca.evPeak{ch}(i)) - timingWin(data.Ca.evPeak{ch}(i-1))) / 1000); end
      end
      data.Ca.nEvents(ch)   = length(data.Ca.evStart{ch});
      data.Ca.frequency(ch) = length(data.Ca.evStart{ch}) / ((timingWin(length(timingWin)) - timingWin(1)) / 1000);
      data.Ca.ampAve(ch)    = mean(data.Ca.amp{ch});
      data.Ca.durAve(ch)    = mean(data.Ca.duration{ch});
      data.Ca.IEIAve(ch)    = mean(data.Ca.IEI{ch});
    end
  end
  fprintf('done\n');
end

%% Calculate standard SWR status (if selected)
if param.useSWRWindowOption && (param.swrCaOption || param.spkCaOption)
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

%% Correlate SWR and Ca events
if param.swrCaOption
  
  % Create new data structures if not already present
  if ~isfield(data.Ca,'SWR') data.Ca.SWR = struct; end
  if ~isfield(data.SWR,'Ca') data.SWR.Ca = struct; end
  
  % Initialize cell arrays:
  data.Ca.SWR.evStatusA = [];
  data.Ca.SWR.evStartA  = [];
  data.Ca.SWR.evPeakA   = [];
  data.Ca.SWR.evEndA    = [];
  
  data.Ca.SWR.evStartA{1, data.Ca.nChannels} = [];
  data.Ca.SWR.evPeakA{1, data.Ca.nChannels}  = [];
  data.Ca.SWR.evEndA{1, data.Ca.nChannels}   = [];
  
  data.SWR.Ca.evStatusA = [];
  data.SWR.Ca.evStartA  = [];
  data.SWR.Ca.evPeakA   = [];
  data.SWR.Ca.evEndA    = [];
  
  % Set alignment range of LFP based on param.skipDetectLim
  lfpRange = find(data.LFP.timing >= 1000 * param.skipDetectLim);
  
  %% Align files
  fprintf(['aligning time arrays for SWR-Calcium coincidence analysis (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    if param.useSWRDurationOption
      
      [data.SWR.Ca.evStatusA, data.SWR.Ca.evStartA, data.SWR.Ca.evEndA, data.Ca.SWR.evStatusA(:,ch), data.Ca.SWR.evStartA{ch}, data.Ca.SWR.evEndA{ch}, data.SWR.Ca.timingA] = ...
        timeAlign(data.SWR.evStatus(lfpRange), data.Ca.evStatus(:,ch), data.LFP.timing(lfpRange), data.Ca.timing(CaRange));
      % Re-calculate peak - truncating if necessary:
      data.SWR.Ca.evPeakA = data.SWR.evPeak(1:length(data.SWR.Ca.evStartA));
      
    elseif param.useSWRWindowOption
      
      [data.SWR.Ca.evStatusA, data.SWR.Ca.evStartA, data.SWR.Ca.evEndA, data.Ca.SWR.evStatusA(:,ch), data.Ca.SWR.evStartA{ch}, data.Ca.SWR.evEndA{ch}, data.SWR.Ca.timingA] = ...
        timeAlign(data.SWR.evStatusStand(lfpRange), data.Ca.evStatus(:,ch), data.LFP.timing(lfpRange), data.Ca.timing(CaRange));
      % Re-calculate SWR peak from mid-point of start and end:
      data.SWR.Ca.evPeakA = round(0.5*(data.SWR.Ca.evStartA + data.SWR.Ca.evEndA));
      
    end
    data.Ca.SWR.evPeakA{ch} = data.Ca.evPeak{ch}(1:length(data.Ca.SWR.evStartA{ch}));
  end
  fprintf('done\n');
  
  %% Calculate overlap of events
  
  % Initialize cell arrays:
  data.Ca.SWR.evStatusC = [];
  data.Ca.SWR.evStartC  = [];
  data.Ca.SWR.evEndC    = [];
  data.Ca.SWR.evIndex   = [];
  
  data.Ca.SWR.evStartC{1, data.Ca.nChannels} = [];
  data.Ca.SWR.evEndC{1, data.Ca.nChannels}   = [];
  data.Ca.SWR.evIndex{1, data.Ca.nChannels}  = [];
  
  data.SWR.Ca.evStatusC = [];
  data.SWR.Ca.evStartC  = [];
  data.SWR.Ca.evEndC    = [];
  data.SWR.Ca.evIndex   = [];
  
  data.SWR.Ca.evStartC{1, data.Ca.nChannels} = [];
  data.SWR.Ca.evEndC{1, data.Ca.nChannels}   = [];
  data.SWR.Ca.evIndex{1, data.Ca.nChannels}  = [];
  
  % Calculate overlap of SWR and Ca events for each cell
  fprintf(['detecting Ca transients coincident with SWRs (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    [data.Ca.SWR.evStatusC(:,ch), data.Ca.SWR.evStartC{ch}, data.Ca.SWR.evEndC{ch}, data.Ca.SWR.evIndex{ch}] = eventOverlap(data.SWR.Ca.evStatusA, data.SWR.Ca.evStartA, data.SWR.Ca.evEndA, ...
      data.Ca.SWR.evStatusA(:,ch), data.Ca.SWR.evStartA{ch}, data.Ca.SWR.evEndA{ch}, data.SWR.Ca.timingA, 2);
  end
  fprintf('done\n');
  
  % Find SWRs that coincide with at least one Ca transient
  fprintf(['detecting SWRs coincident with Ca transients (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    [data.SWR.Ca.evStatusC(:,ch), data.SWR.Ca.evStartC{ch}, data.SWR.Ca.evEndC{ch}, data.SWR.Ca.evIndex{ch}] = eventOverlap(data.SWR.Ca.evStatusA, data.SWR.Ca.evStartA, data.SWR.Ca.evEndA, ...
      data.Ca.SWR.evStatusA(:,ch), data.Ca.SWR.evStartA{ch}, data.Ca.SWR.evEndA{ch}, data.SWR.Ca.timingA, 1);
  end
  fprintf('done\n');
  
  % Calculate summed status of coincident SWRs over all cells (value signifies # active cells):
  data.SWR.Ca.evStatusSumC = sum(data.SWR.Ca.evStatusC, 2); % Summed status of coincident SWRs over all cells
  
  % Parse summed status to get start and end times:
  [data.SWR.Ca.evStartSumC, data.SWR.Ca.evEndSumC] = eventParse(data.SWR.Ca.evStatusSumC);
  
  % Compute total number of events:
  data.SWR.Ca.nEventsA      = size(data.SWR.Ca.evStartA, 1);        % # SWRs
  [data.SWR.Ca.nEventsC, ~] = cellfun(@size, data.SWR.Ca.evStartC); % for each cell, # SWRs with coincident Ca event
  data.SWR.Ca.nEventsSumC   = size(data.SWR.Ca.evStartSumC, 1);     % # SWRs with any coincident cell event
  [data.Ca.SWR.nEventsA, ~] = cellfun(@size, data.Ca.SWR.evStartA); % for each cell, # Ca events
  data.Ca.SWR.nEventsSumA   = sum(data.Ca.SWR.nEventsA);            % Sum of all Ca events
  [data.Ca.SWR.nEventsC, ~] = cellfun(@size, data.Ca.SWR.evStartC); % for each cell, # Ca events with coincident SWR
  data.Ca.SWR.nEventsSumC   = sum(data.Ca.SWR.nEventsC);            % Sum of all Ca events with coincident SWR
  data.Ca.SWR.fracEventsC   = data.Ca.SWR.nEventsC ./ data.Ca.SWR.nEventsA; % Fraction of coicident events
  
  % Calculate simplified event matrix:
  data.SWR.Ca.evMatrix = zeros(data.SWR.Ca.nEventsA, data.Ca.nChannels);
  
  for ch = 1:data.Ca.nChannels
    evStatusC = data.SWR.Ca.evStatusA .* data.Ca.SWR.evStatusA(:,ch);
    
    for ev = 1:data.SWR.Ca.nEventsA
      if (sum(evStatusC(data.SWR.Ca.evStartA(ev) : data.SWR.Ca.evEndA(ev))) > 0)
        data.SWR.Ca.evMatrix(ev, ch) = 1;
      end
    end
  end
  
  data.SWR.Ca.nCellsC = sum(data.SWR.Ca.evMatrix, 2); % # Cells active for each SWR event
  
  % Re-order structure arrays
  data.Ca.SWR = orderfields(data.Ca.SWR);
  data.SWR.Ca = orderfields(data.SWR.Ca);
  data.SWR    = orderfields(data.SWR);
  
end

data       = orderfields(data);
data.param = orderfields(data.param);
data.Ca    = orderfields(data.Ca);
data.Ca.param = orderfields(data.Ca.param);

%% Save file
if all(saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(saveFile,'-struct','data');
  fprintf('done\n');
end

%% Export Calcium event file
if (all(expCaFile) && param.expCaEvOption)
  fprintf(['exporting Ca events (file ' dataFileName ')... ']);
  exportCaEvents(data, saveFile, expCaFile)
  fprintf('done\n');
end

%% Export SWR event file
if (all(expSWRFile) && param.expSWREvOption)
  fprintf(['exporting SWR events (file ' dataFileName ')... ']);
  exportSWREvents(data, saveFile, expSWRFile)
  fprintf('done\n');
end

end