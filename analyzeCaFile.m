function [data, hand] = analyzeCaFile(data, hand, param, saveFile, expCaFile, expSWRFile, expStimFile)
%% [data, hand] = analyzeCaFile(data, hand, param, saveFile, expCaFile, expSWRFile)

%  Script to detect Ca transients above thresholds (previously calculated).
%  Will correlate to previosuly detected SWR events if option selected

%% Handle input arguments - if not entered
if (nargin < 7) expStimFile = []; end
if (nargin < 6) expSWRFile  = []; end
if (nargin < 5) expCaFile   = []; end
if (nargin < 4) saveFile    = []; end
if (nargin < 3) param       = struct; end
if (nargin < 2) hand        = struct; end
if (nargin < 1) data        = struct; end

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

% Check if necessary data structures are present
if ~isfield(data,'Ca')
  [fileName, filePath] = uigetfile('.mat', 'Select *.mat file with imported Ca data');
  dataFile = [filePath fileName];
  if ~all(dataFile) error('No previously analyzed *.mat file selected'); end
  data = load(dataFile);
end

if ~isfield(data,'Ca') error('missing Calcium data, run processCaFile.m first'); end

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

% Select file for stim events (if selected)
if isempty(expStimFile) && param.expStimEvOption
  defaultName = [parentPath dataFileName '_stimEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of stimulation events', defaultName);
  expStimFile = [exportPath exportName];
end

% Set alignment range of Calcium based on param.skipDetectLim
CaRange = find(data.Ca.timing >= 1000 * param.skipDetectLim);

% Set alignment range of LFP based on param.skipDetectLim
if param.swrCaOption || param.spkCaOption || param.stimCaOption 
  lfpRange = find(data.LFP.timing >= 1000 * param.skipDetectLim);
end

%% Calcium Event detection
if param.peakDetectCa
  fprintf('detecting events %4.0f SD above baseline... ', param.sdMult);
  
  % Throw error if thresholds not yet calculated:
  if ~isfield(data.Ca, 'baseMean') || ~isfield(data.Ca, 'baseThresh') || ~isfield(data.Ca, 'peakThresh') 
    error('Missing threholds, first run calcThresh.m');
  end
    
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
  
  %% Align files
  fprintf(['aligning time arrays for SWR-Calcium coincidence analysis (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    if param.useSWRDurationOption
      [data.SWR.Ca.evStatusA, data.SWR.Ca.evStartA, data.SWR.Ca.evEndA, data.Ca.SWR.evStatusA(:,ch), data.Ca.SWR.evStartA{ch}, data.Ca.SWR.evEndA{ch}, data.SWR.Ca.timingA] = ...
        timeAlign(data.SWR.evStatus(lfpRange), data.Ca.evStatus(:,ch), data.LFP.timing(lfpRange), data.Ca.timing(CaRange));
      
    elseif param.useSWRWindowOption
      [data.SWR.Ca.evStatusA, data.SWR.Ca.evStartA, data.SWR.Ca.evEndA, data.Ca.SWR.evStatusA(:,ch), data.Ca.SWR.evStartA{ch}, data.Ca.SWR.evEndA{ch}, data.SWR.Ca.timingA] = ...
        timeAlign(data.SWR.evStatusStand(lfpRange), data.Ca.evStatus(:,ch), data.LFP.timing(lfpRange), data.Ca.timing(CaRange));
      
    end
    
    % Truncate Ca peaks if necessary
    data.Ca.SWR.evPeakA{ch} = data.Ca.evPeak{ch}(1:length(data.Ca.SWR.evStartA{ch}));
  end
  
  % Re-calculate SWR peaks from mid-point of start and end:
  data.SWR.Ca.evPeakA = round(0.5*(data.SWR.Ca.evStartA + data.SWR.Ca.evEndA));
    
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
  
  % Initialize event matrices:
  data.SWR.Ca.evMatrix = zeros(data.SWR.Ca.nEventsA, data.Ca.nChannels);
  data.Ca.SWR.swr   = struct;
  data.Ca.SWR.spont = struct;
  data.Ca.SWR.swr.evMatrix{data.Ca.nChannels}   = [];
  data.Ca.SWR.spont.evMatrix{data.Ca.nChannels} = [];

  % Initialize SWR and spont event structures
  data.Ca.SWR.swr.amp{data.Ca.nChannels} = [];
  data.Ca.SWR.swr.duration{data.Ca.nChannels} = [];
  data.Ca.SWR.swr.ampAve    = zeros(1, data.Ca.nChannels);
  data.Ca.SWR.swr.durAve    = zeros(1, data.Ca.nChannels);
  data.Ca.SWR.swr.frequency = zeros(1, data.Ca.nChannels);
  data.Ca.SWR.swr.nEvents   = zeros(1, data.Ca.nChannels);

  data.Ca.SWR.spont.amp{data.Ca.nChannels} = [];
  data.Ca.SWR.spont.duration{data.Ca.nChannels} = [];
  data.Ca.SWR.spont.ampAve    = zeros(1, data.Ca.nChannels);
  data.Ca.SWR.spont.durAve    = zeros(1, data.Ca.nChannels);
  data.Ca.SWR.spont.frequency = zeros(1, data.Ca.nChannels);
  data.Ca.SWR.spont.nEvents   = zeros(1, data.Ca.nChannels);
  
  % Calculate total time of SWR and spont periods [s]
  if param.useSWRDurationOption
    swrTime = data.SWR.Ca.nEventsA * mean(data.SWR.duration) / 1000;
  else
    swrTime = data.SWR.Ca.nEventsA * param.swrWindow / 1000;
  end
  spontTime = ((data.SWR.Ca.timingA(end) - data.SWR.Ca.timingA(1)) / 1000) - swrTime;

  % Calculate events matrices and SWR/Spont variables
  for ch = 1:data.Ca.nChannels
    
    % Initialize cell arrays
    data.Ca.SWR.swr.evMatrix{ch}   = zeros(data.Ca.SWR.nEventsA(ch), 1);
    data.Ca.SWR.spont.evMatrix{ch} = zeros(data.Ca.SWR.nEventsA(ch), 1);
    
    % Intersection status array:
    evStatusC = data.SWR.Ca.evStatusA .* data.Ca.SWR.evStatusA(:,ch);
    
    % Event matrix, indicating for each SWR whether cells are active:
    for swr = 1:data.SWR.Ca.nEventsA
      if (sum(evStatusC(data.SWR.Ca.evStartA(swr) : data.SWR.Ca.evEndA(swr))) > 0)
        data.SWR.Ca.evMatrix(swr, ch) = 1;
      end
    end

    % Event matrices, indicating for each Ca transient whether it was during SWR or Spont periods:  
    for ev = 1:data.Ca.SWR.nEventsA(ch)
      if (sum(evStatusC(data.Ca.SWR.evStartA{ch}(ev) : data.Ca.SWR.evEndA{ch}(ev))) > 0)
        data.Ca.SWR.swr.evMatrix{ch}(ev)   = 1;
      else
        data.Ca.SWR.spont.evMatrix{ch}(ev) = 1;
      end
    end
    
    % Calculate Ca transient characteristics for SWR events
    data.Ca.SWR.swr.nEvents(ch)   = sum(data.Ca.SWR.swr.evMatrix{ch});
    data.Ca.SWR.swr.amp{ch}       = nonzeros(data.Ca.amp{ch}(1:data.Ca.SWR.nEventsA(ch)) .* data.Ca.SWR.swr.evMatrix{ch}')';
    data.Ca.SWR.swr.duration{ch}  = nonzeros(data.Ca.duration{ch}(1:data.Ca.SWR.nEventsA(ch)) .* data.Ca.SWR.swr.evMatrix{ch}')';
    data.Ca.SWR.swr.ampAve(ch)    = mean(data.Ca.SWR.swr.amp{ch});
    data.Ca.SWR.swr.durAve(ch)    = mean(data.Ca.SWR.swr.duration{ch});
    data.Ca.SWR.swr.frequency(ch) = data.Ca.SWR.swr.nEvents(ch) / swrTime;
    
    % Calculate Ca transient characteristics for spont events
    data.Ca.SWR.spont.nEvents(ch)   = sum(data.Ca.SWR.spont.evMatrix{ch});
    data.Ca.SWR.spont.amp{ch}       = nonzeros(data.Ca.amp{ch}(1:data.Ca.SWR.nEventsA(ch)) .* data.Ca.SWR.spont.evMatrix{ch}')';
    data.Ca.SWR.spont.duration{ch}  = nonzeros(data.Ca.duration{ch}(1:data.Ca.SWR.nEventsA(ch)) .* data.Ca.SWR.spont.evMatrix{ch}')';
    data.Ca.SWR.spont.ampAve(ch)    = mean(data.Ca.SWR.spont.amp{ch});
    data.Ca.SWR.spont.durAve(ch)    = mean(data.Ca.SWR.spont.duration{ch});
    data.Ca.SWR.spont.frequency(ch) = data.Ca.SWR.spont.nEvents(ch) / spontTime;
     
  end
  
  data.SWR.Ca.nCellsC = sum(data.SWR.Ca.evMatrix, 2); % # Cells active for each SWR event
  
  %% Correlation Matrices
  data.SWR.Ca.evMatrixCorr = data.SWR.Ca.evMatrix;
  % Only consider events with >0 active cells
  ev2 = 1;
  for ev1 = 1:length(data.SWR.Ca.nCellsC)
    if data.SWR.Ca.nCellsC(ev1) == 0
      data.SWR.Ca.evMatrixCorr(ev2,:) = [];
    else
      ev2 = ev2 + 1;
    end
  end
    
  % Only compute correlations if sufficient number of cells, otherwise may crash
  if data.Ca.nChannels >= 5
    
    % Calculate correlation matrix between SWR events using Jaccard-Similarity distance
    data.SWR.Ca.corrMatrix = 1 - squareform(pdist(data.SWR.Ca.evMatrixCorr, 'jaccard'));
    data.SWR.Ca.corrMatrix(isnan(data.SWR.Ca.corrMatrix)) = 0; % Replace SWRs with no active cells with zero correlation
    data.SWR.Ca.corrMatrix = triu(data.SWR.Ca.corrMatrix, 1); % Replace diagonal and redundant half with zero
    data.SWR.Ca.corrVector = data.SWR.Ca.corrMatrix(triu(true(size(data.SWR.Ca.corrMatrix)), 1));
    [data.SWR.Ca.cdfF, data.SWR.Ca.cdfX] = ecdf(data.SWR.Ca.corrVector);
    
    % Calculate correlation matrix between cells using Jaccard-Similarity distance
    data.Ca.SWR.corrMatrix = 1 - squareform(pdist(data.SWR.Ca.evMatrixCorr', 'jaccard'));
    data.Ca.SWR.corrMatrix(isnan(data.Ca.SWR.corrMatrix)) = 0; % Replace inactive cells with zero correlation
    data.Ca.SWR.corrMatrix = triu(data.Ca.SWR.corrMatrix, 1); % Replace diagonal and redundant half with zero
    data.Ca.SWR.corrVector = data.Ca.SWR.corrMatrix(triu(true(size(data.Ca.SWR.corrMatrix)), 1));
    [data.Ca.SWR.cdfF, data.Ca.SWR.cdfX] = ecdf(data.Ca.SWR.corrVector);
  
  end
  
%   % Plot Event Matrix:
%   figure
%   imagesc('XData', 1:size(data.SWR.Ca.evMatrixCorr,1), 'YData', 1:size(data.SWR.Ca.evMatrixCorr,2), 'CData', data.SWR.Ca.evMatrixCorr');
%   axis([1 size(data.SWR.Ca.evMatrixCorr,1) 1 size(data.SWR.Ca.evMatrixCorr,2)]);
%   caxis([0 1]);
%   colormap(flipud(hot));
%   
%   % Plot Correlation Matrix:
%   figure
%   imagesc('XData', 1:size(data.SWR.Ca.corrMatrix,1), 'YData', 1:size(data.SWR.Ca.corrMatrix,2), 'CData', data.SWR.Ca.corrMatrix);
%   axis([1 size(data.SWR.Ca.corrMatrix,1) 1 size(data.SWR.Ca.corrMatrix,2)]);
%   caxis([0 1]);
%   colormap(flipud(hot));
%   colorbar
  
  % Re-order structure arrays
  data.Ca.SWR = orderfields(data.Ca.SWR);
  data.SWR.Ca = orderfields(data.SWR.Ca);
  data.SWR    = orderfields(data.SWR);
  
end

%% Calculate stim response
if param.stimCaOption
  
  % Create new data structures if not already present
  if ~isfield(data.Ca,'stim') data.Ca.stim = struct; end
  if ~isfield(data.stim,'Ca') data.stim.Ca = struct; end
  
  % Initialize cell arrays:
  data.Ca.stim.evStatusA = [];
  data.Ca.stim.evStartA  = [];
  data.Ca.stim.evPeakA   = [];
  data.Ca.stim.evEndA    = [];
  
  data.Ca.stim.evStartA{1, data.Ca.nChannels} = [];
  data.Ca.stim.evPeakA{1, data.Ca.nChannels}  = [];
  data.Ca.stim.evEndA{1, data.Ca.nChannels}   = [];
  
  data.stim.Ca.evStatusA = [];
  data.stim.Ca.evStartA  = [];
  data.stim.Ca.evEndA    = [];
  
  % Calculate extended stim evStatus
  data.stim.evStatusExt = zeros(length(data.LFP.timing), 1);
  for ev = 1:length(data.stim.evStart)
    loWin = data.stim.evStart(ev) + round(param.stimCaLim1 / data.LFP.samplingInt);
    hiWin = data.stim.evStart(ev) + round(param.stimCaLim2 / data.LFP.samplingInt);
    data.stim.evStatusExt(loWin : hiWin) = 1;
  end
  
  %% Align files
  fprintf(['aligning time arrays for stim-Calcium coincidence analysis (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    [data.stim.Ca.evStatusA, data.stim.Ca.evStartA, data.stim.Ca.evEndA, data.Ca.stim.evStatusA(:,ch), data.Ca.stim.evStartA{ch}, data.Ca.stim.evEndA{ch}, data.stim.Ca.timingA] = ...
      timeAlign(data.stim.evStatusExt(lfpRange), data.Ca.evStatus(:,ch), data.LFP.timing(lfpRange), data.Ca.timing(CaRange));
    
    % Re-calculate Ca peaks - truncating if necessary:
    data.Ca.stim.evPeakA{ch} = data.Ca.evPeak{ch}(1:length(data.Ca.stim.evStartA{ch}));
    if data.Ca.stim.evPeakA{ch}(end) > length(data.stim.Ca.timingA) data.Ca.stim.evPeakA{ch}(end) = length(data.stim.Ca.timingA); end
    
  end
  fprintf('done\n');
  
  %% Calculate stim response
  data.Ca.stim.evPeakStim = zeros(length(data.stim.Ca.evStartA), data.Ca.nChannels);
  data.Ca.stim.evAreaStim = zeros(length(data.stim.Ca.evStartA), data.Ca.nChannels);
  
  data.Ca.stim.tSeriesA = data.Ca.tSeries(CaRange,:);
  data.Ca.stim.tSeriesA = data.Ca.stim.tSeriesA(1:length(data.stim.Ca.timingA),:);
    
  for ch = 1:data.Ca.nChannels
    for ev = 1:length(data.stim.Ca.evStartA)
      data.Ca.stim.evPeakStim(ev, ch) = max(data.Ca.stim.tSeriesA(data.stim.Ca.evStartA(ev) : data.stim.Ca.evEndA(ev), ch));
      data.Ca.stim.evAreaStim(ev, ch) = sum(data.Ca.stim.tSeriesA(data.stim.Ca.evStartA(ev) : data.stim.Ca.evEndA(ev), ch)) * data.Ca.samplingInt;
    end
  end
  
  %% Calculate overlap of stim and Ca events
  
  % Initialize cell arrays:
  data.Ca.stim.evStatusC = [];
  data.Ca.stim.evStartC  = [];
  data.Ca.stim.evEndC    = [];
  data.Ca.stim.evIndex   = [];
  
  data.Ca.stim.evStartC{1, data.Ca.nChannels} = [];
  data.Ca.stim.evEndC{1, data.Ca.nChannels}   = [];
  data.Ca.stim.evIndex{1, data.Ca.nChannels}  = [];
  
  data.stim.Ca.evStatusC = [];
  data.stim.Ca.evStartC  = [];
  data.stim.Ca.evEndC    = [];
  data.stim.Ca.evIndex   = [];
  
  data.stim.Ca.evStartC{1, data.Ca.nChannels} = [];
  data.stim.Ca.evEndC{1, data.Ca.nChannels}   = [];
  data.stim.Ca.evIndex{1, data.Ca.nChannels}  = [];
  
  % Calculate overlap of stim and Ca events for each cell
  fprintf(['detecting Ca transients coincident with stimulation (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    [data.Ca.stim.evStatusC(:,ch), data.Ca.stim.evStartC{ch}, data.Ca.stim.evEndC{ch}, data.Ca.stim.evIndex{ch}] = eventOverlap(data.stim.Ca.evStatusA, data.stim.Ca.evStartA, data.stim.Ca.evEndA, ...
      data.Ca.stim.evStatusA(:,ch), data.Ca.stim.evStartA{ch}, data.Ca.stim.evEndA{ch}, data.stim.Ca.timingA, 2);
  end
  fprintf('done\n');
  
  % Find stims that coincide with at least one Ca transient
  fprintf(['detecting stimulation events coincident with Ca transients (file ' dataFileName ')... ']);
  for ch = 1:data.Ca.nChannels
    [data.stim.Ca.evStatusC(:,ch), data.stim.Ca.evStartC{ch}, data.stim.Ca.evEndC{ch}, data.stim.Ca.evIndex{ch}] = eventOverlap(data.stim.Ca.evStatusA, data.stim.Ca.evStartA, data.stim.Ca.evEndA, ...
      data.Ca.stim.evStatusA(:,ch), data.Ca.stim.evStartA{ch}, data.Ca.stim.evEndA{ch}, data.stim.Ca.timingA, 1);
  end
  fprintf('done\n');
  
  % Calculate summed status of coincident stims over all cells (value signifies # active cells):
  data.stim.Ca.evStatusSumC = sum(data.stim.Ca.evStatusC, 2); % Summed status of coincident stims over all cells
  
  % Parse summed status to get start and end times:
  [data.stim.Ca.evStartSumC, data.stim.Ca.evEndSumC] = eventParse(data.stim.Ca.evStatusSumC);
  
  % Compute total number of events:
  data.stim.Ca.nEventsA      = size(data.stim.Ca.evStartA, 1);      % # stims
  [data.stim.Ca.nEventsC, ~] = cellfun(@size, data.stim.Ca.evStartC); % for each cell, # stims with coincident Ca event
  data.stim.Ca.nEventsSumC   = size(data.stim.Ca.evStartSumC, 1);     % # stims with any coincident cell event
  [data.Ca.stim.nEventsA, ~] = cellfun(@size, data.Ca.stim.evStartA); % for each cell, # Ca events
  data.Ca.stim.nEventsSumA   = sum(data.Ca.stim.nEventsA);            % Sum of all Ca events
  [data.Ca.stim.nEventsC, ~] = cellfun(@size, data.Ca.stim.evStartC); % for each cell, # Ca events with coincident stim
  data.Ca.stim.nEventsSumC   = sum(data.Ca.stim.nEventsC);            % Sum of all Ca events with coincident stim
  data.Ca.stim.fracEventsC   = data.Ca.stim.nEventsC ./ data.Ca.stim.nEventsA; % Fraction of coicident events
  
  % Initialize event matrices:
  data.stim.Ca.evMatrix = zeros(data.stim.Ca.nEventsA, data.Ca.nChannels);
  data.Ca.stim.stim  = struct;
  data.Ca.stim.spont = struct;
  data.Ca.stim.stim.evMatrix{data.Ca.nChannels}  = [];
  data.Ca.stim.spont.evMatrix{data.Ca.nChannels} = [];

  % Initialize stim and spont event structures
  data.Ca.stim.stim.amp{data.Ca.nChannels} = [];
  data.Ca.stim.stim.duration{data.Ca.nChannels} = [];
  data.Ca.stim.stim.ampAve    = zeros(1, data.Ca.nChannels);
  data.Ca.stim.stim.durAve    = zeros(1, data.Ca.nChannels);
  data.Ca.stim.stim.frequency = zeros(1, data.Ca.nChannels);
  data.Ca.stim.stim.nEvents   = zeros(1, data.Ca.nChannels);

  data.Ca.stim.spont.amp{data.Ca.nChannels} = [];
  data.Ca.stim.spont.duration{data.Ca.nChannels} = [];
  data.Ca.stim.spont.ampAve    = zeros(1, data.Ca.nChannels);
  data.Ca.stim.spont.durAve    = zeros(1, data.Ca.nChannels);
  data.Ca.stim.spont.frequency = zeros(1, data.Ca.nChannels);
  data.Ca.stim.spont.nEvents   = zeros(1, data.Ca.nChannels);
  
  % Calculate total time of stim and spont periods [s]
  stimTime  = data.stim.Ca.nEventsA * (param.stimCaLim2 - param.stimCaLim1) / 1000;
  spontTime = ((data.stim.Ca.timingA(end) - data.stim.Ca.timingA(1)) / 1000) - stimTime;

  % Calculate events matrices and stim/Spont variables
  for ch = 1:data.Ca.nChannels
    
    % Initialize cell arrays
    data.Ca.stim.stim.evMatrix{ch}  = zeros(data.Ca.stim.nEventsA(ch), 1);
    data.Ca.stim.spont.evMatrix{ch} = zeros(data.Ca.stim.nEventsA(ch), 1);
    
    % Intersection status array:
    evStatusC = data.stim.Ca.evStatusA .* data.Ca.stim.evStatusA(:,ch);
    
    % Event matrix, indicating for each stim whether cells are active:
    for ev = 1:data.stim.Ca.nEventsA
      if (sum(evStatusC(data.stim.Ca.evStartA(ev) : data.stim.Ca.evEndA(ev))) > 0)
        data.stim.Ca.evMatrix(ev, ch) = 1;
      end
    end

    % Event matrices, indicating for each Ca transient whether it was during SWR or Spont periods:  
    for ev = 1:data.Ca.stim.nEventsA(ch)
      if (sum(evStatusC(data.Ca.stim.evStartA{ch}(ev) : data.Ca.stim.evEndA{ch}(ev))) > 0)
        data.Ca.stim.stim.evMatrix{ch}(ev)  = 1;
      else
        data.Ca.stim.spont.evMatrix{ch}(ev) = 1;
      end
    end
    
    % Calculate Ca transient characteristics for SWR events
    data.Ca.stim.stim.nEvents(ch)   = sum(data.Ca.stim.stim.evMatrix{ch});
    data.Ca.stim.stim.amp{ch}       = nonzeros(data.Ca.amp{ch}(1:data.Ca.stim.nEventsA(ch)) .* data.Ca.stim.stim.evMatrix{ch}')';
    data.Ca.stim.stim.duration{ch}  = nonzeros(data.Ca.duration{ch}(1:data.Ca.stim.nEventsA(ch)) .* data.Ca.stim.stim.evMatrix{ch}')';
    data.Ca.stim.stim.ampAve(ch)    = mean(data.Ca.stim.stim.amp{ch});
    data.Ca.stim.stim.durAve(ch)    = mean(data.Ca.stim.stim.duration{ch});
    data.Ca.stim.stim.frequency(ch) = data.Ca.stim.stim.nEvents(ch) / stimTime;
    
    % Calculate Ca transient characteristics for spont events
    data.Ca.stim.spont.nEvents(ch)   = sum(data.Ca.stim.spont.evMatrix{ch});
    data.Ca.stim.spont.amp{ch}       = nonzeros(data.Ca.amp{ch}(1:data.Ca.stim.nEventsA(ch)) .* data.Ca.stim.spont.evMatrix{ch}')';
    data.Ca.stim.spont.duration{ch}  = nonzeros(data.Ca.duration{ch}(1:data.Ca.stim.nEventsA(ch)) .* data.Ca.stim.spont.evMatrix{ch}')';
    data.Ca.stim.spont.ampAve(ch)    = mean(data.Ca.stim.spont.amp{ch});
    data.Ca.stim.spont.durAve(ch)    = mean(data.Ca.stim.spont.duration{ch});
    data.Ca.stim.spont.frequency(ch) = data.Ca.stim.spont.nEvents(ch) / spontTime;
     
  end

  data.stim.Ca.nCellsC = sum(data.stim.Ca.evMatrix, 2); % # Cells active for each stim event
  
  %% Correlation Matrices
  data.stim.Ca.evMatrixCorr = data.stim.Ca.evMatrix;
  % Only consider events with >0 active cells
  ev2 = 1;
  for ev1 = 1:length(data.stim.Ca.nCellsC)
    if data.stim.Ca.nCellsC(ev1) == 0
      data.stim.Ca.evMatrixCorr(ev2,:) = [];
    else
      ev2 = ev2 + 1;
    end
  end
  
  % Only compute correlations if sufficient number of cells, otherwise may crash
  if data.Ca.nChannels >= 5
    
    % Calculate correlation matrix between stim events using Jaccard-Similarity distance
    data.stim.Ca.corrMatrix = 1 - squareform(pdist(data.stim.Ca.evMatrixCorr, 'jaccard'));
    data.stim.Ca.corrMatrix(isnan(data.stim.Ca.corrMatrix)) = 0; % Replace stims with no active cells with zero correlation
    data.stim.Ca.corrMatrix = triu(data.stim.Ca.corrMatrix, 1); % Replace diagonal and redundant half with zero
    data.stim.Ca.corrVector = data.stim.Ca.corrMatrix(triu(true(size(data.stim.Ca.corrMatrix)), 1));
    [data.stim.Ca.cdfF, data.stim.Ca.cdfX] = ecdf(data.stim.Ca.corrVector);
    
    % Calculate correlation matrix between cells using Jaccard-Similarity distance
    data.Ca.stim.corrMatrix = 1 - squareform(pdist(data.stim.Ca.evMatrixCorr', 'jaccard'));
    data.Ca.stim.corrMatrix(isnan(data.Ca.stim.corrMatrix)) = 0; % Replace inactive cells with zero correlation
    data.Ca.stim.corrMatrix = triu(data.Ca.stim.corrMatrix, 1); % Replace diagonal and redundant half with zero
    data.Ca.stim.corrVector = data.Ca.stim.corrMatrix(triu(true(size(data.Ca.stim.corrMatrix)), 1));
    [data.Ca.stim.cdfF, data.Ca.stim.cdfX] = ecdf(data.Ca.stim.corrVector);
    
  end
  
  % Re-order structure arrays
  data.Ca.stim = orderfields(data.Ca.stim);
  data.stim.Ca = orderfields(data.stim.Ca);
  data.stim    = orderfields(data.stim);

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

%% Export stim event file
if (all(expStimFile) && param.expStimEvOption)
  fprintf(['exporting stimulation events (file ' dataFileName ')... ']);
  exportStimEvents(data, saveFile, expStimFile)
  fprintf('done\n');
end

end