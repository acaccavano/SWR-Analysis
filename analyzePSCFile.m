function [data, hand] = analyzePSCFile(data, hand, param, saveFile, pscFile, expPSCFile, expSWRFile)
%% [data, hand] = analyzePSCFile(data, hand, param, saveFile, pscFile, expPSCFile, expSWRFile)

%  function to correlate detected PSCs (both IPSCs and EPSCs) with previosuly analyzed SWR events

%% Handle input arguments - if not entered
if (nargin < 7) expSWRFile = []; end
if (nargin < 6) expPSCFile = []; end
if (nargin < 5) pscFile    = []; end
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
if ~isfield(param,'pscEventPolarity')     param.pscEventPolarity     = 0;   end
if ~isfield(param,'swrPSQOption')         param.swrPSQOption         = 1;   end
if ~isfield(param,'importPSCOption')      param.importPSCOption      = 1;   end
if ~isfield(param,'swrPSCOption')         param.swrPSCOption         = 1;   end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 0;   end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption   = 1;   end
if ~isfield(param,'swrWindow')            param.swrWindow            = 100; end
if ~isfield(param,'parsePSCOption')       param.parsePSCOption       = 1;   end
if ~isfield(param,'calcEvMatrixOption')   param.calcEvMatrixOption   = 1;   end
if ~isfield(param,'expPSCEvOption')       param.expPSCEvOption       = 1;   end
if ~isfield(param,'expSWREvOption')       param.expSWREvOption       = 1;   end
if ~isfield(param,'gammaOption')          param.gammaOption          = 1;   end
if ~isfield(param,'gammaLim1')            param.gammaLim1            = 20;  end
if ~isfield(param,'gammaLim2')            param.gammaLim2            = 50;  end
if ~isfield(param,'rOption')              param.rOption              = 1;   end
if ~isfield(param,'rLim1')                param.rLim1                = 120; end
if ~isfield(param,'rLim2')                param.rLim2                = 220; end
if ~isfield(param,'spectOption')          param.spectOption          = 1;   end
if ~isfield(param,'spectLim1')            param.spectLim1            = 1;   end
if ~isfield(param,'spectLim2')            param.spectLim2            = 600; end
if ~isfield(param,'reAnalyzeOption')      param.reAnalyzeOption      = 0;   end
if ~isfield(param,'nBins')                param.nBins                = 100; end % For PSC peak histogram, not currently selectable from UI

if ~isfield(data, 'LFP')
  [fileName, filePath] = uigetfile('.mat', 'Select *.mat file of analyzed LFP + imported cell channel');
  dataFile = [filePath fileName];
  if ~all(dataFile) error('No file data file selected'); end
  data = load(dataFile);
end

% Check if necessary data structures are present
if ~isfield(data,'LFP') error('Must import LFP channel to data structure before proceeding'); end
if ~isfield(data,'C')   error('Must import cell channel to data structure before proceeding'); end
if ~isfield(data.C,'PSC') data.C.spike = struct; end

% Fix in case Fs not calculated correctly:
data.LFP.param.Fs = round(1000 / data.LFP.samplingInt);

if param.swrPSCOption
  if ~isfield(data,'SWR') error('Must analyze SWRs before proceeding'); end
  if ~isfield(data.SWR,'PSC') data.SWR.PSC = struct; end
  if ~isfield(data.C,'SWR') data.C.SWR = struct; end
  
  if param.parsePSCOption
    if ~isfield(data.C.SWR,'PSC') data.C.SWR.PSC = struct; end
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
if param.importPSCOption %% && ~param.reAnalyzeOption %% Commented out want to have option to re-import events but may have unintended consequences
  if isempty(pscFile)
    [fileName, filePath] = uigetfile('.txt', 'Select PSC event *.txt file exported from pClamp', parentPath);
    pscFile = [filePath fileName];
    if ~all(pscFile) error('No PSC file selected'); end
    data.C.PSC.evFile = pscFile;
  end
  [parentPath, pscFileName, ~] = parsePath(pscFile);
end

% Select export PSC event file, if option selected
if param.expPSCEvOption
  if isempty(expPSCFile)
    defaultPath = [parentPath dataFileName '_pscEvents.csv'];
    [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of PSC events', defaultPath);
    expPSCFile = [exportPath exportName];
    if ~all(expPSCFile)
      warning('No PSC events to be exported - no file selected');
    else
      [~, expPSCFileName, ~] = parsePath(expPSCFile);
    end
  else
    [~, expPSCFileName, ~] = parsePath(expPSCFile);
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

%% Calculate peri-SWR areas (charge or PSQ)
if param.swrPSQOption
  fprintf(['integrating cell charge over all SWR events (file ' dataFileName ')... ']);
  
  % Initialize data arrays:
  if ~isfield(data.C,'SWR') data.C.SWR = struct; end
  data.C.SWR.event    = [];
  data.C.SWR.evNorm   = [];
  data.C.SWR.evArea   = [];
  data.C.SWR.baseline = [];
  data.C.SWR.area     = [];
  data.C.SWR.areaQ1   = [];
  data.C.SWR.areaQ2   = [];
  data.C.SWR.areaQ3   = [];
  data.C.SWR.areaQ4   = [];
  
  if ~isnan(data.SWR.evStart)
    
    data.C.SWR.event{length(data.SWR.evStart)}  = [];
    data.C.SWR.evNorm{length(data.SWR.evStart)} = [];
    data.C.SWR.evArea{length(data.SWR.evStart)} = [];
    
    for i = 1:length(data.SWR.evStart)
      
      % Calculate SWR windows:
      loWin = max(round(data.SWR.evPeak(i) - param.swrWindow / data.LFP.samplingInt), 1);
      hiWin = min(round(data.SWR.evPeak(i) + param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
      
      loBaseWin = max(round(data.SWR.evPeak(i) - 0.5 * param.swrWindow / data.LFP.samplingInt), 1);
      hiBaseWin = min(round(data.SWR.evPeak(i) + 0.5 * param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));

      data.C.SWR.event{i} = data.C.tSeries(loWin : hiWin);
      
      % Determine baseline:
      baseCell = vertcat(data.C.tSeries(loWin : loBaseWin-1), data.C.tSeries(hiBaseWin+1 : hiWin));
      
      % Only consider baseline data in lower 50% quantile:
      if param.pscEventPolarity
        baseCell(baseCell > quantile(baseCell, 0.50)) = [];
      else
        baseCell(baseCell < quantile(baseCell, 0.50)) = [];
      end
      data.C.SWR.baseline(i) = mean(baseCell);
      data.C.SWR.evNorm{i}   = data.C.SWR.event{i} - data.C.SWR.baseline(i);
      
      % Calculate area, defined as from -50:+50, and in 4 quarters: Q1=-100:-50, Q2=-50:0,Q3=0:50, Q4=50:100
      data.C.SWR.area(i)   = data.C.samplingInt * sum(sum(data.C.tSeries(loBaseWin : hiBaseWin) - data.C.SWR.baseline(i)));
      data.C.SWR.areaQ1(i) = data.C.samplingInt * sum(sum(data.C.tSeries(loWin : loBaseWin-1) - data.C.SWR.baseline(i)));
      data.C.SWR.areaQ2(i) = data.C.samplingInt * sum(sum(data.C.tSeries(loBaseWin : data.SWR.evPeak(i)-1) - data.C.SWR.baseline(i)));
      data.C.SWR.areaQ3(i) = data.C.samplingInt * sum(sum(data.C.tSeries(data.SWR.evPeak(i) : hiBaseWin) - data.C.SWR.baseline(i)));
      data.C.SWR.areaQ4(i) = data.C.samplingInt * sum(sum(data.C.tSeries(hiBaseWin+1 : hiWin) - data.C.SWR.baseline(i)));
      
      % Calculate running integrated charge over swrWindow (warning: computationally intensive)
      for j = 1 : length(data.C.SWR.event{i})
        loAreaWin = max(round(loWin + (j - 1) - 0.5 * param.swrWindow / data.LFP.samplingInt), 1);
        hiAreaWin = min(round(loWin + (j - 1) + 0.5 * param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
        data.C.SWR.evArea{i}(j) = data.C.samplingInt * sum(sum(data.C.tSeries(loAreaWin : hiAreaWin) - data.C.SWR.baseline(i)));
      end
      data.C.SWR.evArea{i} = data.C.SWR.evArea{i}';
    end

    % Calculate average evArea
    nSWRsValid = length(data.C.SWR.evArea);
    nSamples   = max(cellfun(@length, data.C.SWR.evArea));
    
    if length(data.C.SWR.evArea{1}) == nSamples
      data.C.SWR.evAreaAve = data.C.SWR.evArea{1};
      iStart = 2;
    else
      data.C.SWR.evAreaAve = data.C.SWR.evArea{2};
      nSWRsValid = nSWRsValid - 1;
      iStart = 3;
    end
    
    for i = iStart:length(data.C.SWR.evArea)
      if length(data.C.SWR.evArea{i}) == nSamples
        data.C.SWR.evAreaAve = data.C.SWR.evAreaAve + data.C.SWR.evArea{i}; 
      else
        nSWRsValid = nSWRsValid - 1;
      end
    end
    data.C.SWR.evAreaAve = data.C.SWR.evAreaAve / nSWRsValid;
    
    % Calculate SEM of evArea
    data.C.SWR.evAreaSEM = (data.C.SWR.evArea{iStart - 1} - data.C.SWR.evAreaAve).*(data.C.SWR.evArea{iStart - 1} - data.C.SWR.evAreaAve);
    for i = iStart:length(data.C.SWR.evArea)
      if length(data.C.SWR.evArea{i}) == nSamples
        data.C.SWR.evAreaSEM = data.C.SWR.evAreaSEM + (data.C.SWR.evArea{i} - data.C.SWR.evAreaAve).*(data.C.SWR.evArea{i} - data.C.SWR.evAreaAve);
      end
    end
    data.C.SWR.evAreaSEM = sqrt(data.C.SWR.evAreaSEM/(nSWRsValid * (nSWRsValid - 1)));

    % Transpose arrays
    data.C.SWR.event    = data.C.SWR.event';
    data.C.SWR.evNorm   = data.C.SWR.evNorm';
    data.C.SWR.evArea   = data.C.SWR.evArea';
    data.C.SWR.baseline = data.C.SWR.baseline';
    data.C.SWR.area     = data.C.SWR.area';
    data.C.SWR.areaQ1   = data.C.SWR.areaQ1';
    data.C.SWR.areaQ2   = data.C.SWR.areaQ2';
    data.C.SWR.areaQ3   = data.C.SWR.areaQ3';
    data.C.SWR.areaQ4   = data.C.SWR.areaQ4';
  end
  fprintf('done\n');
end

%% Import and process event file(s)
if param.importPSCOption %% && ~param.reAnalyzeOption %% Commented out want to have option to re-import events but may have unintended consequences
  
  % Import event files into matlab tables
  warning('off')
  fprintf(['importing PSC event file ' pscFileName '... ']);
  pscTable = readtable(pscFile,'Delimiter','\t','ReadRowNames', 1, 'TreatAsEmpty', 'Not found');
  fprintf('done\n');
  warning('on')
  
  % Re-initialize PSC data structures
  data.C.PSC.evStatusPeak = [];
  data.C.PSC.evStatusSum = [];
  data.C.PSC.evStart  = [];
  data.C.PSC.evPeak   = [];
  data.C.PSC.evEnd    = [];
  data.C.PSC.baseline = [];
  data.C.PSC.amp      = [];
  data.C.PSC.riseTau  = [];
  data.C.PSC.decayTau = [];
  data.C.PSC.area     = [];
  
  fprintf(['processing PSC events (file ' dataFileName ')... ']);
  for i = 1:size(pscTable,2)
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"EventStartTime"))
      pscStartTime = pscTable{:,i};
      for psc = 1:length(pscStartTime)
        if ~isempty(find(data.C.timing >= pscStartTime(psc), 1))
          data.C.PSC.evStart(psc) = find(data.C.timing >= pscStartTime(psc), 1);
        end
      end
      data.C.PSC.evStart = data.C.PSC.evStart';
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"EventEndTime"))
      pscEndTime = pscTable{:,i};
      for psc = 1:length(pscEndTime)
        if ~isempty(find(data.C.timing <= pscEndTime(psc), 1, 'last'))
          data.C.PSC.evEnd(psc) = find(data.C.timing <= pscEndTime(psc), 1, 'last');
        end
      end
      data.C.PSC.evEnd = data.C.PSC.evEnd';
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"Baseline"))
      data.C.PSC.baseline = pscTable{:,i};
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"PeakAmp"))
      data.C.PSC.amp = pscTable{:,i};
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"TimeOfPeak"))
      pscPeakTime = pscTable{:,i};
      for psc = 1:length(pscPeakTime)
        if ~isempty(find(data.C.timing >= pscPeakTime(psc), 1, 'last'))
          data.C.PSC.evPeak(psc) = find(data.C.timing >= pscPeakTime(psc), 1);
        end
      end
      data.C.PSC.evPeak = data.C.PSC.evPeak';
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"RiseTau"))
      data.C.PSC.riseTau = pscTable{:,i};
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"DecayTau"))
      data.C.PSC.decayTau = pscTable{:,i};
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"Area"))
      data.C.PSC.area = pscTable{:,i};
    end
    
  end
  fprintf('done\n');

  % Compute PSC status arrays
  data.C.PSC.evStatusPeak = zeros(length(data.C.timing),1);
  data.C.PSC.evStatusSum  = zeros(length(data.C.timing),1);
  
  for ev = 1:length(data.C.PSC.evStart)
    
    % Will range between 0 and 1 at time of PSC peak
    data.C.PSC.evStatusPeak(data.C.PSC.evPeak(ev)) = 1;
    
    % Will range from 0 to # concurrent PSCs
    for i = data.C.PSC.evStart(ev) : data.C.PSC.evEnd(ev)
      data.C.PSC.evStatusSum(i) = data.C.PSC.evStatusSum(i) + 1;
    end
  end
end

%% Calculate standard SWR status (if selected)
if param.useSWRWindowOption && param.swrPSCOption
  data.SWR.evStatusStand  = zeros(length(data.SWR.evStatus),1);
  data.SWR.evStartStand   = [];
  data.SWR.evEndStand     = [];
  
  if ~isempty(data.SWR.evStart)
    data.SWR.evStartStand = zeros(length(data.SWR.evStart),1);
    data.SWR.evEndStand   = zeros(length(data.SWR.evStart),1);
    
    for swr = 1:length(data.SWR.evStart)
      data.SWR.evStartStand(swr) = max(round(data.SWR.evPeak(swr) - 0.5 * param.swrWindow / data.LFP.samplingInt), 1);
      data.SWR.evEndStand(swr) = min(round(data.SWR.evPeak(swr) + 0.5 * param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
      data.SWR.evStatusStand(data.SWR.evStartStand(swr):data.SWR.evEndStand(swr)) = 1;
    end
  end
end

%% Parse events into SWR blocks (currently mandatory)
if param.parsePSCOption
  data.C.SWR.PSC.evPeak = [];
  data.C.SWR.PSC.evStatusPeak = [];
  data.C.SWR.PSC.evStatusSum  = [];
  
  if ~isempty(data.SWR.evStart)
    
    % Initialize event locked data window cell arrays
    data.C.SWR.PSC.evPeak{length(data.SWR.evStart)}       = [];
    data.C.SWR.PSC.evStatusPeak{length(data.SWR.evStart)} = [];
    data.C.SWR.PSC.evStatusSum{length(data.SWR.evStart)}  = [];
    
    % For each SWR event, parse PSC events
    for swr = 1:length(data.SWR.evStart)
      
      % Re-calculate start and end windows for each SWR event:
      loWin = max(round(data.SWR.evPeak(swr) - param.swrWindow / data.LFP.samplingInt), 1);
      hiWin = min(round(data.SWR.evPeak(swr) + param.swrWindow / data.LFP.samplingInt), length(data.C.PSC.evStatusPeak));
      
      % Parse PSC peaks into SWR cell array
      data.C.SWR.PSC.evPeak{swr}       = data.C.PSC.evPeak(data.C.PSC.evPeak>=loWin & data.C.PSC.evPeak<=hiWin) - (loWin - 1);
      data.C.SWR.PSC.evStatusPeak{swr} = data.C.PSC.evStatusPeak(loWin : hiWin);
      data.C.SWR.PSC.evStatusSum{swr}  = data.C.PSC.evStatusSum(loWin : hiWin);
      
    end
    
    data.C.SWR.PSC.evPeak       = data.C.SWR.PSC.evPeak';
    data.C.SWR.PSC.evStatusPeak = data.C.SWR.PSC.evStatusPeak';
    data.C.SWR.PSC.evStatusSum  = data.C.SWR.PSC.evStatusSum';
    
    % Calculate cumulative event PSC status over all SWR events:
    nSamples = max(cellfun(@length, data.C.SWR.PSC.evStatusPeak));
    Bins     = linspace(data.SWR.evTiming(1), data.SWR.evTiming(length(data.SWR.evTiming)), param.nBins)';
    binWidth = Bins(2) - Bins(1);
    
    data.C.SWR.PSC.evStatusPeakSum = zeros(nSamples, 1);
    data.C.SWR.PSC.evHist          = zeros(param.nBins, 1);
    data.C.SWR.PSC.nPSCsHist       = 0;
    
    for swr = 1:length(data.C.SWR.PSC.evStatusPeak)
      
      if length(data.C.SWR.PSC.evStatusPeak{swr}) == nSamples
        data.C.SWR.PSC.evStatusPeakSum = data.C.SWR.PSC.evStatusPeakSum + data.C.SWR.PSC.evStatusPeak{swr};
      end
      
      nPSCSWR = length(data.C.SWR.PSC.evPeak{swr});
      data.C.SWR.PSC.nPSCsHist = data.C.SWR.PSC.nPSCsHist + nPSCSWR;
      
      for psc = 1:nPSCSWR
        peakTime = data.SWR.evTiming(data.C.SWR.PSC.evPeak{swr}(psc));
        data.C.SWR.PSC.evHist = data.C.SWR.PSC.evHist + double(peakTime >= Bins & peakTime < Bins + binWidth);
      end
    end
    
    % Normalize histogram to all SWRs to get prob. of spiking per bin
    data.C.SWR.PSC.evHist = data.C.SWR.PSC.evHist / length(data.C.SWR.PSC.evPeak);
  end
end

%% Calculate overlap of SWRs and PSCs
if param.swrPSCOption
  
  % Calculate overlap for SWRs and PSCs
  fprintf(['detecting PSCs coincident with SWRs (file ' dataFileName ')... ']);
  if param.useSWRDurationOption
    [data.C.PSC.evStatusPeakC, data.C.PSC.evStartC, data.C.PSC.evEndC, data.C.PSC.evIndex] = eventOverlap(data.SWR.evStatus, data.SWR.evStart, data.SWR.evEnd, ...
      data.C.PSC.evStatusPeak, data.C.PSC.evPeak, data.C.PSC.evPeak, data.C.timing, 2);
  elseif param.useSWRWindowOption
    [data.C.PSC.evStatusPeakC, data.C.PSC.evStartC, data.C.PSC.evEndC, data.C.PSC.evIndex] = eventOverlap(data.SWR.evStatusStand, data.SWR.evStartStand, data.SWR.evEndStand, ...
      data.C.PSC.evStatusPeak, data.C.PSC.evPeak, data.C.PSC.evPeak, data.C.timing, 2);
  end
  fprintf('done\n');
  
  % Find SWRs that coincide with at least one PSC
  fprintf(['detecting SWRs coincident with PSCs (file ' dataFileName ')... ']);
  if param.useSWRDurationOption
    [data.SWR.PSC.evStatusC, data.SWR.PSC.evStartC, data.SWR.PSC.evEndC, data.SWR.PSC.evIndex] = eventOverlap(data.SWR.evStatus, data.SWR.evStart, data.SWR.evEnd, ...
      data.C.PSC.evStatusPeak, data.C.PSC.evPeak, data.C.PSC.evPeak, data.C.timing, 1);
  elseif param.useSWRWindowOption
    [data.SWR.PSC.evStatusC, data.SWR.PSC.evStartC, data.SWR.PSC.evEndC, data.SWR.PSC.evIndex] = eventOverlap(data.SWR.evStatusStand, data.SWR.evStartStand, data.SWR.evEndStand, ...
      data.C.PSC.evStatusPeak, data.C.PSC.evPeak, data.C.PSC.evPeak, data.C.timing, 1);
  end
  fprintf('done\n');
  
  data.SWR.PSC.nEventsA = size(data.SWR.evStart,1);
  data.SWR.PSC.nEventsC = size(data.SWR.PSC.evStartC,1);
  data.C.PSC.nEventsA   = size(data.C.PSC.evStart,1);
  data.C.PSC.nEventsC   = size(data.C.PSC.evStartC,1);
  
  % Calculate simplified event matrix (currently mandatory)
  if param.calcEvMatrixOption
    data.C.PSC.swrMatrix  = zeros(data.C.PSC.nEventsA, 1);
    data.SWR.PSC.evMatrix = zeros(data.SWR.PSC.nEventsA, 1);
    
    if param.useSWRDurationOption
      coincStatus = data.C.PSC.evStatusPeak .* data.SWR.evStatus;
    elseif param.useSWRWindowOption
      coincStatus = data.C.PSC.evStatusPeak .* data.SWR.evStatusStand;
    end
    
    for psc = 1:data.C.PSC.nEventsA
      if coincStatus(data.C.PSC.evPeak(psc)) > 0
        data.C.PSC.swrMatrix(psc) = 1;
      end
    end
    
    for swr = 1:data.SWR.PSC.nEventsA
      if (sum(coincStatus(data.SWR.evStart(swr) : data.SWR.evEnd(swr))) > 0)
        data.SWR.PSC.evMatrix(swr) = 1;
      end
    end
  end
  
  
  %% Calculate ripple phase of PSCs
  if ~isempty(data.SWR.evStart) && isfield(data.R.SWR.phase, 'evPhase')
    fprintf(['detecting oscillation phase of PSCs (file ' dataFileName ')... ']);
    
    % Initialize phase variables
    data.C.SWR.PSC.ripplePhase  = [];
    data.C.SWR.PSC.ripplePhase{length(data.C.SWR.PSC.evPeak)} = [];
    data.C.SWR.PSC.ripplePhase = data.C.SWR.PSC.ripplePhase';
    
    data.C.PSC.ripplePhase  = NaN*ones(length(data.C.PSC.evPeak), 1);
    data.C.PSC.ripplePhaseX = NaN*ones(length(data.C.PSC.evPeak), 1);
    data.C.PSC.ripplePhaseY = NaN*ones(length(data.C.PSC.evPeak), 1);
    
    % Calculate gamma phase if available (have to iterate through structure to not throw error)
    if isfield(data, 'gamma')
      if isfield(data.gamma, 'SWR')
        if isfield(data.gamma.SWR, 'phase')
          data.C.SWR.PSC.gammaPhase = [];
          data.C.SWR.PSC.gammaPhase{length(data.C.SWR.PSC.evPeak)} = [];
          data.C.SWR.PSC.gammaPhase = data.C.SWR.PSC.gammaPhase';
          
          data.C.PSC.gammaPhase  = NaN*ones(length(data.C.PSC.evPeak), 1);
          data.C.PSC.gammaPhaseX = NaN*ones(length(data.C.PSC.evPeak), 1);
          data.C.PSC.gammaPhaseY = NaN*ones(length(data.C.PSC.evPeak), 1);
        end
      end
    end
    
    nSamples = max(cellfun(@length, data.SWR.event));
    
    % Thresholds for phase using previously saved param.sdMult from analyzeLFP
    minAmpR = std(data.R.tSeries) * data.LFP.param.sdMult;
    if isfield(data.C.SWR.PSC, 'gammaPhase') minAmpG = std(data.gamma.tSeries) * data.LFP.param.sdMult; end
    
    for psc1 = 1:data.C.PSC.nEventsA
      if (data.C.PSC.swrMatrix(psc1) == 1) % Coincident PSC
        pscInd = data.C.PSC.evIndex(:, 2);
        swr    = data.C.PSC.evIndex(pscInd == psc1, 1);
        loWin  = max(round(data.SWR.evPeak(swr) - param.swrWindow / data.LFP.samplingInt), 1);
        psc2   = find(data.C.SWR.PSC.evPeak{swr} == data.C.PSC.evPeak(psc1) - (loWin - 1));
        
        % Only consider non-truncated SWR windows with PSCs
        if(length(data.SWR.event{swr}) == nSamples) && ~isempty(psc2)
          
          peakTime = data.SWR.evTiming(data.C.SWR.PSC.evPeak{swr}(psc2));
          
          % Only select PSCs for which trough-to-peak amplitude is greater than minAmpR
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
              data.C.SWR.PSC.ripplePhase{swr}(psc2) = data.R.SWR.phase.evPhase{swr}(data.C.SWR.PSC.evPeak{swr}(psc2));
              data.C.PSC.ripplePhase(psc1)  = data.C.SWR.PSC.ripplePhase{swr}(psc2);
              data.C.PSC.ripplePhaseX(psc1) = cos(data.C.PSC.ripplePhase(psc1));
              data.C.PSC.ripplePhaseY(psc1) = sin(data.C.PSC.ripplePhase(psc1));
            end
          end
          
          % Calculate gamma phase if available
          if isfield(data.C.SWR.PSC, 'gammaPhase')
            
            % Only select PSCs for which trough-to-peak amplitude is greater than minAmpG
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
                data.C.SWR.PSC.gammaPhase{swr}(psc2) = data.gamma.SWR.phase.evPhase{swr}(data.C.SWR.PSC.evPeak{swr}(psc2));
                data.C.PSC.gammaPhase(psc1) = data.C.SWR.PSC.gammaPhase{swr}(psc2);
                data.C.PSC.gammaPhaseX(psc1) = cos(data.C.PSC.gammaPhase(psc1));
                data.C.PSC.gammaPhaseY(psc1) = sin(data.C.PSC.gammaPhase(psc1));
              end
            end
          end
        end
      end
    end
    fprintf('done\n');
  end
end

%% Optional filters
if param.gammaOption
  % Apply Gaussian filter to extract gamma signal
  fprintf(['band-pass filtering gamma between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.gammaLim1, param.gammaLim2);
  if ~isfield(data,'gammaC') data.gammaC = struct; end
  data.gammaC.tSeries = gaussianFilt(data.C.tSeries, param.gammaLim1, param.gammaLim2, data.C.samplingInt, 1);
  data.gammaC.tPower  = bandpower(data.gammaC.tSeries);
  data.gammaC.lim1    = param.gammaLim1;
  data.gammaC.lim2    = param.gammaLim2;
  fprintf('done\n');
end

if param.rOption
  % Apply Gaussian filter to extract ripple signal
  fprintf(['band-pass filtering ripple between %4.1f-%4.1fHz (file ' dataFileName ')... '], param.rLim1, param.rLim2);
  if ~isfield(data,'RC') data.RC = struct; end
  data.RC.tSeries = gaussianFilt(data.C.tSeries, param.rLim1, param.rLim2, data.C.samplingInt, 1);
  data.RC.tPower  = bandpower(data.RC.tSeries);
  data.RC.lim1    = param.rLim1;
  data.RC.lim2    = param.rLim2;
  fprintf('done\n');
end

% Calculate SWR event-locked gamma, ripple, and/or spectrograms of PSC signal
if isfield(data,'SWR') && (isfield(data,'gammaC') || isfield(data,'RC') || param.spectOption)
  
  if isfield(data,'gammaC')
    if ~isfield(data.gammaC,'SWR') data.gammaC.SWR = struct; end
    data.gammaC.SWR.event = [];
    data.gammaC.SWR.power = [];
    if isfield(data,'gamma') data.gammaC.SWR.xCorr = []; end
  end
  
  if isfield(data,'RC')
    if ~isfield(data.RC,'SWR') data.RC.SWR = struct; end
    data.RC.SWR.event = [];
    data.RC.SWR.power = [];
    data.RC.SWR.xCorr = [];
  end
  
  if ~isnan(data.SWR.evStart)
    
    % Initialize event-locked data window and correlation cell arrays
    nSamples = max(cellfun(@length, data.SWR.event));
    if isfield(data,'gammaC')  
      data.gammaC.SWR.event{length(data.SWR.evStart)} = []; 
      if isfield(data,'gamma') 
        data.gammaC.SWR.xCorr{length(data.SWR.evStart)} = []; 
        data.gammaC.SWR.oCorr   = NaN * ones(length(data.SWR.event), 2);
        data.gammaC.SWR.minCorr = NaN * ones(length(data.SWR.event), 2);
        data.gammaC.SWR.maxCorr = NaN * ones(length(data.SWR.event), 2);
      end
    end
    
    if isfield(data,'RC') 
      data.RC.SWR.event{length(data.SWR.evStart)} = []; 
      data.RC.SWR.xCorr{length(data.SWR.evStart)} = [];
      data.RC.SWR.oCorr   = NaN * ones(length(data.SWR.event), 2);
      data.RC.SWR.minCorr = NaN * ones(length(data.SWR.event), 2);
      data.RC.SWR.maxCorr = NaN * ones(length(data.SWR.event), 2);
    end
    
    for swr = 1:length(data.SWR.evStart)
      loWin = max(round(data.SWR.evPeak(swr) - param.swrWindow / data.LFP.samplingInt), 1);
      hiWin = min(round(data.SWR.evPeak(swr) + param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
      loBaseWin = max(round(data.SWR.evPeak(swr) - 0.5 * param.swrWindow / data.LFP.samplingInt), 1);
      hiBaseWin = min(round(data.SWR.evPeak(swr) + 0.5 * param.swrWindow / data.LFP.samplingInt), length(data.LFP.tSeries));
      
      if isfield(data,'gammaC')
        data.gammaC.SWR.event{swr} = data.gammaC.tSeries(loWin : hiWin);
        data.gammaC.SWR.power(swr) = bandpower(data.gammaC.tSeries(loBaseWin : hiBaseWin));
        
        % Cross-correlation of gamma between LFP and cell:
        if isfield(data,'gamma') && (length(data.gamma.SWR.event{swr}) == nSamples)
          data.gammaC.SWR.xCorr{swr} = xcorr(data.gamma.SWR.event{swr}, data.gammaC.SWR.event{swr}, round(0.5*(length(data.gammaC.SWR.event{swr}) - 1)), 'normalized');
          
          data.gammaC.SWR.oCorr(swr,1) = data.gammaC.SWR.xCorr{swr}(0.5*(length(data.gammaC.SWR.event{swr}) - 1) + 1);
          data.gammaC.SWR.oCorr(swr,2) = 0;
          
          [data.gammaC.SWR.maxCorr(swr,1), corrInd] = max(data.gammaC.SWR.xCorr{swr});
          data.gammaC.SWR.maxCorr(swr,2) = data.SWR.evTiming(corrInd);
          
          [data.gammaC.SWR.minCorr(swr,1), corrInd] = min(data.gammaC.SWR.xCorr{swr});
          data.gammaC.SWR.minCorr(swr,2) = data.SWR.evTiming(corrInd);
        end
      end
      
      if isfield(data,'RC')
        data.RC.SWR.event{swr} = data.RC.tSeries(loWin : hiWin);
        data.RC.SWR.power(swr) = bandpower(data.RC.tSeries(loBaseWin : hiBaseWin));
        
        % Cross-correlation of ripple between LFP and cell:
        if length(data.R.SWR.event{swr}) == nSamples
          data.RC.SWR.xCorr{swr} = xcorr(data.R.SWR.event{swr}, data.RC.SWR.event{swr}, round(0.5*(length(data.RC.SWR.event{swr}) - 1)), 'normalized');
          
          data.RC.SWR.oCorr(swr,1) = data.RC.SWR.xCorr{swr}(0.5*(length(data.RC.SWR.event{swr}) - 1) + 1);
          data.RC.SWR.oCorr(swr,2) = 0;
          
          [data.RC.SWR.maxCorr(swr,1), corrInd] = max(data.RC.SWR.xCorr{swr});
          data.RC.SWR.maxCorr(swr,2) = data.SWR.evTiming(corrInd);
          
          [data.RC.SWR.minCorr(swr,1), corrInd] = min(data.RC.SWR.xCorr{swr});
          data.RC.SWR.minCorr(swr,2) = data.SWR.evTiming(corrInd);
        end
      end
    end
    
    % Transpose arrays:
    if isfield(data,'gammaC')
      data.gammaC.SWR.event = data.gammaC.SWR.event';
      data.gammaC.SWR.power = data.gammaC.SWR.power';
      if isfield(data,'gamma') data.gammaC.SWR.xCorr = data.gammaC.SWR.xCorr'; end
    end
    
    if isfield(data,'RC')
      data.RC.SWR.event = data.RC.SWR.event';
      data.RC.SWR.power = data.RC.SWR.power';
      if isfield(data,'RC') data.RC.SWR.xCorr = data.RC.SWR.xCorr'; end
    end
    
    % Spectral analysis:
    if param.spectOption
      fprintf(['spectral analysis of SWR-locked events (file ' dataFileName ')... ']);
      fRange = param.spectLim1 : param.spectLim2;
      [data.C.SWR, ~] = calcSpect(data.C.SWR, [], fRange, data.LFP.param.Fs, 3, 0);
      data.C.SWR = calcEvFFT(data.C.SWR, data.LFP.param, param.spectLim1, param.spectLim2);
      
      if isfield(data,'gammaC')
        data.gammaC.SWR = calcEvFFT(data.gammaC.SWR, data.LFP.param, data.gammaC.lim1, data.gammaC.lim2);
        data.gammaC.SWR = calcEvPhase(data.gammaC.SWR, data.SWR, data.gammaC.lim1, data.gammaC.lim2);
      end
      
      if isfield(data,'RC')
        data.RC.SWR = calcEvFFT(data.RC.SWR, data.LFP.param, data.RC.lim1, data.RC.lim2);
        data.RC.SWR = calcEvPhase(data.RC.SWR, data.SWR, data.RC.lim1, data.RC.lim2);
      end
      
      fprintf('done\n');
    end
    
    if isfield(data,'gammaC')
      data.gammaC.SWR = orderStruct(data.gammaC.SWR);
      data.gammaC = orderStruct(data.gammaC);
    end
    
    if isfield(data,'RC')
      data.RC.SWR = orderStruct(data.RC.SWR);
      data.RC = orderStruct(data.RC);
    end
    
  end
end

% Re-order structure arrays:
data         = orderStruct(data);
data.C       = orderStruct(data.C);
data.C.PSC   = orderStruct(data.C.PSC);
data.C.param = orderStruct(data.C.param);
data.SWR     = orderStruct(data.SWR);
data.SWR.PSC = orderStruct(data.SWR.PSC);

if param.parsePSCOption
  data.C.SWR.PSC = orderStruct(data.C.SWR.PSC);
end

%% Save matlab file
if all(saveFile)
  fprintf(['saving file ' dataFileName '... ']);
  save(saveFile,'-struct','data');
  fprintf('done\n');
end

%% Export PSC event table
if all(expPSCFile) && param.expPSCEvOption
  fprintf(['exporting PSC event file ' expPSCFileName '... ']);
  exportPSCEvents(data, saveFile, expPSCFile);
  fprintf('done\n');
end

%% Export SWR event file
if (all(expSWRFile) && param.expSWREvOption)
  fprintf(['exporting updated SWR events (file ' expSWRFileName ')... ']);
  exportSWREvents(data, saveFile, expSWRFile)
  fprintf('done\n');
end

end
