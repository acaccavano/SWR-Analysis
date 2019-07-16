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
if ~isfield(param,'fileNum')              param.fileNum = 1; end
if ~isfield(param,'importPSCOption')      param.importPSCOption = 1; end
if ~isfield(param,'swrPSCOption')         param.swrPSCOption = 1; end
if ~isfield(param,'useSWRDurationOption') param.useSWRDurationOption = 0; end
if ~isfield(param,'useSWRWindowOption')   param.useSWRWindowOption = 1; end
if ~isfield(param,'swrWindow')            param.swrWindow = 100; end
if ~isfield(param,'parsePSCOption')       param.parsePSCOption = 1; end
if ~isfield(param,'calcEvMatrixOption')   param.calcEvMatrixOption = 1; end
if ~isfield(param,'expPSCEvOption')       param.expPSCEvOption = 1; end
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
if ~isfield(data.C,'PSC') data.C.spike = struct; end

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
if param.importPSCOption && ~param.reAnalyzeOption
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
    [exportName, exportPath] = uiputfile('.txt','Select *.csv file to export table of PSC events', defaultPath);
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

%% Import and process event file(s)
if param.importPSCOption && ~param.reAnalyzeOption
  
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
  
  fprintf(['processing PSC events (file ' dataFileName ')... ']);
  for i = 1:size(pscTable,2)
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"EventStartTime"))
      pscStartTime = pscTable{:,i};
      for ev = 1:length(pscStartTime)
        if ~isempty(find(data.C.timing >= pscStartTime(ev), 1))
          data.C.PSC.evStart(ev) = find(data.C.timing >= pscStartTime(ev), 1);
        end
      end
      data.C.PSC.evStart = data.C.PSC.evStart';
    end
    
    if ~isempty(strfind(pscTable.Properties.VariableNames{i},"EventEndTime"))
      pscEndTime = pscTable{:,i};
      for ev = 1:length(pscEndTime)
        if ~isempty(find(data.C.timing <= pscEndTime(ev), 1, 'last'))
          data.C.PSC.evEnd(ev) = find(data.C.timing <= pscEndTime(ev), 1, 'last');
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
      for ev = 1:length(pscPeakTime)
        if ~isempty(find(data.C.timing >= pscPeakTime(ev), 1, 'last'))
          data.C.PSC.evPeak(ev) = find(data.C.timing >= pscPeakTime(ev), 1);
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
    for psc1 = 1:data.C.PSC.nEventsA
      if (data.C.PSC.swrMatrix(psc1) == 1) % Coincident PSC
        pscInd = data.C.PSC.evIndex(:, 2);
        swr    = data.C.PSC.evIndex(pscInd == psc1, 1);
        loWin  = max(round(data.SWR.evPeak(swr) - param.swrWindow / data.LFP.samplingInt), 1);
        psc2   = find(data.C.SWR.PSC.evPeak{swr} == data.C.PSC.evPeak(psc1) - (loWin - 1));
        
        % Only consider non-truncated SWR windows with PSCs
        if(length(data.SWR.event{swr}) == nSamples) && ~isempty(psc2)
          data.C.SWR.PSC.ripplePhase{swr}(psc2) = data.R.SWR.phase.evPhase{swr}(data.C.SWR.PSC.evPeak{swr}(psc2));
          data.C.PSC.ripplePhase(psc1)  = data.C.SWR.PSC.ripplePhase{swr}(psc2);
          data.C.PSC.ripplePhaseX(psc1) = cos(data.C.PSC.ripplePhase(psc1));
          data.C.PSC.ripplePhaseY(psc1) = sin(data.C.PSC.ripplePhase(psc1));
          
          % Calculate gamma phase if available
          if isfield(data.C.SWR.PSC, 'gammaPhase')
            data.C.SWR.PSC.gammaPhase{swr}(psc2) = data.gamma.SWR.phase.evPhase{swr}(data.C.SWR.PSC.evPeak{swr}(psc2));
            data.C.PSC.gammaPhase(psc1) = data.C.SWR.PSC.gammaPhase{swr}(psc2);
            data.C.PSC.gammaPhaseX(psc1) = cos(data.C.PSC.gammaPhase(psc1));
            data.C.PSC.gammaPhaseY(psc1) = sin(data.C.PSC.gammaPhase(psc1));
          end
          
        end
      end
    end
    fprintf('done\n');
  end
  
end

% Re-order structure arrays:
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
