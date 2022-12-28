function processEvFile(param, evFile, expFile, statFile, duration)
%% processevFile(param, evFile, expFile, statFile, duration)
%
%  Function to clean and process csv file of events exported from Clampfit.
%  Removes duplicate events (with same peak time), and calculates stats for
%  each trace. Not called from GUI, must be ran as standalone on csv file
%  or from processEvBatch
%
%  Inputs: (all optional - will be prompted)
%   param    = structure containing all parameters - used for manual column selection (CAUTION - may change for different pClamp versions)
%   evFile   = full path to pClamp event file to import
%   expFile  = full path to pClamp event file to export (required, can be same to overwrite)
%   statFile = full path to file to save averaged trace stats (optional)
%   duration = duration of trace in s

%% Handle input arguments - if not entered
if (nargin < 5); duration = []; end
if (nargin < 4); statFile = []; end
if (nargin < 3); expFile  = []; end
if (nargin < 2); evFile   = []; end
if (nargin < 1); param    = struct; end

% Handle case in which empty variable is supplied:
if isempty(param); param  = struct; end

% Set default parameters if not specified:
if ~isfield(param,'traceCol');   param.traceCol   =  1;  end
if ~isfield(param,'categCol');   param.categCol   =  3;  end
if ~isfield(param,'pkAmpCol');   param.pkAmpCol   =  8;  end
if ~isfield(param,'pkTimeCol');  param.pkTimeCol  = 10;  end
if ~isfield(param,'halfWdCol');  param.halfWdCol  = 14;  end
if ~isfield(param,'rTauCol');    param.rTauCol    = 18;  end
if ~isfield(param,'dTauCol');    param.dTauCol    = 19;  end
if ~isfield(param,'rTimeCol');   param.rTimeCol   = 25;  end
if ~isfield(param,'dTimeCol');   param.dTimeCol   = 27;  end
if ~isfield(param,'areaCol');    param.areaCol    = 28;  end
if ~isfield(param,'IEICol');     param.IEICol     = 30;  end

if isempty(evFile)
  [evName, evPath] = uigetfile('.csv', 'Select *.csv file of exported events from Clampfit');
  evFile = [evPath evName];
  if ~all(evFile); error('No event file selected'); end
end

% Parse evFile to determine default save name
[parentPath, evFileName, ~] = parsePath(evFile);

% Prompt for change in expFile
if isempty(expFile)
  [expName, expPath] = uiputfile('.csv','Select same or alternate *.csv file to export processed events', evFile);
  expFile = [expPath expName];
  if ~all(expFile)
    warning('No processed event file to be exported - no file selected');
  else
    [parentPath, ~, ~] = parsePath(expFile);
  end
end

% Prompt for statFile
if isempty(statFile)
  defaultPath = [parentPath evFileName '_Stats.csv'];
  [statName, statPath] = uiputfile('.csv','Select *.csv file to save averaged stats', defaultPath);
  statFile = [statPath statName];
  if ~all(statFile); warning('No stat file to be saved - no file selected'); end
end

% Define duration
if isempty(duration)
  answer = inputdlg({'Enter trace duration in s'}, 'Input', [1 35], {'60.000'});
  duration = str2double(answer{1});
end

% Force all imported columns to be numeric:
tabOptions = detectImportOptions(evFile);
for i = 1:length(tabOptions.VariableTypes)
  tabOptions.VariableTypes{i} = 'double';
end
evTable = readtable(evFile, tabOptions);

% Remove duplicates with same time of peak:
evTable = sortrows(evTable,[param.traceCol param.pkTimeCol param.categCol]); % Sort based on trace, then time of peak, then category, ensures lower category events will have preference.
[~, ia, ~] = unique(horzcat(evTable{:, 1}, evTable{:, param.pkTimeCol}), 'rows');  % Unique value index of concatenation of trace and time of peak
evTable = evTable(ia, :); % Remove duplicates

% Calculate average statistics
nTrace    = max(evTable{:,1});
trace     = 1:nTrace;
trace     = trace';
bins      = 1 : nTrace + 1;
bins      = bins';
nEvent    = histcounts(evTable{:,1}, bins)';
freq      = nEvent / duration;
amp       = zeros(nTrace, 1);
halfWidth = zeros(nTrace, 1);
riseTau   = zeros(nTrace, 1);
decayTau  = zeros(nTrace, 1);
riseTime  = zeros(nTrace, 1);
decayTime = zeros(nTrace, 1);
area      = zeros(nTrace, 1);  
IEI       = zeros(nTrace, 1);  

for i = 1:nTrace
  amp(i)       = mean(evTable{evTable{:,1} == i, param.pkAmpCol},  'omitnan');
  halfWidth(i) = mean(evTable{evTable{:,1} == i, param.halfWdCol}, 'omitnan');
  riseTau(i)   = mean(evTable{evTable{:,1} == i, param.rTauCol},   'omitnan');
  decayTau(i)  = mean(evTable{evTable{:,1} == i, param.dTauCol},   'omitnan');
  riseTime(i)  = mean(evTable{evTable{:,1} == i, param.rTimeCol},  'omitnan');
  decayTime(i) = mean(evTable{evTable{:,1} == i, param.dTimeCol},  'omitnan');
  area(i)      = mean(evTable{evTable{:,1} == i, param.areaCol},   'omitnan');
  IEI(i)       = mean(evTable{evTable{:,1} == i, param.IEICol},    'omitnan');
end

% Export processed table, replacing NaN values with blanks:
if all(expFile)
  varNames = evTable.Properties.VariableNames;
  tmpCell  = table2cell(evTable);
  tmpCell(cell2mat(cellfun(@(x)any(isnan(x)),tmpCell,'UniformOutput',false))) = {[]};
  evTable = cell2table(tmpCell);
  evTable.Properties.VariableNames = varNames;
  writetable(evTable, expFile, 'Delimiter', ',');
end

% Save average statistics, replacing NaN values with blanks:
if all(statFile)
  varNames = {'Trace', 'nEvent', 'Freq_Hz', 'Amp_pA', 'HalfWidth_ms', 'RiseTau_ms', 'DecayTau_ms', 'RiseTime_ms', 'DecayTime_ms', 'Area_pAs', 'IEI_s'};
  statTable = table(trace, nEvent, freq, amp, halfWidth, riseTau, decayTau, riseTime, decayTime, area/1000, IEI/1000, 'VariableNames', varNames);
  tmpCell  = table2cell(statTable);
  tmpCell(cell2mat(cellfun(@(x)any(isnan(x)),tmpCell,'UniformOutput',false))) = {[]};
  statTable = cell2table(tmpCell);
  statTable.Properties.VariableNames = varNames;
  writetable(statTable, statFile, 'Delimiter', ',');
end

end
