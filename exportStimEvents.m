function exportStimEvents(data, saveFile, exportFile)
%% exportStimEvents(data, saveFile, exportFile)
%
%  Function to export csv file of all stim event stats

% Handle input arguments - if not entered
if (nargin < 3); exportFile = []; end
if (nargin < 2); saveFile   = []; end
if (nargin < 1); data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportStimEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_StimEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of stimulation events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile); error('No stimulation events to be exported - no file selected'); end
end

colNames{data.stim.Ca.nEventsA} = [];

for i = 1:data.stim.Ca.nEventsA
  colNames{i} = strcat('stim', num2str(i));
end

rowNames = {'stimStart'; 'stimEnd'; 'stimPeak'; 'lfpAmp'; 'lfpBL'; 'lfpSlope'; 'nCellsActive'};

peakNames{data.Ca.nChannels, 1} = [];
areaNames{data.Ca.nChannels, 1} = [];
boolNames{data.Ca.nChannels, 1} = [];

for i = 1:data.Ca.nChannels
  peakNames{i} = strcat('Cell', num2str(i), '_Peak');
  areaNames{i} = strcat('Cell', num2str(i), '_Area');
  boolNames{i} = strcat('Cell', num2str(i));
end
  
rowNames = [rowNames; peakNames; areaNames; boolNames];
varArray = vertcat(data.stim.evStart', data.stim.evEnd', data.stim.evPeak', data.stim.lfpAmp', data.stim.lfpBL', data.stim.lfpSlope', data.stim.Ca.nCellsC', data.Ca.stim.evPeakStim', data.Ca.stim.evAreaStim', data.stim.Ca.evMatrix');
outTable = array2table(varArray, 'VariableNames', colNames, 'RowNames', rowNames);

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp, 'VariableNames', outTable.Properties.VariableNames, 'RowNames', outTable.Properties.RowNames);

writetable(outTable, exportFile, 'Delimiter', ',', 'WriteVariableNames', 1, 'WriteRowNames', 1);

end

