function exportStimEvents(data, saveFile, exportFile)

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportStimEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_StimEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of stimulation events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('No stimulation events to be exported - no file selected'); end
end

varNames  = {'stimStart', 'stimEnd', 'stimPeak', 'lfpAmp', 'lfpBL', 'lfpSlope'};
cellInd   = 1:data.Ca.nChannels;

peakNames{data.Ca.nChannels} = [];
areaNames{data.Ca.nChannels} = [];
fireNames{data.Ca.nChannels} = [];

for i = 1:data.Ca.nChannels
  peakNames{i} = strcat('dFoF_Peak', num2str(cellInd(i)));
  areaNames{i} = strcat('dFoF_Area', num2str(cellInd(i)));
  fireNames{i} = strcat('Fire', num2str(cellInd(i)));
end
  
varNames = [varNames peakNames areaNames fireNames];
varArray = horzcat(data.stim.evStart, data.stim.evEnd, data.stim.evPeak, data.stim.lfpAmp, data.stim.lfpBL, data.stim.lfpSlope, data.Ca.stim.evPeakStim, data.Ca.stim.evAreaStim, data.stim.Ca.evMatrix);
outTable = array2table(varArray, 'VariableNames', varNames);

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

