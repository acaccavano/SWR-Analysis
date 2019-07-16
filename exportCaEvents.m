function exportCaEvents(data, saveFile, exportFile)

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_CaEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of Calcium events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('No Calcium events to be exported - no file selected'); end
end

varNames = {'Cell', 'nEvents', 'frequency_Hz', 'aveIEI_s', 'aveAmplitude_dFoF', 'aveDuration_s', };
cellInd  = 1:length(data.Ca.nEvents);
outTable = table(cellInd', data.Ca.nEvents', data.Ca.frequency', data.Ca.IEIAve', data.Ca.ampAve', data.Ca.durAve'/1000, 'VariableNames', varNames);

% IF SWR-Ca correlation has been performed:
if isfield(data.Ca, 'SWR')
  outTable = [outTable table(data.Ca.SWR.nEventsA', data.Ca.SWR.nEventsC', data.Ca.SWR.fracEventsC', 'VariableNames', {'nEvents_Align', 'nEvents_Coinc', 'fracEvents_Coinc'})];
end

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

