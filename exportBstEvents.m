function exportBstEvents(data, saveFile, exportFile)

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_bstEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of burst events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('No burst events to be exported - no file selected'); end
end

varNames = {'evStart', 'evEnd', 'swrCoinc', 'nSpikesInBurst', 'intraBurstInterval_ms'};
outTable = table(data.C.burst.evStartA, data.C.burst.evEndA, data.C.burst.swrMatrix, data.C.burst.nSpike, data.C.burst.intraBI, 'VariableNames', varNames);

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

