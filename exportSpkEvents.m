function exportSpkEvents(data, saveFile, exportFile)

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_spkEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of spike events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('No spike events to be exported - no file selected'); end
end

varNames = {'evStart', 'evPeak', 'evEnd', 'swrCoinc'};
outTable = table(data.C.spike.evStartA, data.C.spike.evPeakA, data.C.spike.evEndA, data.C.spike.swrMatrix, 'VariableNames', varNames);

if isfield(data.C.spike, 'gamma')
  outTable = [outTable table(data.C.spike.gamma.phase, data.C.spike.gamma.phaseX, data.C.spike.gamma.phaseY, 'VariableNames', {'gamma.phase', 'gamma.phaseX', 'gamma.phaseY'})];
end

if isfield(data.C.spike, 'R')
  outTable = [outTable table(data.C.spike.R.phase, data.C.spike.R.phaseX, data.C.spike.R.phaseY, 'VariableNames', {'R.phase', 'R.phaseX', 'R.phaseY'})];
end

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

