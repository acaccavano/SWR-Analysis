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

if isfield(data.C.spike, 'gammaPhase')
  outTable = [outTable table(data.C.spike.gammaPhase, data.C.spike.gammaPhaseX, data.C.spike.gammaPhaseY, 'VariableNames', {'gammaPhase', 'gammaPhaseX', 'gammaPhaseY'})];
end

if isfield(data.C.spike, 'ripplePhase')
  outTable = [outTable table(data.C.spike.ripplePhase, data.C.spike.ripplePhaseX, data.C.spike.ripplePhaseY, 'VariableNames', {'ripplePhase', 'ripplePhaseX', 'ripplePhaseY'})];
end

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

