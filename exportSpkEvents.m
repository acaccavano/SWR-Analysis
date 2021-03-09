function exportSpkEvents(data, saveFile, exportFile)
%% exportSpkEvents(data, saveFile, exportFile)
%
%  Function to export csv file of all spike event stats that are available

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

varNames = {'evStart', 'evPeak', 'evEnd'};
outTable = table(data.C.spike.evStartA, data.C.spike.evPeakA, data.C.spike.evEndA, 'VariableNames', varNames);

if isfield(data.C.spike, 'swrMatrix')
  outTable = [outTable table(data.C.spike.swrMatrix, 'VariableNames', {'swrCoinc'})];
end

if isfield(data.C.spike, 'theta')
  outTable = [outTable table(data.C.spike.theta.phase, data.C.spike.theta.phaseX, data.C.spike.theta.phaseY, 'VariableNames', {'theta.phase', 'theta.phaseX', 'theta.phaseY'})];
end

if isfield(data.C.spike, 'beta')
  outTable = [outTable table(data.C.spike.beta.phase, data.C.spike.beta.phaseX, data.C.spike.beta.phaseY, 'VariableNames', {'beta.phase', 'beta.phaseX', 'beta.phaseY'})];
end

if isfield(data.C.spike, 'gamma')
  outTable = [outTable table(data.C.spike.gamma.phase, data.C.spike.gamma.phaseX, data.C.spike.gamma.phaseY, 'VariableNames', {'gamma.phase', 'gamma.phaseX', 'gamma.phaseY'})];
end

if isfield(data.C.spike, 'hgamma')
  outTable = [outTable table(data.C.spike.hgamma.phase, data.C.spike.hgamma.phaseX, data.C.spike.hgamma.phaseY, 'VariableNames', {'hgamma.phase', 'hgamma.phaseX', 'hgamma.phaseY'})];
end

if isfield(data.C.spike, 'R')
  outTable = [outTable table(data.C.spike.R.phase, data.C.spike.R.phaseX, data.C.spike.R.phaseY, 'VariableNames', {'R.phase', 'R.phaseX', 'R.phaseY'})];
end

if isfield(data.C.spike, 'fR')
  outTable = [outTable table(data.C.spike.fR.phase, data.C.spike.fR.phaseX, data.C.spike.fR.phaseY, 'VariableNames', {'fR.phase', 'fR.phaseX', 'fR.phaseY'})];
end

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

