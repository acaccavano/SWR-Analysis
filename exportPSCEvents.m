function exportPSCEvents(data, saveFile, exportFile)

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_pscEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of PSC events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('No PSC events to be exported - no file selected'); end
end

varNames = {'evStart', 'evPeak', 'evEnd', 'swrCoinc', 'Baseline_pA', 'PeakAmp_pA', 'RiseTau_ms', 'DecayTau_ms', 'area'};
outTable = table(data.C.PSC.evStart, data.C.PSC.evPeak, data.C.PSC.evEnd, data.C.PSC.swrMatrix, data.C.PSC.baseline, data.C.PSC.amp, data.C.PSC.riseTau, data.C.PSC.decayTau, data.C.PSC.area, 'VariableNames', varNames);

if isfield(data.C.PSC, 'gammaPhase') 
  outTable = [outTable table(data.C.PSC.gammaPhase, data.C.PSC.gammaPhaseX, data.C.PSC.gammaPhaseY, 'VariableNames', {'gammaPhase', 'gammaPhaseX', 'gammaPhaseY'})];
end

if isfield(data.C.PSC, 'ripplePhase') 
  outTable = [outTable table(data.C.PSC.ripplePhase, data.C.PSC.ripplePhaseX, data.C.PSC.ripplePhaseY, 'VariableNames', {'ripplePhase', 'ripplePhaseX', 'ripplePhaseY'})];
end

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end

