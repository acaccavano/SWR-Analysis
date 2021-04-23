function exportCaEvents(data, saveFile, exportFile)
%% exportCaEvents(data, saveFile, exportFile)
%
%  Function to export csv file of all calcium event stats, averaged per cell, that are available

% Handle input arguments - if not entered
if (nargin < 3); exportFile = []; end
if (nargin < 2); saveFile   = []; end
if (nargin < 1); data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_CaEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export table of Calcium events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile); error('No Calcium events to be exported - no file selected'); end
end

cellInd  = 1:length(data.Ca.nEvents);
outTable = table(cellInd', 'VariableNames', {'Cell'});

if isfield(data.Ca, 'cellType'); outTable = [outTable table(data.Ca.cellType', 'VariableNames', {'CellType'})]; end

varNames = {'nEvents', 'frequency_Hz', 'aveIEI_s', 'aveAmplitude_dFoF', 'aveDuration_s'};
outTable = [outTable table(data.Ca.nEvents', data.Ca.frequency', data.Ca.IEIAve', data.Ca.ampAve', data.Ca.durAve'/1000, 'VariableNames', varNames)];

% If SWR-Ca correlation has been performed:
if isfield(data.Ca, 'SWR')
  varNames = {'nEvents_Align', 'nEvents_Coinc', 'fracEvents_Coinc'};
  outTable = [outTable table(data.Ca.SWR.nEventsA', data.Ca.SWR.nEventsC', data.Ca.SWR.fracEventsC', 'VariableNames', varNames)];
  
  if isfield(data.Ca.SWR, 'swr')
    varNames = {'nEventsSWR', 'freqSWR_Hz', 'aveAmpSWR_dFoF', 'aveDurSWR_s'};
    outTable = [outTable table(data.Ca.SWR.swr.nEvents', data.Ca.SWR.swr.frequency', data.Ca.SWR.swr.ampAve', data.Ca.SWR.swr.durAve'/1000, 'VariableNames', varNames)]; 
  end
  
  if isfield(data.Ca.SWR, 'spont') 
    varNames = {'nEventsSpont', 'freqSpont_Hz', 'aveAmpSpont_dFoF', 'aveDurSpont_s'};
    outTable = [outTable table(data.Ca.SWR.spont.nEvents', data.Ca.SWR.spont.frequency', data.Ca.SWR.spont.ampAve', data.Ca.SWR.spont.durAve'/1000, 'VariableNames', varNames)]; 
  end
end

% If Stim-Ca correlation has been performed:
if isfield(data.Ca, 'stim')
  varNames = {'nEvents_Align', 'nEvents_Coinc', 'fracEvents_Coinc'};
  outTable = [outTable table(data.Ca.stim.nEventsA', data.Ca.stim.nEventsC', data.Ca.stim.fracEventsC', 'VariableNames', varNames)];
  
  if isfield(data.Ca.stim, 'stim')  
    varNames = {'nEventsStim', 'freqStim_Hz', 'aveAmpStim_dFoF', 'aveDurStim_s'};
    outTable = [outTable table(data.Ca.stim.stim.nEvents', data.Ca.stim.stim.frequency', data.Ca.stim.stim.ampAve', data.Ca.stim.stim.durAve'/1000, 'VariableNames', varNames)]; 
  end
  
  if isfield(data.Ca.stim, 'spont') 
    varNames = {'nEventsSpont', 'freqSpont_Hz', 'aveAmpSpont_dFoF', 'aveDurSpont_s'};
    outTable = [outTable table(data.Ca.stim.spont.nEvents', data.Ca.stim.spont.frequency', data.Ca.stim.spont.ampAve', data.Ca.stim.spont.durAve'/1000, 'VariableNames', varNames)]; 
  end
end

% Replace NaN values with blanks
if isfield(data.Ca, 'cellType')
  tmp = table2cell(outTable(:, 3:end));
  tmp(isnan(cell2mat(tmp))) = {[]};
  outTable = [outTable(:, 1:2) array2table(tmp,'VariableNames',outTable.Properties.VariableNames(3:end))];
else
  tmp = table2cell(outTable);
  tmp(isnan(outTable.Variables)) = {[]};
  outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);
end

writetable(outTable, exportFile, 'Delimiter', ',');

end