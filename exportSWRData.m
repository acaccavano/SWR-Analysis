function exportSWRData(data, param, saveFile, expDataFile)

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSWRData');
end

if (nargin < 3)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_swrData.txt'];
  [exportName, exportPath] = uiputfile('.txt','Select txt file to export episodic SWR data', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('Please select valid file'); end
end

for i = 1:length(data.SWR.event)
  dataOutSize(i) = length(data.SWR.event{i});
end

% Check for partial first and last events, and remove if so:
if (dataOutSize(1) < max(dataOutSize))
  data.SWR.event(1) = [];
  data.R.SWR.event(1) = [];
  dataOutSize(1) = [];
  if param.cellOption
    data.C.SWR.event(1) = [];
    data.C.SWR.evNorm(1) = [];
  end
end
if (dataOutSize(length(data.SWR.event)) < max(dataOutSize))
  data.SWR.event(length(data.SWR.event)) = [];
  data.R.SWR.event(length(data.R.SWR.event)) = [];
  if param.cellOption
    data.C.SWR.event(length(data.C.SWR.event)) = [];
    data.C.SWR.evNorm(length(data.C.SWR.evNorm)) = [];
  end
end

dataOut = (0: data.LFP.samplingInt : 2*param.swrWindow)';

% Output table names
if param.cellOption
  tableVarNames = cell(1, 1 + 4*length(data.SWR.event));
else
  tableVarNames = cell(1, 1 + 2*length(data.SWR.event));
end

tableVarNames{1} = 'time';
nameInd = 2;

% Sort cell array into output table
for i = 1:length(data.SWR.event)
  % LFP event-locked data:
  dataOut = horzcat(dataOut, data.SWR.event{i});
  tableVarNames{nameInd} = ['LFP_' num2str(i)];
  nameInd = nameInd + 1;
  % ripple event-locked data:
  dataOut = horzcat(dataOut, data.R.SWR.event{i});
  tableVarNames{nameInd} = ['R_' num2str(i)];
  nameInd = nameInd + 1;
  
  if param.cellOption
    % cell event-locked data:
    dataOut = horzcat(dataOut, data.C.SWR.event{i});
    tableVarNames{nameInd} = ['Cell_' num2str(i)];
    nameInd = nameInd + 1;
    % cell (normalized) event-locked data:
    dataOut = horzcat(dataOut, data.C.SWR.evNorm{i});
    tableVarNames{nameInd} = ['CellNorm_' num2str(i)];
    nameInd = nameInd + 1;
  end
end

% Export data table
dataOutTable = array2table(dataOut, 'VariableNames', tableVarNames);
writetable(dataOutTable, expDataFile, 'Delimiter', '\t');

end

