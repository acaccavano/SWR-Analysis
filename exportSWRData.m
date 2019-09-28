function exportSWRData(data, param, expDataFile)

% Handle inputs
if (nargin < 3) expDataFile = []; end
if (nargin < 2) param = struct; end
if (nargin < 1) error('Supply data structure to use exportSWRData'); end

% Select export file if not supplied
if isempty(expDataFile)
  [parentPath, saveFileName, ~] = parsePath(data.saveFile);
  defaultName = [parentPath saveFileName '_swrData.txt'];
  [exportName, exportPath] = uiputfile('.txt','Select txt file to export episodic SWR data', defaultName);
  expDataFile = [exportPath exportName];
  if ~all(expDataFile) error('Please select valid file'); end
end

% Set default parameters
if ~isfield(param,'gammaOption')      param.gammaOption      = 1; end
if ~isfield(param,'rOption')          param.rOption          = 1; end
if ~isfield(param,'cellOption')       param.cellOption       = 1; end
if ~isfield(param,'cellRawOption')    param.cellRawOption    = 0; end
if ~isfield(param,'cellGammaOption')  param.cellGammaOption  = 1; end
if ~isfield(param,'cellRippleOption') param.cellRippleOption = 1; end
if ~isfield(param,'truncateEvs')      param.truncateEvs      = 1; end
if ~isfield(param,'maxNumEvs')        param.maxNumEvs        =  50; end
if ~isfield(param,'swrWindow')        param.swrWindow        = 100; end

if param.truncateEvs
  nEvs = param.maxNumEvs;
else
  nEvs = length(data.SWR.event);
end

for i = 1:nEvs
  dataOutSize(i) = length(data.SWR.event{i});
end

% Check for partial first and last events, and remove if so:
if (dataOutSize(1) < max(dataOutSize))
  data.SWR.event(1) = [];
  if param.gammaOption   data.gamma.SWR.event(1) = []; end
  if param.rOption       data.R.SWR.event(1)     = []; end
  if param.cellOption
    if param.cellRawOption    data.C.SWR.event(1)      = []; end
    data.C.SWR.evNorm(1) = [];
    if param.cellGammaOption  data.gammaC.SWR.event(1) = []; end
    if param.cellRippleOption data.RC.SWR.event(1)     = []; end
  end
  dataOutSize(1) = [];
  nEvs = nEvs - 1;
end

if (dataOutSize(nEvs) < max(dataOutSize))
  data.SWR.event(nEvs) = [];
  if param.gammaOption   data.gamma.SWR.event(nEvs) = []; end
  if param.rOption       data.R.SWR.event(nEvs)     = []; end
  if param.cellOption
    if param.cellRawOption    data.C.SWR.event(nEvs)      = []; end
    data.C.SWR.evNorm(nEvs) = [];
    if param.cellGammaOption  data.gammaC.SWR.event(nEvs) = []; end
    if param.cellRippleOption data.RC.SWR.event(nEvs)     = []; end
  end
  dataOutSize(nEvs) = [];
  nEvs = nEvs - 1;
end

dataOut = (0: data.LFP.samplingInt : 2*param.swrWindow)';
nSignals = 1 + param.gammaOption + param.rOption;
if param.cellOption nSignals = nSignals + 1 + param.cellRawOption + param.cellGammaOption + param.cellRippleOption; end

% Output table names
tableVarNames = cell(1, 1 + nSignals*nEvs);
tableVarNames{1} = 'time';
nameInd = 2;

% Sort cell array into output table
for i = 1:nEvs
  
  % LFP event-locked data:
  dataOut = horzcat(dataOut, data.SWR.event{i});
  tableVarNames{nameInd} = ['LFP_' num2str(i)];
  nameInd = nameInd + 1;
  
  % Gamma event-locked data
  if param.gammaOption
    dataOut = horzcat(dataOut, data.gamma.SWR.event{i});
    tableVarNames{nameInd} = ['G_' num2str(i)];
    nameInd = nameInd + 1;
  end
  
  % Ripple event-locked data:
  if param.rOption
    dataOut = horzcat(dataOut, data.R.SWR.event{i});
    tableVarNames{nameInd} = ['R_' num2str(i)];
    nameInd = nameInd + 1;
  end
  
  if param.cellOption
    
    % cell event-locked data:
    if param.cellRawOption
      dataOut = horzcat(dataOut, data.C.SWR.event{i});
      tableVarNames{nameInd} = ['C_' num2str(i)];
      nameInd = nameInd + 1;
    end
    
    % cell (normalized) event-locked data:
    dataOut = horzcat(dataOut, data.C.SWR.evNorm{i});
    tableVarNames{nameInd} = ['CNorm_' num2str(i)];
    nameInd = nameInd + 1;
    
    % cell gamma-filtered event-locked data:
    if param.cellGammaOption
      dataOut = horzcat(dataOut, data.gammaC.SWR.event{i});
      tableVarNames{nameInd} = ['GC_' num2str(i)];
      nameInd = nameInd + 1;
    end
    
    % cell ripple-filtered  event-locked data:
    if param.cellRippleOption
      dataOut = horzcat(dataOut, data.RC.SWR.event{i});
      tableVarNames{nameInd} = ['RC_' num2str(i)];
      nameInd = nameInd + 1;
    end
  end
end

% Export data table
dataOutTable = array2table(dataOut, 'VariableNames', tableVarNames);
writetable(dataOutTable, expDataFile, 'Delimiter', '\t');

end

