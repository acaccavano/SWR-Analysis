function exportSWREvents(data, saveFile, exportFile)

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_swrEvents.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export events', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('Please select valid file'); end
end

% LFP SWR event data
data.SWR.IEI = vertcat(0, data.SWR.IEI); % Always one less than other parameters, add zero for 1st event
varNames = {'evStart', 'evPeak', 'evEnd', 'duration_ms', 'IEI_s', 'amplitude_uV'};
outTable = table(data.SWR.evStart, data.SWR.evPeak, data.SWR.evEnd, data.SWR.duration, data.SWR.IEI, 10^3 * data.SWR.amp, 'VariableNames', varNames);
if isfield(data.SWR, 'area') outTable = [outTable table(data.SWR.area, 'VariableNames', {'area_uVs'})]; end
if isfield(data.SWR, 'power') outTable = [outTable table(10^6 * data.SWR.power, 'VariableNames', {'SWR_power_uV2'})]; end

% LFP Sharp Wave
if isfield(data, 'SW')
  if isfield(data.SW, 'SWR')
    if isfield(data.SW.SWR, 'power') outTable = [outTable table(10^6 * data.SW.SWR.power, 'VariableNames', {'SW_Power_uV2'})]; end
  end
end

% LFP Gamma
if isfield(data, 'gamma')
  if isfield(data.gamma, 'SWR')
    if isfield(data.gamma.SWR, 'power') outTable = [outTable table(10^6 * data.gamma.SWR.power, 'VariableNames', {'Gamma_Power_uV2'})]; end
    if isfield(data.gamma.SWR, 'FFT')
      if isfield(data.gamma.SWR.FFT, 'pkFreq') outTable = [outTable table(data.gamma.SWR.FFT.pkFreq, 'VariableNames', {'Gamma_pkFreq_Hz'})]; end
    end
    if isfield(data.gamma.SWR, 'phase')
      if isfield(data.gamma.SWR.phase, 'nCycle') outTable = [outTable table(data.gamma.SWR.phase.nCycle, 'VariableNames', {'Gamma_nCycle'})]; end
      if isfield(data.gamma.SWR.phase, 'phFreq') outTable = [outTable table(data.gamma.SWR.phase.phFreq, 'VariableNames', {'Gamma_phFreq_Hz'})]; end
    end
  end
end

% LFP Ripple
if isfield(data, 'R')
  if isfield(data.R, 'SWR')
    if isfield(data.R.SWR, 'power') outTable = [outTable table(10^6 * data.R.SWR.power, 'VariableNames', {'Ripple_Power_uV2'})]; end
    if isfield(data.R.SWR, 'FFT')
      if isfield(data.R.SWR.FFT, 'pkFreq') outTable = [outTable table(data.R.SWR.FFT.pkFreq, 'VariableNames', {'Ripple_pkFreq_Hz'})]; end
    end
    if isfield(data.R.SWR, 'phase')
      if isfield(data.R.SWR.phase, 'nCycle') outTable = [outTable table(data.R.SWR.phase.nCycle, 'VariableNames', {'Ripple_nCycle'})]; end
      if isfield(data.R.SWR.phase, 'phFreq') outTable = [outTable table(data.R.SWR.phase.phFreq, 'VariableNames', {'Ripple_phFreq_Hz'})]; end
    end
  end
end

% LFP Fast Ripple
if isfield(data, 'fR')
  if isfield(data.fR, 'SWR')
    if isfield(data.fR.SWR, 'power') outTable = [outTable table(10^6 * data.fR.SWR.power, 'VariableNames', {'fastRipple_Power_uV2'})]; end
    if isfield(data.fR.SWR, 'FFT')
      if isfield(data.fR.SWR.FFT, 'pkFreq') outTable = [outTable table(data.fR.SWR.FFT.pkFreq, 'VariableNames', {'fastRipple_pkFreq_Hz'})]; end
    end
    if isfield(data.fR.SWR, 'phase')
      if isfield(data.fR.SWR.phase, 'nCycle') outTable = [outTable table(data.fR.SWR.phase.nCycle, 'VariableNames', {'fastRipple_nCycle'})]; end
      if isfield(data.fR.SWR.phase, 'phFreq') outTable = [outTable table(data.fR.SWR.phase.phFreq, 'VariableNames', {'fastRipple_phFreq_Hz'})]; end
    end
  end
end

if isfield(data, 'C')
  if isfield(data.C, 'SWR')
    if isfield(data.C.SWR, 'baseline') outTable = [outTable table(data.C.SWR.baseline, 'VariableNames', {'cellBaseline_pA'})]; end
    if isfield(data.C.SWR, 'area')     outTable = [outTable table(10^-3 * data.C.SWR.area, 'VariableNames', {'cellArea_pC'})]; end
    if isfield(data.C.SWR, 'areaQ1')   outTable = [outTable table(10^-3 * data.C.SWR.areaQ1, 'VariableNames', {'cellAreaQ1_pC'})]; end
    if isfield(data.C.SWR, 'areaQ2')   outTable = [outTable table(10^-3 * data.C.SWR.areaQ2, 'VariableNames', {'cellAreaQ2_pC'})]; end
    if isfield(data.C.SWR, 'areaQ3')   outTable = [outTable table(10^-3 * data.C.SWR.areaQ3, 'VariableNames', {'cellAreaQ3_pC'})]; end
    if isfield(data.C.SWR, 'areaQ4')   outTable = [outTable table(10^-3 * data.C.SWR.areaQ4, 'VariableNames', {'cellAreaQ4_pC'})]; end
  end
end


% Cell channel gamma
if isfield(data, 'gammaC')
  if isfield(data.gammaC, 'SWR')
    if isfield(data.gammaC.SWR, 'power') outTable = [outTable table(data.gammaC.SWR.power, 'VariableNames', {'cellGamma_Power_pA2'})]; end
    if isfield(data.gammaC.SWR, 'FFT')
      if isfield(data.gammaC.SWR.FFT, 'pkFreq') outTable = [outTable table(data.gammaC.SWR.FFT.pkFreq, 'VariableNames', {'cellGamma_pkFreq_Hz'})]; end
    end
    if isfield(data.gammaC.SWR, 'phase')
      if isfield(data.gammaC.SWR.phase, 'nCycle') outTable = [outTable table(data.gammaC.SWR.phase.nCycle, 'VariableNames', {'cellGamma_nCycle'})]; end
      if isfield(data.gammaC.SWR.phase, 'phFreq') outTable = [outTable table(data.gammaC.SWR.phase.phFreq, 'VariableNames', {'cellGamma_phFreq_Hz'})]; end
    end
  end
end

% LFP-Cell Gamma Correlation
if isfield(data, 'gamma') && isfield(data, 'gammaC')
  if isfield(data.gammaC, 'SWR')
    if isfield(data.gammaC.SWR, 'oCorr') 
      outTable = [outTable table(data.gammaC.SWR.oCorr(:,1), data.gammaC.SWR.minCorr(:,1), data.gammaC.SWR.minCorr(:,2), data.gammaC.SWR.maxCorr(:,1), data.gammaC.SWR.maxCorr(:,2), 'VariableNames', {'Gamma_Corr', 'Gamma_minCorr', 'Gamma_minTime_ms', 'Gamma_maxCorr', 'Gamma_maxTime_ms'})]; 
    end
  end
end

% Cell channel ripple
if isfield(data, 'RC')
  if isfield(data.RC, 'SWR')
    if isfield(data.RC.SWR, 'power') outTable = [outTable table(data.RC.SWR.power, 'VariableNames', {'cellRipple_Power_pA2'})]; end
    if isfield(data.RC.SWR, 'FFT')
      if isfield(data.RC.SWR.FFT, 'pkFreq') outTable = [outTable table(data.RC.SWR.FFT.pkFreq, 'VariableNames', {'cellRipple_pkFreq_Hz'})]; end
    end
    if isfield(data.RC.SWR, 'phase')
      if isfield(data.RC.SWR.phase, 'nCycle') outTable = [outTable table(data.RC.SWR.phase.nCycle, 'VariableNames', {'cellRipple_nCycle'})]; end
      if isfield(data.RC.SWR.phase, 'phFreq') outTable = [outTable table(data.RC.SWR.phase.phFreq, 'VariableNames', {'cellRipple_phFreq_Hz'})]; end
    end
  end
end

% LFP-Cell Ripple Correlation
if isfield(data, 'R') && isfield(data, 'RC')
  if isfield(data.RC, 'SWR')
    if isfield(data.RC.SWR, 'oCorr') 
      outTable = [outTable table(data.RC.SWR.oCorr(:,1), data.RC.SWR.minCorr(:,1), data.RC.SWR.minCorr(:,2), data.RC.SWR.maxCorr(:,1), data.RC.SWR.maxCorr(:,2), 'VariableNames', {'Ripple_Corr', 'Ripple_minCorr', 'Ripple_minTime_ms', 'Ripple_maxCorr', 'Ripple_maxTime_ms'})]; 
    end
  end
end

% If SWR-Spk correlation has been performed:
if isfield(data.SWR, 'spike')
  if (length(data.SWR.evStart) > length(data.SWR.spike.evMatrix))
    outTable = [outTable table(vertcat(data.SWR.spike.evMatrix, NaN * (1: length(data.SWR.evStart) - length(data.SWR.spike.evMatrix))'), 'VariableNames', {'swrSpkC'})];
  else
    outTable = [outTable table(data.SWR.spike.evMatrix, 'VariableNames', {'swrSpkC'})];
  end
end

% If SWR-Bst correlation has been performed:
if isfield(data.SWR, 'burst')
  if (length(data.SWR.evStart) > length(data.SWR.burst.evMatrix))
    outTable = [outTable table(vertcat(data.SWR.burst.evMatrix, NaN * (1: length(data.SWR.evStart) - length(data.SWR.burst.evMatrix))'), 'VariableNames', {'swrBstC'})];
  else
    outTable = [outTable table(data.SWR.burst.evMatrix, 'VariableNames', {'swrBstC'})];
  end
end

% If SWR-PSC correlation has been performed:
if isfield(data.SWR, 'PSC')
  if (length(data.SWR.evStart) > length(data.SWR.PSC.evMatrix))
    outTable = [outTable table(vertcat(data.SWR.PSC.evMatrix, NaN * (1: length(data.SWR.evStart) - length(data.SWR.PSC.evMatrix))'), 'VariableNames', {'swrPSCC'})];
  else
    outTable = [outTable table(data.SWR.PSC.evMatrix, 'VariableNames', {'swrPSCC'})];
  end
end

% If SWR-Ca correlation has been performed:
if isfield(data.SWR, 'Ca')
  if (length(data.SWR.evStart) > data.SWR.Ca.nEventsA)
    outTable = [outTable table(vertcat(data.SWR.Ca.nCellsC, NaN * (1: length(data.SWR.evStart) - data.SWR.Ca.nEventsA)'), 'VariableNames', {'nCellsCaC'})];
  else
    outTable = [outTable table(data.SWR.Ca.nCellsC, 'VariableNames', {'nCellsCaC'})];
  end
end

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end