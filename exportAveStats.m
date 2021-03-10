function exportAveStats(data, saveFile, exportFile)
%% exportAveStats(data, saveFile, exportFile)
%
%  Function to export csv file of averages of all LFP data available

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportAveStats');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_aveStats.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export average statistics', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('Please select valid file'); end
end

% Initialize Table:
outTable = table;

%% LFP
% LFP SWR event data
if isfield(data, 'SWR')
  varNames = {'nSWRs', 'SWR_frequency_Hz', 'SWR_duration_ms', 'SWR_IEI_s', 'SWR_amplitude_uV', 'SWR_area_uVs', 'SWR_power_uV2'};
  outTable = [outTable table(length(data.SWR.amp), data.SWR.frequency, mean(data.SWR.duration,'omitnan'), mean(data.SWR.IEI,'omitnan'), 10^3 * mean(data.SWR.amp,'omitnan'), mean(data.SWR.area,'omitnan'), 10^6 * mean(data.SWR.power,'omitnan'), 'VariableNames', varNames)];
end

% LFP Sharp Wave
if isfield(data, 'SW')
  % Total Stats:
  if isfield(data.SW, 'tPower') outTable = [outTable table(10^6 * data.SW.tPower, 'VariableNames', {'SW_Total_Power_uV2'})]; end
  % SWR Stats:
  if isfield(data.SW, 'SWR')
    if isfield(data.SW.SWR, 'power') outTable = [outTable table(10^6 * mean(data.SW.SWR.power,'omitnan'), 'VariableNames', {'SW_SWR_Power_uV2'})]; end
  end
end

% LFP Theta
if isfield(data, 'theta')
  % Total Stats:
  if isfield(data.theta, 'tPower') outTable = [outTable table(10^6 * data.theta.tPower, 'VariableNames', {'Theta_Tot_Power_uV2'})]; end
  if isfield(data.theta, 'FFT')
    if isfield(data.theta.FFT, 'pkFreq')  outTable = [outTable table(data.theta.FFT.pkFreq, 'VariableNames', {'Theta_FFT_pkFreq_Hz'})]; end
    if isfield(data.theta.FFT, 'fitMean') outTable = [outTable table(data.theta.FFT.fitMean, 'VariableNames', {'Theta_FFT_fitMean_Hz'})]; end
    if isfield(data.theta.FFT, 'fitSD')   outTable = [outTable table(data.theta.FFT.fitSD, 'VariableNames', {'Theta_FFT_fitSD_Hz'})]; end
    if isfield(data.theta.FFT, 'fitFWHM') outTable = [outTable table(data.theta.FFT.fitFWHM, 'VariableNames', {'Theta_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.theta, 'phase')
    if isfield(data.theta.phase, 'nCycle') outTable = [outTable table(data.theta.phase.nCycle, 'VariableNames', {'Theta_Tot_nCycle'})]; end
    if isfield(data.theta.phase, 'phFreq') outTable = [outTable table(data.theta.phase.phFreq, 'VariableNames', {'Theta_Tot_phFreq_Hz'})]; end
  end
end

% LFP Beta
if isfield(data, 'beta')
  % Total Stats:
  if isfield(data.beta, 'tPower') outTable = [outTable table(10^6 * data.beta.tPower, 'VariableNames', {'Beta_Tot_Power_uV2'})]; end
  if isfield(data.beta, 'FFT')
    if isfield(data.beta.FFT, 'pkFreq')  outTable = [outTable table(data.beta.FFT.pkFreq, 'VariableNames', {'Beta_FFT_pkFreq_Hz'})]; end
    if isfield(data.beta.FFT, 'fitMean') outTable = [outTable table(data.beta.FFT.fitMean, 'VariableNames', {'Beta_FFT_fitMean_Hz'})]; end
    if isfield(data.beta.FFT, 'fitSD')   outTable = [outTable table(data.beta.FFT.fitSD, 'VariableNames', {'Beta_FFT_fitSD_Hz'})]; end
    if isfield(data.beta.FFT, 'fitFWHM') outTable = [outTable table(data.beta.FFT.fitFWHM, 'VariableNames', {'Beta_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.beta, 'phase')
    if isfield(data.beta.phase, 'nCycle') outTable = [outTable table(data.beta.phase.nCycle, 'VariableNames', {'Beta_Tot_nCycle'})]; end
    if isfield(data.beta.phase, 'phFreq') outTable = [outTable table(data.beta.phase.phFreq, 'VariableNames', {'Beta_Tot_phFreq_Hz'})]; end
  end
end

% LFP Gamma
if isfield(data, 'gamma')
  % Total Stats:
  if isfield(data.gamma, 'tPower') outTable = [outTable table(10^6 * data.gamma.tPower, 'VariableNames', {'Gamma_Tot_Power_uV2'})]; end
  if isfield(data.gamma, 'FFT') 
    if isfield(data.gamma.FFT, 'pkFreq')  outTable = [outTable table(data.gamma.FFT.pkFreq, 'VariableNames', {'Gamma_FFT_pkFreq_Hz'})]; end
  	if isfield(data.gamma.FFT, 'fitMean') outTable = [outTable table(data.gamma.FFT.fitMean, 'VariableNames', {'Gamma_FFT_fitMean_Hz'})]; end
    if isfield(data.gamma.FFT, 'fitSD')   outTable = [outTable table(data.gamma.FFT.fitSD, 'VariableNames', {'Gamma_FFT_fitSD_Hz'})]; end
  	if isfield(data.gamma.FFT, 'fitFWHM') outTable = [outTable table(data.gamma.FFT.fitFWHM, 'VariableNames', {'Gamma_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.gamma, 'phase') 
    if isfield(data.gamma.phase, 'nCycle') outTable = [outTable table(data.gamma.phase.nCycle, 'VariableNames', {'Gamma_Tot_nCycle'})]; end
  	if isfield(data.gamma.phase, 'phFreq') outTable = [outTable table(data.gamma.phase.phFreq, 'VariableNames', {'Gamma_Tot_phFreq_Hz'})]; end
  end
  % SWR Stats:
  if isfield(data.gamma, 'SWR')
    if isfield(data.gamma.SWR, 'power') outTable = [outTable table(10^6 * mean(data.gamma.SWR.power,'omitnan'), 'VariableNames', {'Gamma_SWR_Power_uV2'})]; end
    if isfield(data.gamma.SWR, 'FFT')
      if isfield(data.gamma.SWR.FFT, 'pkFreq') outTable = [outTable table(mean(data.gamma.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'Gamma_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.gamma.SWR, 'phase')
      if isfield(data.gamma.SWR.phase, 'nCycle') outTable = [outTable table(mean(data.gamma.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'Gamma_SWR_nCycle'})]; end
      if isfield(data.gamma.SWR.phase, 'phFreq') outTable = [outTable table(mean(data.gamma.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'Gamma_SWR_phFreq_Hz'})]; end
    end
  end
end

% LFP High Gamma
if isfield(data, 'hgamma')
  % Total Stats:
  if isfield(data.hgamma, 'tPower') outTable = [outTable table(10^6 * data.hgamma.tPower, 'VariableNames', {'HGamma_Tot_Power_uV2'})]; end
  if isfield(data.hgamma, 'FFT') 
    if isfield(data.hgamma.FFT, 'pkFreq')  outTable = [outTable table(data.hgamma.FFT.pkFreq, 'VariableNames', {'HGamma_FFT_pkFreq_Hz'})]; end
  	if isfield(data.hgamma.FFT, 'fitMean') outTable = [outTable table(data.hgamma.FFT.fitMean, 'VariableNames', {'HGamma_FFT_fitMean_Hz'})]; end
    if isfield(data.hgamma.FFT, 'fitSD')   outTable = [outTable table(data.hgamma.FFT.fitSD, 'VariableNames', {'HGamma_FFT_fitSD_Hz'})]; end
  	if isfield(data.hgamma.FFT, 'fitFWHM') outTable = [outTable table(data.hgamma.FFT.fitFWHM, 'VariableNames', {'HGamma_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.hgamma, 'phase') 
    if isfield(data.hgamma.phase, 'nCycle') outTable = [outTable table(data.hgamma.phase.nCycle, 'VariableNames', {'HGamma_Tot_nCycle'})]; end
  	if isfield(data.hgamma.phase, 'phFreq') outTable = [outTable table(data.hgamma.phase.phFreq, 'VariableNames', {'HGamma_Tot_phFreq_Hz'})]; end
  end
  % SWR Stats:
  if isfield(data.hgamma, 'SWR')
    if isfield(data.hgamma.SWR, 'power') outTable = [outTable table(10^6 * mean(data.hgamma.SWR.power,'omitnan'), 'VariableNames', {'HGamma_SWR_Power_uV2'})]; end
    if isfield(data.hgamma.SWR, 'FFT')
      if isfield(data.hgamma.SWR.FFT, 'pkFreq') outTable = [outTable table(mean(data.hgamma.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'HGamma_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.hgamma.SWR, 'phase')
      if isfield(data.hgamma.SWR.phase, 'nCycle') outTable = [outTable table(mean(data.hgamma.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'HGamma_SWR_nCycle'})]; end
      if isfield(data.hgamma.SWR.phase, 'phFreq') outTable = [outTable table(mean(data.hgamma.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'HGamma_SWR_phFreq_Hz'})]; end
    end
  end
end

% LFP Ripple
if isfield(data, 'R')
  % Total Stats:
  if isfield(data.R, 'tPower') outTable = [outTable table(10^6 * data.R.tPower, 'VariableNames', {'Ripple_Tot_Power_uV2'})]; end
  % SWR Stats:
  if isfield(data.R, 'SWR')
    if isfield(data.R.SWR, 'power') outTable = [outTable table(10^6 * mean(data.R.SWR.power,'omitnan'), 'VariableNames', {'Ripple_SWR_Power_uV2'})]; end
    if isfield(data.R.SWR, 'FFT')
      if isfield(data.R.SWR.FFT, 'pkFreq') outTable = [outTable table(mean(data.R.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'Ripple_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.R.SWR, 'phase')
      if isfield(data.R.SWR.phase, 'nCycle') outTable = [outTable table(mean(data.R.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'Ripple_SWR_nCycle'})]; end
      if isfield(data.R.SWR.phase, 'phFreq') outTable = [outTable table(mean(data.R.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'Ripple_SWR_phFreq_Hz'})]; end
    end
  end
end

% LFP Fast Ripple
if isfield(data, 'fR')
  % Total Stats:
  if isfield(data.fR, 'tPower') outTable = [outTable table(10^6 * data.fR.tPower, 'VariableNames', {'fastRipple_Tot_Power_uV2'})]; end
  % SWR Stats:
  if isfield(data.fR, 'SWR')
    if isfield(data.fR.SWR, 'power') outTable = [outTable table(10^6 * mean(data.fR.SWR.power,'omitnan'), 'VariableNames', {'fastRipple_SWR_Power_uV2'})]; end
    if isfield(data.fR.SWR, 'FFT')
      if isfield(data.fR.SWR.FFT, 'pkFreq') outTable = [outTable table(mean(data.fR.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'fastRipple_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.fR.SWR, 'phase')
      if isfield(data.fR.SWR.phase, 'nCycle') outTable = [outTable table(mean(data.fR.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'fastRipple_SWR_nCycle'})]; end
      if isfield(data.fR.SWR.phase, 'phFreq') outTable = [outTable table(mean(data.fR.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'fastRipple_SWR_phFreq_Hz'})]; end
    end
  end
end

%% Cell-Attached
% if isfield(data, 'C')
%   
%   
% end

% Convert to cell array and replace NaN values with blanks:
varNames = outTable.Properties.VariableNames;
tmpCell  = table2cell(outTable);
tmpCell(isnan(outTable.Variables)) = {[]};

% Convert to array and transpose table:
tmpArray = cell2mat(tmpCell);
outTable = array2table(tmpArray.');
outTable.Properties.RowNames = varNames;

writetable(outTable, exportFile, 'Delimiter', ',', 'WriteVariableNames', false, 'WriteRowNames', true);

end