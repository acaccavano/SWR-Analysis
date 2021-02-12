function exportLFPAverages(data, saveFile, exportFile)
%% exportLFPAverages(data, saveFile, exportFile)
%
%  Function to export csv file of averages of all LFP data available

% Handle input arguments - if not entered
if (nargin < 3) exportFile = []; end
if (nargin < 2) saveFile   = []; end
if (nargin < 1) data       = []; end

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportSpkEvents');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_lfpAve.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export averages', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile) error('Please select valid file'); end
end

% LFP SWR event data
varNames = {'nSWRs', 'SWR_frequency_Hz', 'SWR_duration_ms', 'SWR_IEI_s', 'SWR_amplitude_uV'};
outTable = table(length(data.SWR.amp), data.SWR.frequency, mean(data.SWR.duration,'omitnan'), mean(data.SWR.IEI,'omitnan'), 10^3 * mean(data.SWR.amp,'omitnan'), 'VariableNames', varNames);
if isfield(data.SWR, 'area') outTable = [outTable table(mean(data.SWR.area,'omitnan'), 'VariableNames', {'SWR_area_uVs'})]; end
if isfield(data.SWR, 'power') outTable = [outTable table(10^6 * mean(data.SWR.power,'omitnan'), 'VariableNames', {'SWR_power_uV2'})]; end

% LFP Sharp Wave
if isfield(data, 'SW')
  if isfield(data.SW, 'tPower') outTable = [outTable table(10^6 * data.SW.tPower, 'VariableNames', {'SW_Total_Power_uV2'})]; end
  if isfield(data.SW, 'SWR')
    if isfield(data.SW.SWR, 'power') outTable = [outTable table(10^6 * mean(data.SW.SWR.power,'omitnan'), 'VariableNames', {'SW_SWR_Power_uV2'})]; end
  end
end

% LFP Theta
if isfield(data, 'theta')
  if isfield(data.theta, 'tPower') outTable = [outTable table(10^6 * data.theta.tPower, 'VariableNames', {'Theta_Total_Power_uV2'})]; end
end

% LFP Theta
if isfield(data, 'beta')
  if isfield(data.beta, 'tPower') outTable = [outTable table(10^6 * data.beta.tPower, 'VariableNames', {'Beta_Total_Power_uV2'})]; end
end

% LFP Gamma
if isfield(data, 'gamma')
  if isfield(data.gamma, 'tPower') outTable = [outTable table(10^6 * data.gamma.tPower, 'VariableNames', {'Gamma_Total_Power_uV2'})]; end
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
  if isfield(data.hgamma, 'tPower') outTable = [outTable table(10^6 * data.hgamma.tPower, 'VariableNames', {'HGamma_Total_Power_uV2'})]; end
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
  if isfield(data.R, 'tPower') outTable = [outTable table(10^6 * data.R.tPower, 'VariableNames', {'Ripple_Total_Power_uV2'})]; end
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
  if isfield(data.fR, 'tPower') outTable = [outTable table(10^6 * data.fR.tPower, 'VariableNames', {'fastRipple_Total_Power_uV2'})]; end
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

% Replace NaN values with blanks
tmp = table2cell(outTable);
tmp(isnan(outTable.Variables)) = {[]};
outTable = array2table(tmp,'VariableNames',outTable.Properties.VariableNames);

writetable(outTable, exportFile, 'Delimiter', ',');

end