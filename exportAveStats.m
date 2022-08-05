function exportAveStats(data, saveFile, exportFile)
%% exportAveStats(data, saveFile, exportFile)
%
%  Function to export csv file of averages of all data available

% Handle input arguments - if not entered
if (nargin < 3); exportFile = []; end
if (nargin < 2); saveFile   = []; end
if (nargin < 1); data       = []; end

transposeOption = false; % Set to true to get one data column with headings as row titles

if isempty(data) || isempty(saveFile)
  error('Enter sufficient inputs to use function exportAveStats');
end

if isempty(exportFile)
  [parentPath, saveFileName, ~] = parsePath(saveFile);
  defaultName = [parentPath saveFileName '_aveStats.csv'];
  [exportName, exportPath] = uiputfile('.csv','Select *.csv file to export average statistics', defaultName);
  exportFile = [exportPath exportName];
  if ~all(exportFile); error('Please select valid file'); end
end

% Initialize Table:
outTable = table;

%% LFP
% LFP SWR event data
if isfield(data, 'SWR')
  varNames = {'nSWRs', 'SWR_frequency_Hz', 'SWR_duration_ms', 'SWR_IEI_s', 'SWR_amplitude_uV', 'SWR_area_uVs', 'SWR_power_uV2'};
  outTable = [outTable table(length(data.SWR.evStart), data.SWR.frequency, mean(data.SWR.duration,'omitnan'), mean(data.SWR.IEI,'omitnan'), 10^3 * mean(data.SWR.amp,'omitnan'), mean(data.SWR.area,'omitnan'), 10^6 * mean(data.SWR.power,'omitnan'), 'VariableNames', varNames)];
end

% LFP Sharp Wave
if isfield(data, 'SW')
  
  % SW Event Stats:
  varNames = {'nSWs', 'SW_frequency_Hz', 'SW_duration_ms', 'SW_IEI_s', 'SW_Event_power_uV2'};
  outTable = [outTable table(length(data.SW.evStart), data.SW.frequency, mean(data.SW.duration,'omitnan'), mean(data.SW.IEI,'omitnan'), 10^6 * mean(data.SW.power,'omitnan'), 'VariableNames', varNames)];
  
  % Total Stats:
  if isfield(data.SW, 'tPower'); outTable = [outTable table(10^6 * data.SW.tPower, 'VariableNames', {'SW_Total_Power_uV2'})]; end
  
  % SWR Stats:
  if isfield(data.SW, 'SWR')
    if isfield(data.SW.SWR, 'power'); outTable = [outTable table(10^6 * mean(data.SW.SWR.power,'omitnan'), 'VariableNames', {'SW_SWR_Power_uV2'})]; end
  end
end

% LFP Theta
if isfield(data, 'theta')
  
  % Total Stats:
  if isfield(data.theta, 'tPower'); outTable = [outTable table(10^6 * data.theta.tPower, 'VariableNames', {'Theta_Tot_Power_uV2'})]; end
  if isfield(data.theta, 'FFT')
    if isfield(data.theta.FFT, 'pkFreq');  outTable = [outTable table(data.theta.FFT.pkFreq, 'VariableNames', {'Theta_FFT_pkFreq_Hz'})]; end
    if isfield(data.theta.FFT, 'fitMean'); outTable = [outTable table(data.theta.FFT.fitMean, 'VariableNames', {'Theta_FFT_fitMean_Hz'})]; end
    if isfield(data.theta.FFT, 'fitSD');   outTable = [outTable table(data.theta.FFT.fitSD, 'VariableNames', {'Theta_FFT_fitSD_Hz'})]; end
    if isfield(data.theta.FFT, 'fitFWHM'); outTable = [outTable table(data.theta.FFT.fitFWHM, 'VariableNames', {'Theta_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.theta, 'phase')
    if isfield(data.theta.phase, 'nCycle'); outTable = [outTable table(data.theta.phase.nCycle, 'VariableNames', {'Theta_Tot_nCycle'})]; end
    if isfield(data.theta.phase, 'phFreq'); outTable = [outTable table(data.theta.phase.phFreq, 'VariableNames', {'Theta_Tot_phFreq_Hz'})]; end
  end
end

% LFP Beta
if isfield(data, 'beta')
  
  % Total Stats:
  if isfield(data.beta, 'tPower'); outTable = [outTable table(10^6 * data.beta.tPower, 'VariableNames', {'Beta_Tot_Power_uV2'})]; end
  if isfield(data.beta, 'FFT')
    if isfield(data.beta.FFT, 'pkFreq');  outTable = [outTable table(data.beta.FFT.pkFreq, 'VariableNames', {'Beta_FFT_pkFreq_Hz'})]; end
    if isfield(data.beta.FFT, 'fitMean'); outTable = [outTable table(data.beta.FFT.fitMean, 'VariableNames', {'Beta_FFT_fitMean_Hz'})]; end
    if isfield(data.beta.FFT, 'fitSD');   outTable = [outTable table(data.beta.FFT.fitSD, 'VariableNames', {'Beta_FFT_fitSD_Hz'})]; end
    if isfield(data.beta.FFT, 'fitFWHM'); outTable = [outTable table(data.beta.FFT.fitFWHM, 'VariableNames', {'Beta_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.beta, 'phase')
    if isfield(data.beta.phase, 'nCycle'); outTable = [outTable table(data.beta.phase.nCycle, 'VariableNames', {'Beta_Tot_nCycle'})]; end
    if isfield(data.beta.phase, 'phFreq'); outTable = [outTable table(data.beta.phase.phFreq, 'VariableNames', {'Beta_Tot_phFreq_Hz'})]; end
  end
end

% LFP Gamma
if isfield(data, 'gamma')
  
  % Total Stats:
  if isfield(data.gamma, 'tPower'); outTable = [outTable table(10^6 * data.gamma.tPower, 'VariableNames', {'Gamma_Tot_Power_uV2'})]; end
  if isfield(data.gamma, 'FFT')
    if isfield(data.gamma.FFT, 'pkFreq');  outTable = [outTable table(data.gamma.FFT.pkFreq, 'VariableNames', {'Gamma_FFT_pkFreq_Hz'})]; end
    if isfield(data.gamma.FFT, 'fitMean'); outTable = [outTable table(data.gamma.FFT.fitMean, 'VariableNames', {'Gamma_FFT_fitMean_Hz'})]; end
    if isfield(data.gamma.FFT, 'fitSD');   outTable = [outTable table(data.gamma.FFT.fitSD, 'VariableNames', {'Gamma_FFT_fitSD_Hz'})]; end
    if isfield(data.gamma.FFT, 'fitFWHM'); outTable = [outTable table(data.gamma.FFT.fitFWHM, 'VariableNames', {'Gamma_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.gamma, 'phase')
    if isfield(data.gamma.phase, 'nCycle'); outTable = [outTable table(data.gamma.phase.nCycle, 'VariableNames', {'Gamma_Tot_nCycle'})]; end
    if isfield(data.gamma.phase, 'phFreq'); outTable = [outTable table(data.gamma.phase.phFreq, 'VariableNames', {'Gamma_Tot_phFreq_Hz'})]; end
  end
  if isfield(data.gamma, 'xFreq')
    if isfield(data.gamma.xFreq, 'pacMI_Len');       outTable = [outTable table(data.gamma.xFreq.pacMI_Len, 'VariableNames', {'Gamma_pacMI_Len'})]; end
    if isfield(data.gamma.xFreq, 'pacMI_Phase');     outTable = [outTable table(data.gamma.xFreq.pacMI_Phase, 'VariableNames', {'Gamma_pacMI_Phase'})]; end
    if isfield(data.gamma.xFreq, 'pacMI_LenZ');      outTable = [outTable table(data.gamma.xFreq.pacMI_LenZ, 'VariableNames', {'Gamma_pacMI_LenZ'})]; end
    if isfield(data.gamma.xFreq, 'pacMIShf_LenAve'); outTable = [outTable table(data.gamma.xFreq.pacMIShf_LenAve, 'VariableNames', {'Gamma_pacMIShf_LenAve'})]; end
    if isfield(data.gamma.xFreq, 'pacMIShf_LenSTD'); outTable = [outTable table(data.gamma.xFreq.pacMIShf_LenSTD, 'VariableNames', {'Gamma_pacMIShf_LenSTD'})]; end
  end
  
  % SWR Stats:
  if isfield(data.gamma, 'SWR')
    if isfield(data.gamma.SWR, 'power'); outTable = [outTable table(10^6 * mean(data.gamma.SWR.power,'omitnan'), 'VariableNames', {'Gamma_SWR_Power_uV2'})]; end
    if isfield(data.gamma.SWR, 'FFT')
      if isfield(data.gamma.SWR.FFT, 'pkFreq'); outTable = [outTable table(mean(data.gamma.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'Gamma_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.gamma.SWR, 'phase')
      if isfield(data.gamma.SWR.phase, 'nCycle'); outTable = [outTable table(mean(data.gamma.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'Gamma_SWR_nCycle'})]; end
      if isfield(data.gamma.SWR.phase, 'phFreq'); outTable = [outTable table(mean(data.gamma.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'Gamma_SWR_phFreq_Hz'})]; end
    end
  end
end

% LFP High Gamma
if isfield(data, 'hgamma')
  
  % Total Stats:
  if isfield(data.hgamma, 'tPower'); outTable = [outTable table(10^6 * data.hgamma.tPower, 'VariableNames', {'HGamma_Tot_Power_uV2'})]; end
  if isfield(data.hgamma, 'FFT') 
    if isfield(data.hgamma.FFT, 'pkFreq');  outTable = [outTable table(data.hgamma.FFT.pkFreq, 'VariableNames', {'HGamma_FFT_pkFreq_Hz'})]; end
  	if isfield(data.hgamma.FFT, 'fitMean'); outTable = [outTable table(data.hgamma.FFT.fitMean, 'VariableNames', {'HGamma_FFT_fitMean_Hz'})]; end
    if isfield(data.hgamma.FFT, 'fitSD');   outTable = [outTable table(data.hgamma.FFT.fitSD, 'VariableNames', {'HGamma_FFT_fitSD_Hz'})]; end
  	if isfield(data.hgamma.FFT, 'fitFWHM'); outTable = [outTable table(data.hgamma.FFT.fitFWHM, 'VariableNames', {'HGamma_FFT_fitFWHM_Hz'})]; end
  end
  if isfield(data.hgamma, 'phase') 
    if isfield(data.hgamma.phase, 'nCycle'); outTable = [outTable table(data.hgamma.phase.nCycle, 'VariableNames', {'HGamma_Tot_nCycle'})]; end
  	if isfield(data.hgamma.phase, 'phFreq'); outTable = [outTable table(data.hgamma.phase.phFreq, 'VariableNames', {'HGamma_Tot_phFreq_Hz'})]; end
  end
  if isfield(data.hgamma, 'xFreq')
    if isfield(data.hgamma.xFreq, 'pacMI_Len');       outTable = [outTable table(data.hgamma.xFreq.pacMI_Len, 'VariableNames', {'HGamma_pacMI_Len'})]; end
    if isfield(data.hgamma.xFreq, 'pacMI_Phase');     outTable = [outTable table(data.hgamma.xFreq.pacMI_Phase, 'VariableNames', {'HGamma_pacMI_Phase'})]; end
    if isfield(data.hgamma.xFreq, 'pacMI_LenZ');      outTable = [outTable table(data.hgamma.xFreq.pacMI_LenZ, 'VariableNames', {'HGamma_pacMI_LenZ'})]; end
    if isfield(data.hgamma.xFreq, 'pacMIShf_LenAve'); outTable = [outTable table(data.hgamma.xFreq.pacMIShf_LenAve, 'VariableNames', {'HGamma_pacMIShf_LenAve'})]; end
    if isfield(data.hgamma.xFreq, 'pacMIShf_LenSTD'); outTable = [outTable table(data.hgamma.xFreq.pacMIShf_LenSTD, 'VariableNames', {'HGamma_pacMIShf_LenSTD'})]; end
  end
  
  % SWR Stats:
  if isfield(data.hgamma, 'SWR')
    if isfield(data.hgamma.SWR, 'power'); outTable = [outTable table(10^6 * mean(data.hgamma.SWR.power,'omitnan'), 'VariableNames', {'HGamma_SWR_Power_uV2'})]; end
    if isfield(data.hgamma.SWR, 'FFT')
      if isfield(data.hgamma.SWR.FFT, 'pkFreq'); outTable = [outTable table(mean(data.hgamma.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'HGamma_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.hgamma.SWR, 'phase')
      if isfield(data.hgamma.SWR.phase, 'nCycle'); outTable = [outTable table(mean(data.hgamma.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'HGamma_SWR_nCycle'})]; end
      if isfield(data.hgamma.SWR.phase, 'phFreq'); outTable = [outTable table(mean(data.hgamma.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'HGamma_SWR_phFreq_Hz'})]; end
    end
  end
end

% LFP Ripple
if isfield(data, 'R')
  
  % R Event Stats:
  varNames = {'nRipples', 'Ripple_frequency_Hz', 'Ripple_duration_ms', 'Ripple_IEI_s', 'Ripple_Event_power_uV2'};
  outTable = [outTable table(length(data.R.evStart), data.R.frequency, mean(data.R.duration,'omitnan'), mean(data.R.IEI,'omitnan'), 10^6 * mean(data.R.power,'omitnan'), 'VariableNames', varNames)];

  % Total Stats:
  if isfield(data.R, 'tPower'); outTable = [outTable table(10^6 * data.R.tPower, 'VariableNames', {'Ripple_Tot_Power_uV2'})]; end
  
  % SWR Stats:
  if isfield(data.R, 'SWR')
    if isfield(data.R.SWR, 'power'); outTable = [outTable table(10^6 * mean(data.R.SWR.power,'omitnan'), 'VariableNames', {'Ripple_SWR_Power_uV2'})]; end
    if isfield(data.R.SWR, 'FFT')
      if isfield(data.R.SWR.FFT, 'pkFreq'); outTable = [outTable table(mean(data.R.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'Ripple_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.R.SWR, 'phase')
      if isfield(data.R.SWR.phase, 'nCycle'); outTable = [outTable table(mean(data.R.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'Ripple_SWR_nCycle'})]; end
      if isfield(data.R.SWR.phase, 'phFreq'); outTable = [outTable table(mean(data.R.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'Ripple_SWR_phFreq_Hz'})]; end
    end
  end
  
end

% LFP Fast Ripple
if isfield(data, 'fR')
  
  % Total Stats:
  if isfield(data.fR, 'tPower'); outTable = [outTable table(10^6 * data.fR.tPower, 'VariableNames', {'fastRipple_Tot_Power_uV2'})]; end
  
  % SWR Stats:
  if isfield(data.fR, 'SWR')
    if isfield(data.fR.SWR, 'power'); outTable = [outTable table(10^6 * mean(data.fR.SWR.power,'omitnan'), 'VariableNames', {'fastRipple_SWR_Power_uV2'})]; end
    if isfield(data.fR.SWR, 'FFT')
      if isfield(data.fR.SWR.FFT, 'pkFreq'); outTable = [outTable table(mean(data.fR.SWR.FFT.pkFreq,'omitnan'), 'VariableNames', {'fastRipple_SWR_pkFreq_Hz'})]; end
    end
    if isfield(data.fR.SWR, 'phase')
      if isfield(data.fR.SWR.phase, 'nCycle'); outTable = [outTable table(mean(data.fR.SWR.phase.nCycle,'omitnan'), 'VariableNames', {'fastRipple_SWR_nCycle'})]; end
      if isfield(data.fR.SWR.phase, 'phFreq'); outTable = [outTable table(mean(data.fR.SWR.phase.phFreq,'omitnan'), 'VariableNames', {'fastRipple_SWR_phFreq_Hz'})]; end
    end
  end
end

%% Cell-Attached
if isfield(data, 'C')
  
  % Spike Stats:
  if isfield(data.C, 'spike')
    
    % Total Stats:
    if isfield(data.C.spike, 'nEvents'); outTable = [outTable table(data.C.spike.nEvents, 'VariableNames', {'nSpikes'})]; end
    if isfield(data.C.spike, 'frequency'); outTable = [outTable table(data.C.spike.frequency, 'VariableNames', {'Spike_frequency_Hz'})]; end
    
    % SWR Stats:
    if isfield(data.C.spike, 'nEventsA'); outTable = [outTable table(data.C.spike.nEventsA, 'VariableNames', {'nSpikes_Align'})]; end
    if isfield(data.C.spike, 'nEventsC'); outTable = [outTable table(data.C.spike.nEventsC, 'VariableNames', {'nSpikes_Coinc'})]; end
    
    % Theta Stats:
    if isfield(data.C.spike, 'theta')
      if isfield(data.C.spike.theta, 'phaseAve') 
        varNames = {'Spike-Theta_Ave_rad', 'Spike-Theta_R', 'Spike-Theta_P', 'Spike-Theta_Z'};
        outTable = [outTable table(data.C.spike.theta.phaseAve, data.C.spike.theta.phaseR, data.C.spike.theta.phaseP, data.C.spike.theta.phaseZ, 'VariableNames', varNames)];
      end
    end
    
    % Beta Stats:
    if isfield(data.C.spike, 'beta')
      if isfield(data.C.spike.beta, 'phaseAve') 
        varNames = {'Spike-Beta_Ave_rad', 'Spike-Beta_R', 'Spike-Beta_P', 'Spike-Beta_Z'};
        outTable = [outTable table(data.C.spike.beta.phaseAve, data.C.spike.beta.phaseR, data.C.spike.beta.phaseP, data.C.spike.beta.phaseZ, 'VariableNames', varNames)];
      end
    end
    
    % Gamma Stats:
    if isfield(data.C.spike, 'gamma')
      if isfield(data.C.spike.gamma, 'phaseAve') 
        varNames = {'Spike-Gamma_Ave_rad', 'Spike-Gamma_R', 'Spike-Gamma_P', 'Spike-Gamma_Z'};
        outTable = [outTable table(data.C.spike.gamma.phaseAve, data.C.spike.gamma.phaseR, data.C.spike.gamma.phaseP, data.C.spike.gamma.phaseZ, 'VariableNames', varNames)];
      end
    end
    
    % High Gamma Stats:
    if isfield(data.C.spike, 'hgamma')
      if isfield(data.C.spike.hgamma, 'phaseAve') 
        varNames = {'Spike-highGamma_Ave_rad', 'Spike-highGamma_R', 'Spike-highGamma_P', 'Spike-highGamma_Z'};
        outTable = [outTable table(data.C.spike.hgamma.phaseAve, data.C.spike.hgamma.phaseR, data.C.spike.hgamma.phaseP, data.C.spike.hgamma.phaseZ, 'VariableNames', varNames)];
      end
    end

    % Ripple Stats:
    if isfield(data.C.spike, 'R')
      if isfield(data.C.spike.R, 'phaseAve') 
        varNames = {'Spike-Ripple_Ave_rad', 'Spike-Ripple_R', 'Spike-Ripple_P', 'Spike-Ripple_Z'};
        outTable = [outTable table(data.C.spike.R.phaseAve, data.C.spike.R.phaseR, data.C.spike.R.phaseP, data.C.spike.R.phaseZ, 'VariableNames', varNames)];
      end
    end
    
    % Fast Ripple Stats:
    if isfield(data.C.spike, 'fR')
      if isfield(data.C.spike.fR, 'phaseAve') 
        varNames = {'Spike-fastRipple_Ave_rad', 'Spike-fastRipple_R', 'Spike-fastRipple_P', 'Spike-fastRipple_Z'};
        outTable = [outTable table(data.C.spike.fR.phaseAve, data.C.spike.fR.phaseR, data.C.spike.fR.phaseP, data.C.spike.fR.phaseZ, 'VariableNames', varNames)];
      end
    end
  end
  
  % Burst Stats:
  if isfield(data.C, 'burst')
    % SWR Stats:
    if isfield(data.C.burst, 'nEventsA'); outTable = [outTable table(data.C.burst.nEventsA, 'VariableNames', {'nBursts_Align'})]; end
    if isfield(data.C.burst, 'nEventsC'); outTable = [outTable table(data.C.burst.nEventsC, 'VariableNames', {'nBursts_Coinc'})]; end
    if isfield(data.C.burst, 'nSpike'); outTable = [outTable table(mean(data.C.burst.nSpike,'omitnan'), 'VariableNames', {'nSpikesinBurst'})]; end
    if isfield(data.C.burst, 'intraBI'); outTable = [outTable table(mean(data.C.burst.intraBI,'omitnan'), 'VariableNames', {'intraBurstInt_ms'})]; end
  end
end

%% Calcium
if isfield(data, 'Ca')
  varNames = {'nCells', 'nCellsActive', 'nEvents', 'frequency_Hz', 'aveIEI_s', 'aveAmplitude_dFoF', 'aveDuration_s', 'aveArea_dFoFs'};
  outTable = [outTable table(mean(data.Ca.nEvents,'omitnan'), mean(data.Ca.frequency,'omitnan'), mean(data.Ca.IEIAve,'omitnan'), mean(data.Ca.ampAve,'omitnan'), mean(data.Ca.durAve/1000,'omitnan'), mean(data.Ca.areaAve/1000,'omitnan'), 'VariableNames', varNames)];
end

% Convert to cell array and replace NaN values with blanks:
varNames = outTable.Properties.VariableNames;
tmpCell  = table2cell(outTable);
tmpCell(isnan(outTable.Variables)) = {[]};

% Convert to array and transpose table:
if transposeOption
  tmpArray = cell2mat(tmpCell);
  outTable = array2table(tmpArray.');
  outTable.Properties.RowNames = varNames;
  writetable(outTable, exportFile, 'Delimiter', ',', 'WriteVariableNames', false, 'WriteRowNames', true);
else
  outTable = cell2table(tmpCell);
  outTable.Properties.VariableNames = varNames;
  writetable(outTable, exportFile, 'Delimiter', ',');
end



end