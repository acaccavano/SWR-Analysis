function S = calcTotPhase(S, LFP, param)
%% S = calcTotPhase(S, LFP)
%
%  Function to calculate the phase for specific frequency range of signal S using a
%  piece-wise linear fit to maxima and minima of filtered signal, with some
%  error handling to ensure abberently high or low results are discarded.
%  Requires LFP structure for timing

if ~isfield(S,'phase') S.phase = struct; end
if ~isfield(param,'minPkDiffFrac') param.minPkDiffFrac = 1.5;  end % Fraction of max frequency (lim2) to set min peak difference for exclusion of events. 1 = max frequency of band, is pretty conservative, default = 1.5, so two successive oscillations can be slightly faster than allowed
if ~isfield(param,'minPkPromFrac') param.minPkPromFrac = 0.0;  end % Fraction of SD of signal to set Min peak prominence (zero to ignore as criteria)
if ~isfield(param,'samplingInt')   param.samplingInt = LFP.timing(2) - LFP.timing(1); end

% (Re)Initialize data arrays:
S.phase.maxVal  = [];
S.phase.maxLoc  = [];
S.phase.minVal  = [];
S.phase.minLoc  = [];
S.phase.tPhase  = [];
S.phase.nCycle  = [];
S.phase.phFreq  = [];

%% Compute oscillation phase:
minPkDiff = 1000/(param.minPkDiffFrac * S.lim2 * param.samplingInt);
minPkProm = param.minPkPromFrac * std(S.tSeries);

[S.phase.maxVal, S.phase.maxLoc] = findpeaks(S.tSeries, 'MinPeakDistance', minPkDiff, 'MinPeakProminence', minPkProm);
S.phase.minVal = NaN * ones(length(S.phase.maxVal) - 1, 1);
S.phase.minLoc = NaN * ones(length(S.phase.maxVal) - 1, 1);

% Piecewise Linear Interpolation
S.phase.tPhase = zeros(length(S.tSeries), 1);

formatSpec = fprintfInline('interpolating piecewise linear phase function for peak =  ', length(S.phase.maxLoc)-1);
for i = 1 : length(S.phase.maxLoc) - 1
  
  [S.phase.minVal(i), S.phase.minLoc(i)] = min(S.tSeries(S.phase.maxLoc(i) : S.phase.maxLoc(i+1)));
  S.phase.minLoc(i) = S.phase.minLoc(i) + S.phase.maxLoc(i);
  
  tSubArray1 = S.phase.maxLoc(i) : S.phase.minLoc(i);
  tSubArray2 = S.phase.minLoc(i) : S.phase.maxLoc(i+1);
  S.phase.tPhase(tSubArray1) = linspace(0, pi, length(tSubArray1));
  S.phase.tPhase(tSubArray2) = linspace(pi, 2*pi, length(tSubArray2));
  fprintf(formatSpec,i);
end
fprintf('...done\n');

% Calculate number of cycles and frequency
S.phase.nCycle = length(S.phase.maxLoc) - 1;
S.phase.phFreq = S.phase.nCycle / ((LFP.timing(S.phase.maxLoc(end)) - LFP.timing(S.phase.maxLoc(1)))/1000);

% Error handling: limit phase frequency in case it exceeds frequency bounds
% if S.phase.phFreq < S.lim1
%   S.phase.phFreq = S.lim1;
%   S.phase.nCycle = S.phase.phFreq * ((LFP.timing(end) - LFP.timing(1))/1000);
% elseif S.phase.phFreq > S.lim2
%   S.phase.phFreq = S.lim2;
%   S.phase.nCycle = S.phase.phFreq * ((LFP.timing(end) - LFP.timing(1))/1000);
% end

S.phase = orderStruct(S.phase);

end