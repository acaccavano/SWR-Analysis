function S = calcEvPhase(S, SWR, param, lim1, lim2)
%% S = calcEvPhase(S, SWR, param, lim1, lim2)
%
%  Function to calculate the phase for specific a frequency range of SWR events using a
%  piece-wise linear fit to maxima and minima of filtered signal, with some
%  error handling to ensure abberently high or low results are discarded. 

if ~isfield(S,'phase'); S.phase = struct; end
if ~isfield(param,'minPkDiffFrac'); param.minPkDiffFrac = 1.5;  end % Fraction of max frequency (lim2) to set min peak difference for exclusion of events. 1 = max frequency of band, is pretty conservative, default = 1.5, so two successive oscillations can be slightly faster than allowed
if ~isfield(param,'minPkPromFrac'); param.minPkPromFrac = 0.0;  end % Fraction of SD of signal to set Min peak prominence (zero to ignore as criteria)
if ~isfield(param,'samplingInt');   param.samplingInt = SWR.evTiming(2) - SWR.evTiming(1); end

% (Re)Initialize data arrays:
S.phase.maxVal   = [];
S.phase.maxLoc   = [];
S.phase.minVal   = [];
S.phase.minLoc   = [];
S.phase.evPhase  = [];

S.phase.maxVal{length(S.event)}  = [];
S.phase.maxLoc{length(S.event)}  = [];
S.phase.minVal{length(S.event)}  = [];
S.phase.minLoc{length(S.event)}  = [];
S.phase.evPhase{length(S.event)} = [];

S.phase.maxVal  = S.phase.maxVal';
S.phase.maxLoc  = S.phase.maxLoc';
S.phase.minVal  = S.phase.minVal';
S.phase.minLoc  = S.phase.minLoc';
S.phase.evPhase = S.phase.evPhase';

S.phase.nCycle = NaN * zeros(length(S.event), 1);
S.phase.phFreq = NaN * zeros(length(S.event), 1);

%% Compute oscillation phase for valid events:
nSamples  = max(cellfun(@length,S.event));
minPkDiff = 1000/(param.minPkDiffFrac * lim2 * param.samplingInt);

for ev = 1:length(S.event)
  if ~isempty(S.event{ev}) && (length(S.event{ev}) == nSamples)
    
    minPkProm = param.minPkPromFrac * std(S.event{ev});
    
    [S.phase.maxVal{ev}, S.phase.maxLoc{ev}] = findpeaks(S.event{ev}, 'MinPeakDistance', minPkDiff, 'MinPeakProminence', minPkProm);
    S.phase.minVal{ev} = NaN * ones(length(S.phase.maxVal) - 1, 1);
    S.phase.minLoc{ev} = NaN * ones(length(S.phase.maxVal) - 1, 1);

    S.phase.evPhase{ev} = zeros(length(S.event{ev}), 1);

    for i = 1 : length(S.phase.maxLoc{ev}) - 1
      
      [S.phase.minVal{ev}(i), S.phase.minLoc{ev}(i)] = min(S.event{ev}(S.phase.maxLoc{ev}(i) : S.phase.maxLoc{ev}(i+1)));
      S.phase.minLoc{ev}(i) = S.phase.maxLoc{ev}(i) + S.phase.minLoc{ev}(i);
      
      tSubArray1 = S.phase.maxLoc{ev}(i) : S.phase.minLoc{ev}(i);
      tSubArray2 = S.phase.minLoc{ev}(i) : S.phase.maxLoc{ev}(i+1);
      S.phase.evPhase{ev}(tSubArray1) = linspace(0, pi, length(tSubArray1));
      S.phase.evPhase{ev}(tSubArray2) = linspace(pi, 2*pi, length(tSubArray2));
    end
    
    %% Calculate number of cycles and frequency
    
    if isempty(S.phase.maxLoc{ev}(1)); S.phase.maxLoc{ev}(1) = 1; end
    if isempty(S.phase.maxLoc{ev}(end)); S.phase.maxLoc{ev}(end) = length(SWR.evTiming); end
    
    midPoint = floor(length(SWR.evTiming)/2) + 1;
    evStartStand = max(midPoint - (SWR.evPeak(ev) - SWR.evStart(ev)), S.phase.maxLoc{ev}(1));
    evEndStand   = min(midPoint + (SWR.evEnd(ev)  - SWR.evPeak(ev)), S.phase.maxLoc{ev}(end));
    
    contPhase = S.phase.evPhase{ev}(evStartStand:evEndStand);
    nCycle = 0;
    
    for i = 1:length(contPhase)
      if (contPhase(i) == 0); nCycle = nCycle + 1; end
      contPhase(i) = contPhase(i) + pi*2*nCycle;
    end
    
    S.phase.nCycle(ev) = (contPhase(end) - contPhase(1)) / (2*pi);
    S.phase.phFreq(ev) = S.phase.nCycle(ev) / (SWR.duration(ev)/1000);
    
    % Error handling: limit phase frequency in case it exceeds frequency bounds
    if S.phase.phFreq(ev) < lim1
      S.phase.phFreq(ev) = lim1;
      S.phase.nCycle(ev) = S.phase.phFreq(ev) * (SWR.duration(ev)/1000);
    elseif S.phase.phFreq(ev) > lim2
      S.phase.phFreq(ev) = lim2;
      S.phase.nCycle(ev) = S.phase.phFreq(ev) * (SWR.duration(ev)/1000);
    end
    
  end
end

S.phase = orderStruct(S.phase);

end