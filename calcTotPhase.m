function S = calcTotPhase(S, SWR, minFreq, maxFreq)
%% S = calcEvPhase(S, SWR, minFreq, maxFreq)
%
%  Function to calculate the phase for specific frequency range of signal using a
%  piece-wise linear fit to maxima and minima of filtered signal, with some
%  error handling to ensure abberently high or low results are discarded. 

if ~isfield(S,'phase') S.phase = struct; end

splineOption = false;

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
nSamples = max(cellfun(@length,S.event));
for ev = 1:length(S.event)
  if ~isempty(S.event{ev}) && (length(S.event{ev}) == nSamples)
    
    [S.phase.maxVal{ev}, S.phase.maxLoc{ev}] = findpeaks(S.event{ev}, SWR.evTiming, 'MinPeakDistance', 1/maxFreq);
    [S.phase.minVal{ev}, S.phase.minLoc{ev}] = findpeaks(-S.event{ev}, SWR.evTiming, 'MinPeakDistance', 1/maxFreq);
    S.phase.minVal{ev} = -S.phase.minVal{ev};
    
    if ~splineOption % Use linear interpolation
      S.phase.evPhase{ev} = zeros(length(S.event{ev}), 1);
      
      if S.phase.minLoc{ev}(1) < S.phase.maxLoc{ev}(1)
        S.phase.minLoc{ev}(1) = [];
        S.phase.minVal{ev}(1) = [];
      end
      
      for i = 1 : length(S.phase.maxLoc{ev})-1
        tSubArray1 = find(SWR.evTiming>=S.phase.maxLoc{ev}(i) & SWR.evTiming<=S.phase.minLoc{ev}(i));
        tSubArray2 = find(SWR.evTiming>=S.phase.minLoc{ev}(i) & SWR.evTiming<=S.phase.maxLoc{ev}(i+1));
        S.phase.evPhase{ev}(tSubArray1) = linspace(0, pi, length(tSubArray1));
        S.phase.evPhase{ev}(tSubArray2) = linspace(pi, 2*pi, length(tSubArray2));
      end
      
    else % Use cubic spline interpolation - CAUTION: only works for very well-defined oscillations
      if S.phase.maxLoc{ev}(1) < S.phase.minLoc{ev}(1)
        phaseMax = pi*2*(0:length(S.phase.maxLoc{ev})-1);
      else
        phaseMax = pi*2*(1:length(S.phase.maxLoc{ev}));
      end
      phaseMin = pi*(2*(1:length(S.phase.minLoc{ev}))-1);
      phase1 = sort([phaseMax, phaseMin])';
      
      S.phase.evPhase{ev} = mod(spline(sort([phase.maxLoc{ev};phase.minLoc{ev}]), phase1, SWR.evTiming),2*pi);
    end
    
    %% Calculate number of cycles and frequency
    firstPeak = find(S.phase.maxLoc{ev}(1) == SWR.evTiming);
    lastPeak = find(S.phase.maxLoc{ev}(end) == SWR.evTiming);
    
    if isempty(firstPeak) firstPeak = 1; end
    if isempty(lastPeak) lastPeak = length(SWR.evTiming); end
    
    midPoint = floor(length(SWR.evTiming)/2) + 1;
    evStartStand = max(midPoint - (SWR.evPeak(ev) - SWR.evStart(ev)), firstPeak);
    evEndStand   = min(midPoint + (SWR.evEnd(ev)  - SWR.evPeak(ev)), lastPeak);
    
    contPhase = S.phase.evPhase{ev}(evStartStand:evEndStand);
    nCycle = 0;
    
    for i = 1:length(contPhase)
      if (contPhase(i) == 0) nCycle = nCycle + 1; end
      contPhase(i) = contPhase(i) + pi*2*nCycle;
    end
    
    S.phase.nCycle(ev) = (contPhase(end) - contPhase(1)) / (2*pi);
    S.phase.phFreq(ev) = S.phase.nCycle(ev) / (SWR.duration(ev)/1000);
    
    % Error handling: limit phase frequency in case it exceeds frequency bounds
    if S.phase.phFreq(ev) < minFreq
      S.phase.phFreq(ev) = minFreq;
      S.phase.nCycle(ev) = S.phase.phFreq(ev) * (SWR.duration(ev)/1000);
    elseif S.phase.phFreq(ev) > maxFreq
      S.phase.phFreq(ev) = maxFreq;
      S.phase.nCycle(ev) = S.phase.phFreq(ev) * (SWR.duration(ev)/1000);
    end
    
  end
end

S.phase = orderStruct(S.phase);

end