function S = calcTotPhase(S, LFP)
%% S = calcTotPhase(S, LFP)
%
%  Function to calculate the phase for specific frequency range of signal S using a
%  piece-wise linear fit to maxima and minima of filtered signal, with some
%  error handling to ensure abberently high or low results are discarded.
%  Requires LFP structure for timing

if ~isfield(S,'phase') S.phase = struct; end

splineOption = false;

% (Re)Initialize data arrays:
S.phase.maxVal  = [];
S.phase.maxLoc  = [];
S.phase.minVal  = [];
S.phase.minLoc  = [];
S.phase.tPhase  = [];
S.phase.nCycle  = [];
S.phase.phFreq  = [];

%% Compute oscillation phase:
[S.phase.maxVal, S.phase.maxLoc] = findpeaks(S.tSeries, LFP.timing, 'MinPeakDistance', 1/S.lim2);
[S.phase.minVal, S.phase.minLoc] = findpeaks(-S.tSeries, LFP.timing, 'MinPeakDistance', 1/S.lim2);
S.phase.minVal = -S.phase.minVal;

if ~splineOption % Use linear interpolation
  S.phase.tPhase = zeros(length(S.tSeries), 1);
  
  if S.phase.minLoc(1) < S.phase.maxLoc(1)
    S.phase.minLoc(1) = [];
    S.phase.minVal(1) = [];
  end
  
  for i = 1 : length(S.phase.maxLoc)-1
    tSubArray1 = find(LFP.timing>=S.phase.maxLoc(i) & LFP.timing<=S.phase.minLoc(i));
    tSubArray2 = find(LFP.timing>=S.phase.minLoc(i) & LFP.timing<=S.phase.maxLoc(i+1));
    S.phase.tPhase(tSubArray1) = linspace(0, pi, length(tSubArray1));
    S.phase.tPhase(tSubArray2) = linspace(pi, 2*pi, length(tSubArray2));
  end
  
else % Use cubic spline interpolation - CAUTION: only works for very well-defined oscillations
  if S.phase.maxLoc(1) < S.phase.minLoc(1)
    phaseMax = pi*2*(0:length(S.phase.maxLoc)-1);
  else
    phaseMax = pi*2*(1:length(S.phase.maxLoc));
  end
  phaseMin = pi*(2*(1:length(S.phase.minLoc))-1);
  phase1 = sort([phaseMax, phaseMin])';
  
  S.phase.tPhase = mod(spline(sort([S.phase.maxLoc;S.phase.minLoc]), phase1, LFP.timing),2*pi);
end

%% Calculate number of cycles and frequency
firstPeak = find(S.phase.maxLoc(1) == LFP.timing);
lastPeak = find(S.phase.maxLoc(end) == LFP.timing);

if isempty(firstPeak) firstPeak = 1; end
if isempty(lastPeak) lastPeak = length(LFP.timing); end

contPhase = S.phase.tPhase(firstPeak:lastPeak);
nCycle = 0;

for i = 1:length(contPhase)
  if (contPhase(i) == 0) nCycle = nCycle + 1; end
  contPhase(i) = contPhase(i) + pi*2*nCycle;
end

S.phase.nCycle = (contPhase(end) - contPhase(1)) / (2*pi);
S.phase.phFreq = S.phase.nCycle / ((LFP.timing(end) - LFP.timing(1))/1000);

% Error handling: limit phase frequency in case it exceeds frequency bounds
if S.phase.phFreq < S.lim1
  S.phase.phFreq = S.lim1;
  S.phase.nCycle = S.phase.phFreq * ((LFP.timing(end) - LFP.timing(1))/1000);
elseif S.phase.phFreq > S.lim2
  S.phase.phFreq = S.lim2;
  S.phase.nCycle = S.phase.phFreq * ((LFP.timing(end) - LFP.timing(1))/1000);
end

S.phase = orderStruct(S.phase);

end