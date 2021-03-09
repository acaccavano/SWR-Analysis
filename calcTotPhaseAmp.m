function phaseAmp = calcTotPhaseAmp(phase, peak)
% Function to determine amplitude of oscillation surrounding spike based on phase analysis

precMin  = find(phase.minLoc < peak, 1, 'last');
procMin  = find(phase.minLoc >= peak, 1, 'first');
precMax  = find(phase.maxLoc < peak, 1, 'last');
procMax  = find(phase.maxLoc >= peak, 1, 'first');

if ~isempty(precMin) && ~isempty(precMax) && ~isempty(procMin) && ~isempty(procMax)
  if precMin > precMax
    phaseAmp = phase.maxVal(procMax) - phase.minVal(precMin);
  else
    phaseAmp = phase.maxVal(precMax) - phase.minVal(procMin);
  end
else
  phaseAmp = 0;
end