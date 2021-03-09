function phaseAmp = calcEvPhaseAmp(phase, ev, peak)
% Function to determine amplitude of oscillation surrounding spike based on phase analysis

precMin  = find(phase.minLoc{ev} < peak, 1, 'last');
procMin  = find(phase.minLoc{ev} >= peak, 1, 'first');
precMax  = find(phase.maxLoc{ev} < peak, 1, 'last');
procMax  = find(phase.maxLoc{ev} >= peak, 1, 'first');

if ~isempty(precMin) && ~isempty(precMax) && ~isempty(procMin) && ~isempty(procMax)
  if precMin > precMax
    phaseAmp = phase.maxVal{ev}(procMax) - phase.minVal{ev}(precMin);
  else
    phaseAmp = phase.maxVal{ev}(precMax) - phase.minVal{ev}(procMin);
  end
else
  phaseAmp = 0;
end