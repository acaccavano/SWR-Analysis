function S = calcPhaseStats(S)
% Calculate circular statistics

phase = S.phase(~isnan(S.phase));

if ~isempty(phase)
  S.phaseAve = circ_mean(phase);
  if (S.phaseAve < 0) 
    S.phaseAve = S.phaseAve + 2*pi; 
  end
  S.phaseR   = circ_r(phase);
  if length(phase) > 1
    [S.phaseP, S.phaseZ] = circ_rtest(phase);
  end
end