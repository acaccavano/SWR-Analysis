function S = calcEvPhase(S, tArray, maxFreq)
  
  if ~isfield(S,'phase') S.phase = struct; end

  splineOption = false;
  
  % (Re)Initialize data arrays:
  S.phase.maxVal   = [];
  S.phase.maxLoc   = [];
  S.phase.minVal   = [];
  S.phase.minLoc   = [];
  S.phase.evPhase  = [];
  S.phase.evPhaseX = [];
  S.phase.evPhaseY = [];
  
  S.phase.maxVal{length(S.event)}   = [];
  S.phase.maxLoc{length(S.event)}   = [];
  S.phase.minVal{length(S.event)}   = [];
  S.phase.minLoc{length(S.event)}   = [];
  S.phase.evPhase{length(S.event)}  = [];
  S.phase.evPhase{length(S.event)}  = [];
  S.phase.evPhaseX{length(S.event)} = [];
  S.phase.evPhaseY{length(S.event)} = [];
  
  S.phase.maxVal  = S.phase.maxVal';
  S.phase.maxLoc  = S.phase.maxLoc';
  S.phase.minVal  = S.phase.minVal';
  S.phase.minLoc  = S.phase.minLoc';
  S.phase.evPhase = S.phase.evPhase';

  % Compute ripple phase for valid events:
  nSamples = max(cellfun(@length,S.event));
  for ev = 1:length(S.event)
    if ~isempty(S.event{ev}) && (length(S.event{ev}) == nSamples)

      [S.phase.maxVal{ev}, S.phase.maxLoc{ev}] = findpeaks(S.event{ev}, tArray, 'MinPeakDistance', 1/maxFreq);
      [S.phase.minVal{ev}, S.phase.minLoc{ev}] = findpeaks(-S.event{ev}, tArray, 'MinPeakDistance', 1/maxFreq);
      S.phase.minVal{ev} = -S.phase.minVal{ev};
      
      if ~splineOption % Use linear interpolation
        S.phase.evPhase{ev} = zeros(length(S.event{ev}), 1);
        
        if S.phase.minLoc{ev}(1) < S.phase.maxLoc{ev}(1)
          S.phase.minLoc{ev}(1) = [];
          S.phase.minVal{ev}(1) = [];
        end
        
        for i = 1 : length(S.phase.maxLoc{ev})-1
          tSubArray1 = find(tArray>=S.phase.maxLoc{ev}(i) & tArray<=S.phase.minLoc{ev}(i));
          tSubArray2 = find(tArray>=S.phase.minLoc{ev}(i) & tArray<=S.phase.maxLoc{ev}(i+1));
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
        
        S.phase.evPhase{ev} = mod(spline(sort([phase.maxLoc{ev};phase.minLoc{ev}]), phase1, tArray),2*pi);
      end
    end
  end

  S.phase = orderStruct(S.phase);
  
end