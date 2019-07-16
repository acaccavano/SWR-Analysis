function data = rejectSWREvent(data, ev)
  % Ensures all data structures are removed, handling the optional
  % structures for gamma and cell data

  %% Remove SWR values
  % Recalculate SWR Status
  for i = data.SWR.evStart(ev):data.SWR.evEnd(ev)
    data.SWR.evStatus(i) = 0;
  end

  % Remove cell and num array values 
  data.SWR.event(ev)     = [];
  data.SWR.evStart(ev)   = [];
  data.SWR.evPeak(ev)    = [];
  data.SWR.evEnd(ev)     = [];
  data.SWR.evIndex(ev,:) = [];
  data.SWR.amp(ev)       = [];
  data.SWR.area(ev)      = [];
  data.SWR.duration(ev)  = [];
  data.SWR.power(ev)     = [];
  
  % IEI recalculation
  if (ev == 1) % Special case for first event
    % Delete this IEI
    data.SWR.IEI(ev) = [];
  elseif ev == length(data.SWR.IEI) + 1 % Special case for last event
    % Delete last IEI
    data.SWR.IEI(ev-1) = [];
  else % All in between events
    % Sum this and following IEI, and delete subsequent one
    data.SWR.IEI(ev-1) = data.SWR.IEI(ev-1) + data.SWR.IEI(ev);
    data.SWR.IEI(ev)   = [];
  end
  
  data.SWR.frequency = length(data.SWR.evStart) / ((data.LFP.timing(length(data.LFP.timing)) - data.LFP.timing(1)) / 1000);

  
  %% Remove SW values
  if isfield(data.SW,'SWR')
    data.SW.SWR.event(ev) = [];
    data.SW.SWR.area(ev)  = [];
    data.SW.SWR.power(ev) = [];
  end
  
  
  %% Remove R values
  if isfield(data.R,'SWR')
    data.R.SWR.event(ev) = [];
    data.R.SWR.power(ev) = [];
    if isfield(data.R.SWR,'FFT')
      data.R.SWR.FFT.F(ev) = [];
      data.R.SWR.FFT.F1(ev) = [];
      data.R.SWR.FFT.F2(ev) = [];
      data.R.SWR.FFT.pkFreq(ev) = [];
    end
    if isfield(data.R.SWR,'phase')
      data.R.SWR.phase.evPhase(ev) = [];
      data.R.SWR.phase.maxLoc(ev) = [];
      data.R.SWR.phase.maxVal(ev) = [];
      data.R.SWR.phase.minLoc(ev) = [];
      data.R.SWR.phase.minVal(ev) = [];
    end
  end
  
  
  %% Remove gamma values (if present)
  if isfield(data,'gamma')
    if isfield(data.gamma,'SWR')
      data.gamma.SWR.event(ev) = [];
      data.gamma.SWR.power(ev) = [];
      if isfield(data.gamma.SWR,'FFT')
        data.gamma.SWR.FFT.F(ev) = [];
        data.gamma.SWR.FFT.F1(ev) = [];
        data.gamma.SWR.FFT.F2(ev) = [];
        data.gamma.SWR.FFT.pkFreq(ev) = [];
      end
      if isfield(data.gamma.SWR,'phase')
        data.gamma.SWR.phase.evPhase(ev) = [];
        data.gamma.SWR.phase.maxLoc(ev) = [];
        data.gamma.SWR.phase.maxVal(ev) = [];
        data.gamma.SWR.phase.minLoc(ev) = [];
        data.gamma.SWR.phase.minVal(ev) = [];
      end
    end
  end
  
  
  %% Remove Cell values (if present)
  if isfield(data,'C')
    if isfield(data.C,'SWR')
      data.C.SWR.event(ev)    = [];
      data.C.SWR.evNorm(ev)   = [];
      data.C.SWR.area(ev)     = [];
      data.C.SWR.areaQ1(ev)   = [];
      data.C.SWR.areaQ2(ev)   = [];
      data.C.SWR.areaQ3(ev)   = [];
      data.C.SWR.areaQ4(ev)   = [];
      data.C.SWR.baseline(ev) = [];
    end
  end
  

end

