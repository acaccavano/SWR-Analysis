function [evStatus, evStart, evEnd, evIndex] = eventOverlap(evStatus1, evStart1, evEnd1, evStatus2, evStart2, evEnd2, timeArray, overlapType)
  % [evStatus, evStart, evEnd] = eventOverlap(evStatus1, evStart1, evEnd1, evStatus2, evStart2, evEnd2, timeArray, overlapType)
  % 
  % Function takes two different event statuses and outputs the status of overlapping events
  % overlapType = 0 (or(evStatus1,evStatus2), 1 (evStatus1), 2 (evStatus2)
  
  if (nargin < 8) overlapType = 0; end

  evConj = evStatus1 .* evStatus2;
  evDisj = or(evStatus1, evStatus2);
  
  [evConjStart, evConjEnd] = eventParse(evConj);
  [evDisjStart, evDisjEnd] = eventParse(evDisj);
  
  evStatus = zeros(size(timeArray,1),1);
  evStart  = zeros(size(evConjStart,1),1);
  evEnd    = zeros(size(evConjStart,1),1);

  % Select correct array for final timing based on overlapType
  if (overlapType == 0)
    evTimingStart = evDisjStart;
    evTimingEnd   = evDisjEnd;
  elseif (overlapType == 1)
    evTimingStart = evStart1;
    evTimingEnd   = evEnd1;
  elseif (overlapType == 2)
    evTimingStart = evStart2;
    evTimingEnd   = evEnd2;
  end
  
  eventk = 1;
  for eventi = 1:size(evConjStart,1)
    eventj = find(evTimingStart > evConjStart(eventi), 1);

    if isempty(eventj) % Special case for last event
      eventj = size(evTimingStart,1);
    else % Normal case
      eventj = eventj - 1;
    end

    evStart(eventk,1) = evTimingStart(eventj);
    evEnd(eventk,1) = evTimingEnd(eventj);
    
    % Check if duplicate event
    if (eventk>1) && (evStart(eventk) == evStart(eventk-1)) && (evEnd(eventk) == evEnd(eventk-1))
      evStart = evStart(1:eventk-1,1);
      evEnd = evEnd(1:eventk-1,1);
      for i = evStart(eventk-1):evEnd(eventk-1)
        evStatus(i) = 1;
      end
    else
      for i = evStart(eventk):evEnd(eventk)
        evStatus(i) = 1;
      end
      eventk = eventk + 1;
    end
  end
  
  % Calculate Index:
  evIndex = zeros(length(evStart),2);
  for ev = 1:length(evStart)
    evIndex(ev,1) = find(evStart1 <= evEnd(ev) & evEnd1 >= evStart(ev), 1);
    evIndex(ev,2) = find(evStart2 <= evEnd(ev) & evEnd2 >= evStart(ev), 1);
  end
  
end