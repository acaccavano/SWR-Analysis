function [evStart, evEnd] = eventParse(evStatus)

  evStart = [];
  evEnd = [];

  i = 1;
  j = 1;
  
  while i < size(evStatus,1)
    % Search for start of event
    evSearch = find(evStatus(i : size(evStatus,1)), 1);
    
    % No more events
    if isempty(evSearch)
      i = size(evStatus,1);
      
    % At least one more event
    else
      evStart(j) = i + evSearch - 1;
      evSearch = find(~evStatus(evStart(j) : size(evStatus,1)), 1);
      
      % If no end found - must be last event
      if isempty(evSearch)
        evEnd(j) = size(evStatus,1);
      else
        evEnd(j) = evStart(j) + evSearch - 2;
      end
      
      % Increment indeces
      i = evEnd(j) + 1;        
      j = j + 1;
    end
  end
  
  % Transpose to column vectors
  evStart = evStart';
  evEnd = evEnd';
end