function [evStatus, evStart, evPeak, evEnd] = peakFindUnique(dataIn, timeArray, peakThresh, baseThresh, peakPolarity, minEvDiff)
%% [evStatus, evStart, evPeak, evEnd] = peakFindUnique(dataIn, timeArray, peakThresh, baseThresh, peakPolarity)
% 
%  Function to find unique peaks given two thresholds, one for peak
%  detection and a second optionally lower one for start and end times.
%  peakPolarity = +1 for maxima, -1 for minima

if (nargin < 6); minEvDiff   = 0;   end

% Initialize array for raster plots:
evStatus = zeros(size(timeArray,1),1);
samplingInt = timeArray(2) - timeArray(1);

% Find all peaks > peakThresh
[peakVal, peakLoc] = findpeaks(peakPolarity * dataIn,'MinPeakHeight', peakThresh);

if isempty(peakLoc)
  evStart = NaN;
  evPeak  = NaN;
  evEnd   = NaN;
else
  % Flip data for finding start of peak
  dataRev = flipud(dataIn);
  peakRev = size(timeArray,1) - peakLoc;
  ev2 = 1; % index for final non-duplicated array
  
  % Find start and end of each event based on baseThresh
  for ev1 = 1:size(peakLoc,1)
    startWidth = find(peakPolarity * dataRev(peakRev(ev1) : size(timeArray,1)) <= baseThresh, 1) + 1;
    endWidth   = find(peakPolarity * dataIn(peakLoc(ev1) : size(timeArray,1)) <= baseThresh, 1) - 1;
    
    % Find start and end locations, with logic for when event runs through beginning or end of dataIn
    if isempty(startWidth)
      startTemp = 1;
    else
      startTemp = max(peakLoc(ev1) - startWidth, 1);
    end
    if isempty(endWidth)
      endTemp = size(timeArray,1);
    else
      endTemp = min(peakLoc(ev1) + endWidth, size(timeArray,1));
    end
    
    % Check if duplicate event
    if (ev1~=1) && (startTemp==evStart(ev2)) && (endTemp==evEnd(ev2))
      % Update peak if greater
      if peakVal(ev1) > peakPolarity * peakVal(ev2)
        evPeak(ev2) = peakLoc(ev1);
      end
      
    % Check if min difference between events is violated, if so combine events
    elseif (ev1~=1) && (startTemp - evEnd(ev2)) * samplingInt < minEvDiff
      % Update peak if greater
      if peakVal(ev1) > peakPolarity * peakVal(ev2)
        evPeak(ev2) = peakLoc(ev1);
      end
      % Set end of last peak to end of current:
      evEnd(ev2)   = endTemp;
      % Fill arrays for raster plot
      for i = evStart(ev2):evEnd(ev2)
        evStatus(i) = 1;
      end
      
    % Non- duplicate event - update everything
    else
      % First iteration
      if (ev1~=1)
        ev2 = ev2 + 1;
      end
      evStart(ev2) = startTemp;
      evPeak(ev2)  = peakLoc(ev1);
      evEnd(ev2)   = endTemp;
      % Fill arrays for raster plot
      for i = evStart(ev2):evEnd(ev2)
        evStatus(i) = 1;
      end
    end
  end
  evStart = evStart';
  evPeak = evPeak';
  evEnd = evEnd';
end
end