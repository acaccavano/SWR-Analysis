function [evStatus, evStart, evPeak, evEnd] = peakFindUnique(dataInp, timeArray, peakThresh, baseThresh, peakPolarity)
%% [evStatus, evStart, evPeak, evEnd] = peakFindUnique(dataInp, timeArray, peakThresh, baseThresh, peakPolarity)
% 
%  Function to find unique peaks given two thresholds, one for peak
%  detection and a second optionally lower one for start and end times.
%  peakPolarity = +1 for maxima, -1 for minima

% Initialize array for raster plots:
evStatus = zeros(size(timeArray,1),1);

% Find all peaks > peakThresh
[peakVal, peakLoc] = findpeaks(peakPolarity * dataInp,'MinPeakHeight', peakThresh);

if isempty(peakLoc)
  evStart = NaN;
  evPeak  = NaN;
  evEnd   = NaN;
else
  % Flip data for finding start of peak
  dataRev = flipud(dataInp);
  peakRev = size(timeArray,1) - peakLoc;
  eventj = 1; % index for final non-duplicated array
  
  % Find start and end of each event based on baseThresh
  for eventi = 1:size(peakLoc,1)
    startWidth = find(peakPolarity * dataRev(peakRev(eventi) : size(timeArray,1)) <= baseThresh, 1) + 1;
    endWidth   = find(peakPolarity * dataInp(peakLoc(eventi) : size(timeArray,1)) <= baseThresh, 1) - 1;
    
    % In case event runs through beginning or end of dataInp
    if isempty(startWidth)
      startTemp = 1;
    else
      startTemp = max(peakLoc(eventi) - startWidth, 1);
    end
    if isempty(endWidth)
      endTemp = size(timeArray,1);
    else
      endTemp = min(peakLoc(eventi) + endWidth, size(timeArray,1));
    end
    
    % Check if duplicate event
    if (eventi~=1) && (startTemp==evStart(eventj)) && (endTemp==evEnd(eventj))
      % Update peak if greater
      if peakVal(eventi) > peakPolarity * dataInp(evPeak(eventj))
        evPeak(eventj) = peakLoc(eventi);
      end
      % Non- duplicate event - update everything
    else
      % First iteration
      if (eventi~=1)
        eventj = eventj + 1;
      end
      evStart(eventj) = startTemp;
      evPeak(eventj)  = peakLoc(eventi);
      evEnd(eventj)   = endTemp;
      % Fill arrays for raster plot
      for i = evStart(eventj):evEnd(eventj)
        evStatus(i) = 1;
      end
    end
  end
  evStart = evStart';
  evPeak = evPeak';
  evEnd = evEnd';
end
end