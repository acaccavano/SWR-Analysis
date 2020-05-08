function [evStatusA1, evStartA1, evEndA1, evStatusA2, evStartA2, evEndA2, timeArray] = timeAlign(evStatus1, evStatus2, timeArray1, timeArray2, alignEndOption)
%% [evStatusA1, evStartA1, evEndA1, evStatusA2, evStartA2, evEndA2, timeArray] = timeAlign(evStatus1, evStatus2, timeArray1, timeArray2)
%
%  Function to align arrays with two different timing, downsampling file
%  with higher sampling rate and trimming as necessary

if (nargin < 5) alignEndOption = 0; end

% If aligning from start, assume nonzero start time is meaningful (i.e. skipping beginning of dFoF files)
% However, if aligning from end, reset to common end-time:
if alignEndOption
  maxTime    = max(timeArray1(end), timeArray2(end));
  timeArray1 = timeArray1 + (maxTime - timeArray1(end));
  timeArray2 = timeArray2 + (maxTime - timeArray2(end));
end

% Deterine largest samplingInt:
samplingInt1 = timeArray1(2) - timeArray1(1);
samplingInt2 = timeArray2(2) - timeArray2(1);
samplingInt  = max(samplingInt1, samplingInt2);

% Determine most restrictive bounds for common time array:
minTime      = max(timeArray1(1), timeArray2(1));
maxTime      = min(timeArray1(end), timeArray2(end));
nSamples     = floor((maxTime - minTime) / samplingInt) + 1;
timeArray    = (minTime:samplingInt:maxTime)';

% Match Timing:
if (samplingInt == samplingInt1) % evStatus2 to be downsampled
 
  if ~alignEndOption % Align from start of files and trim end:
    evStatusA1 = evStatus1(1:nSamples);
  else % Align from end of files and trim beginning:
    evStatusA1 = evStatus1(end-nSamples+1:end);
  end
  
  % Downsample evStatus2
  [evStatusA2, ~] = downsampleMax(evStatus2, timeArray2, timeArray);

elseif (samplingInt == samplingInt2) % evStatus1 to be downsampled
  
  % Downsample evStatus1
  [evStatusA1, ~] = downsampleMax(evStatus1, timeArray1, timeArray);
  
  if alignEndOption % Align from end of files:
    evStatusA2 = evStatus2(1:nSamples);
  else % Align from start of files:
    evStatusA2 = evStatus2(1:nSamples);
  end
end

% update start and end times in case event runs through end and has been truncated
[evStartA1, evEndA1] = eventParse(evStatusA1);
[evStartA2, evEndA2] = eventParse(evStatusA2);

end