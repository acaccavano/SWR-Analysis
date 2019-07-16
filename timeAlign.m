function [evStatusA1, evStartA1, evEndA1, evStatusA2, evStartA2, evEndA2, timeArray] = timeAlign(evStatus1, evStatus2, timeArray1, timeArray2)

  % Match timing of arrays
  samplingInt1 = timeArray1(2) - timeArray1(1);
  samplingInt2 = timeArray2(2) - timeArray2(1);
  samplingInt = max(samplingInt1, samplingInt2);
  minTime   = max(timeArray1(1), timeArray2(1));
  maxTime   = min(timeArray1(length(timeArray1)), timeArray2(length(timeArray2)));
  nSamples  = floor((maxTime - minTime) / samplingInt) + 1;
  timeArray = (minTime:samplingInt:maxTime)';
    
   % Downsample to match timing
  if (samplingInt == samplingInt1)
    downsampleFactor = samplingInt1 / samplingInt2;
    evStatusA1 = evStatus1(1:nSamples);
    evStatusA2 = downsampleMax(evStatus2, downsampleFactor);
    evStatusA2 = evStatusA2(1:nSamples);
  else
    downsampleFactor = samplingInt2 / samplingInt1;
    evStatusA1 = downsampleMax(evStatus1, downsampleFactor);
    evStatusA1 = evStatusA1(1:nSamples);
    evStatusA2 = evStatus2(1:nSamples);
  end

  % update start and end times in case event runs through end and has been truncated
  [evStartA1, evEndA1] = eventParse(evStatusA1);
  [evStartA2, evEndA2] = eventParse(evStatusA2);
  
end



% Upsample to match timing
%     if (samplingInt == samplingInt1)
%       stretchFactor = samplingInt2 / samplingInt1;
%       evStatus1 = evStatus1(1:nSamples);
%       evStatus2 = reshape(evStatus2(:,ones(1,stretchFactor))',[],1);
%       evStatus2 = evStatus2(1:nSamples);
%     else
%       stretchFactor = samplingInt1 / samplingInt2;
%       evStatus1 = reshape(evStatus1(:,ones(1,stretchFactor))',[],1);
%       evStatus1 = evStatus1(1:nSamples);
%       evStatus2 = evStatus2(1:nSamples);
%     end