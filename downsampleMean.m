function outArray  = downsampleMean(inArray, dsFactor)
%% outArray  = downsampleMean(inArray, dsFactor)
%
%  Function to downsample to mean value of interval
%  inArray  = input column vector
%  outArray = downsampled column vector
%  dsFactor = factor with which to downsample

if dsFactor >= 2
  dsFactor = round(dsFactor);
  inArray  = inArray(1:size(inArray,1) - mod(size(inArray,1),dsFactor));
  outArray = mean(reshape(inArray,dsFactor,[]))';
else
  outArray = inArray;
end

end

      