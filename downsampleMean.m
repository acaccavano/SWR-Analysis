function varargout = downsampleMean(inArray, varargin)
%% [outArray, outTiming] = downsampleMean(inArray, inTiming, dsTiming)
%
%  Function to downsample function to mean value of interval
%  inArray   = input column vector of data to downsample
%  inTiming  = input column vector of timing of inArray
%  dsTiming  = input column vector of desired timing of downsampled data
%  outArray  = downsampled column vector
%  outTiming = output column vector of outArray timing (Quality control - should match dsTiming, may deviate slightly for non-integer downsample factors)

% Input dsFactor
if length(varargin) == 1
  dsFactor = varargin{1};
  
% Input inTiming and dsTiming
elseif length(varargin) == 2
  inTiming   = varargin{1};
  dsTiming   = varargin{2};
  samplingIn = inTiming(2) - inTiming(1);
  samplingDS = dsTiming(2) - dsTiming(1);
  dsFactor   = samplingDS / samplingIn;
end

% Simple downsampling for integer dsFactor:
if (round(dsFactor, 3) == round(dsFactor))
  if (dsFactor >= 2)
    inArray      = inArray(1:size(inArray,1) - mod(size(inArray,1), round(dsFactor)));
    varargout{1} = mean(reshape(inArray, round(dsFactor), []))';
    if (length(varargin) == 2); varargout{2} = dsTiming; end
  else % Simple pass-through
    varargout{1} = inArray;
    if (length(varargin) == 2); varargout{2} = inTiming; end
  end
  
% Complex and computationally expensive non-integer downsample
else 
  if (dsFactor > 1)
    
    % Interpolate input array to sampling rate resulting in integer dsFactor:
    samplingInterp = samplingIn * (dsFactor / round(dsFactor)); % Sampling interval with spacing divisible by integer dsFactor
    nInterp        = ceil(1 + (inTiming(end) - inTiming(1))/samplingInterp); % Round up to capture full range
    endTimeInterp  = inTiming(1) + (nInterp - 1) * samplingInterp; % Recalculate interpolated end in case greater than inTiming(end)
    interpTiming   = linspace(inTiming(1), endTimeInterp, nInterp)';
    inArray        = interp1(inTiming, inArray, interpTiming, 'linear', 'extrap'); % Linear interpolation to samplingInt
    inTiming       = interpTiming; % Update timing array to match interpolation
    
    dsRange{length(dsTiming), 1} = []; % Cell array specifying range of inTiming to take mean of for each down-sampled point
    outArray = zeros(length(dsTiming), 1);
    outTiming = zeros(length(dsTiming), 1);
    
    % 1st range is truncated at start
    dsRange{1} = find(inTiming >= dsTiming(1) & inTiming < dsTiming(1) + 0.5*samplingDS);
    outArray(1) = mean(inArray(dsRange{1}));
    outTiming(1) = inTiming(dsRange{1}(1));
    
    for i = 2 : length(dsTiming) - 1
      dsRange{i} = find(inTiming >= dsTiming(i) - 0.5*samplingDS & inTiming < dsTiming(i) + 0.5*samplingDS);
      outArray(i) = mean(inArray(dsRange{i}));
      outTiming(i) = inTiming(dsRange{i}(floor(end/2)+1));
    end
    
    % last range is truncated at end
    dsRange{end} = find(inTiming >= dsTiming(end) - 0.5*samplingDS & inTiming <= dsTiming(end));
    outArray(end) = mean(inArray(dsRange{end}));
    outTiming(end) = inTiming(dsRange{end}(end));
    
    varargout{1} = outArray;
    varargout{2} = outTiming;
    
  else % Simple pass-through
    varargout{1} = inArray;
    varargout{2} = inTiming;
  end
end

end