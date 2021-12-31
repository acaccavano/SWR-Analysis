function dataOut = gaussianFilt(dataIn, Fc1, Fc2, samplingInt, Order, Alpha)
%% dataOut = gaussianFilt(dataIn, Fc1, Fc2, samplingInt, Order, Alpha)
% 
% Gaussian FIR filter with constant and corrected phase delay

if (nargin < 6); Alpha = 2.5; end
if (nargin < 5)
  if isempty(Fc1)
    Order = 1;
  else
    if Fc1 < 1
      Order = 20;
    elseif Fc1 >= 1 && Fc1 < 2
      Order = 10;
    elseif Fc1 >= 2 && Fc1 < 5
      Order = 5;
    elseif Fc1 >= 5 && Fc1 < 10
      Order = 2;
    else
      Order = 1;
    end
  end
end

Fs = 2 * round(1000 / (2 * samplingInt));
if ~isempty(Fc2); Fc2 = min(Fs/2 - 1, Fc2); end % Filter will fail if upper limit is >= half the sampling frequency

N = round(Order * Fs);  % Order
flag = 'scale';   % Sampling Flag

% Create the window vector for the design algorithm
win = gausswin(N+1, Alpha);

% Calculate the coefficients using the FIR1 function
if ~isempty(Fc1) && ~isempty(Fc2) % Band-pass:
  b = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
elseif isempty(Fc1) % Lo-pass:
  b = fir1(N, Fc2/(Fs/2), 'low', win, flag);
elseif isempty(Fc2) % Hi-pass:
  b = fir1(N, Fc1/(Fs/2), 'high', win, flag);
end
Hd = dfilt.dffir(b);

% calculate filter delay
filtDelay = grpdelay(Hd, Fs/2, Fs);
if ~isempty(Fc1) && ~isempty(Fc2) % Band-pass:
  filtDelay = round(mean(filtDelay(ceil(Fc1):ceil(Fc2))));
elseif isempty(Fc1) % Lo-pass:
  filtDelay = round(mean(filtDelay(1:ceil(Fc2))));
elseif isempty(Fc2) % Hi-pass:
  filtDelay = round(mean(filtDelay(ceil(Fc1):ceil(Fs/2 - 1))));
end

% pad data with mirrors of beginning and end
if filtDelay < length(dataIn)
  startPad = flipud(dataIn(1 : filtDelay));
  endPad   = flipud(dataIn(length(dataIn) - filtDelay : length(dataIn)));
  dataPad  = [startPad; dataIn; endPad];
  nPads    = 0;
else
  nPads    = floor(filtDelay / length(dataIn));
  remPad   = mod(filtDelay, length(dataIn));
  oddInd   = logical(mod(1:nPads, 2));
  dataPad  = dataIn;
  
  for i = 1:nPads
    if oddInd(i)
      dataPad = [flipud(dataIn); dataPad; flipud(dataIn)];
    else
      dataPad = [dataIn; dataPad; dataIn];
    end
  end
  
  if remPad > 0
    if oddInd(i)
      startPad = dataIn(length(dataIn) - filtDelay : length(dataIn));
      endPad   = dataIn(1 : filtDelay);
    else
      startPad = flipud(dataIn(1 : filtDelay));
      endPad   = flipud(dataIn(length(dataIn) - filtDelay : length(dataIn)));
    end
    dataPad = [startPad; dataPad; endPad];
  end
end

dataOut = filter(Hd, dataPad);

if (nPads == 0)
  dataOut = dataOut(2 * (filtDelay + 1) : length(dataOut));
else
  dataOut = dataOut(2 * filtDelay + 1 : length(dataOut)); 
end

end