function xFreq = calcPACMI(modPhase, pacAmp, nSample, Fs, nShuffle)
% Computes normalized phase-amplitude coupled modulation index
% Adapted from Canolty et al Science 2006
% See also Onslow et al 2011

z = pacAmp .* exp(1i * modPhase); % Create composite signal
pacRaw = mean(z);  % Compute the mean length of composite signal
xFreq.pacMI_Raw = abs(pacRaw);
pacShuffle = zeros(nShuffle, 1);

% Compute shuffled values:
for i = 1 : nShuffle
  skip = shuffleInd(nSample, Fs, nShuffle);
  pacAmpShuffle = [pacAmp(skip(i):end) pacAmp(1:skip(s)-1)];
  pacShuffle(i) = abs(mean(pacAmpShuffle .* exp(1i * modPhase)));
end

% Fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox 
[xFreq.shuffleMean, xFreq.shuffleSTD] = normfit(pacShuffle); 

% Normalize length using shuffled data (z-score) 
xFreq.pacMI_Length = (xFreq.pacMI_Raw - xFreq.shuffleMean) / xFreq.shuffleSTD; 
xFreq.pacMI_Phase = angle(pacRaw);
xFreq.pacMI_Norm = pacMI_Length * exp(1i * pacMI_Phase);