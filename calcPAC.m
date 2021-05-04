function [miPAC, miPACWin, morlFreqP, morlFreqA, modPhase, pacAmp] = calcPAC(data, param)

% Below is simplification of code written by Author: Angela Onslow, May
% 2010, to handle scalar instead of vector analysis.
% Test 

param.morlWidth  = 7;
param.winLength  = 0.5; % [s]
param.winOverlap = 0.2; % [s]

% Shorten signals to get integer number of time windows
nSample  = length(data.LFP.tSeries);
nWin     = ceil(param.winLength * param.Fs);
nOverlap = ceil(param.winOverlap * param.Fs);
tSeries  = data.LFP.tSeries(1 : nSample - mod(nSample, nWin));

% Update nSample
nSample = length(tSeries);
idx     = bsxfun(@plus, (1:nWin)', 1+(0:(fix((nSample - nOverlap)/(nWin - nOverlap)) - 1))*(nWin - nOverlap)) - 1;
%   xbins   = size(idx, 2);
%   alpha = alpha/(xbins*ybins); % Uncomment to use Bonferonni Correction

% Calculate modulating phase (lower frequency) via Morlet wavelet:
morlFreqP = data.theta.lim1 + floor((data.theta.lim2 - data.theta.lim1)/2);
modPhase  = morletPhase(morlFreqP, tSeries, param.Fs, param.morlWidth);

% Calculate amplitude (higher frequency) via Morlet wavelet:
morlFreqA = data.gamma.lim1 + floor((data.gamma.lim2 - data.gamma.lim1)/2);
pacAmp    = morletAmp(morlFreqA, tSeries, param.Fs, param.morlWidth);

% Calculate total PAC measure:
z = pacAmp .* exp(1i * modPhase); % Create composite signal
miRaw = mean(z);  % Compute the mean length of composite signal
miPAC = abs(miRaw);

% Calculate windowed time-series PAC:
miPACWin = zeros(size(idx,2), 1);
for i = 1:size(idx, 2)
  
  z = pacAmp(idx(:,i)) .* exp(1i * modPhase(idx(:,1))); % Create composite signal
  miRaw = mean(z);  % Compute the mean length of composite signal
  miPACWin(i) = abs(miRaw);
  
end