function [LFP, S1, S2] = calcPAC(LFP, S1, S2, param)

% Below is an adaptation of code written by Author: Angela Onslow, May 2010

% Parameters
if (nargin < 4); param = struct; end
if ~isfield(param,'morlWidth');  param.morlWidth  = 7;   end
if ~isfield(param,'winLength');  param.winLength  = 0.5; end
if ~isfield(param,'winOverlap'); param.winOverlap = 0.2; end

tSeries   = LFP.tSeries;

% Calculate modulating phase (lower frequency) via Morlet wavelet:
morlFreqP = S1.lim1 + floor((S1.lim2 - S1.lim1)/2);
modPhase  = morletPhase(morlFreqP, tSeries, param.Fs, param.morlWidth);

% Calculate amplitude (higher frequency) via Morlet wavelet:
morlFreqA = S2.lim1 + floor((S2.lim2 - S2.lim1)/2);
pacAmp    = morletAmp(morlFreqA, tSeries, param.Fs, param.morlWidth);

% Calculate total PAC measure:
z = pacAmp .* exp(1i * modPhase); % Create composite signal
miRaw = mean(z);  % Compute the mean length of composite signal
miPAC = abs(miRaw);

% Truncate signals to get integer number of time windows
nSample   = length(LFP.tSeries);
nWin      = ceil(param.winLength * param.Fs);
nOverlap  = ceil(param.winOverlap * param.Fs);
remSample = mod(nSample, nWin);
tSeries   = tSeries(1 : nSample - remSample);
modPhase  = modPhase(1 : nSample - remSample);
pacAmp    = pacAmp(1 : nSample - remSample);

% Update nSample
nSample = length(tSeries);
idx     = bsxfun(@plus, (1:nWin)', 1+(0:(fix((nSample - nOverlap)/(nWin - nOverlap)) - 1))*(nWin - nOverlap)) - 1;
%   xbins   = size(idx, 2);
%   alpha = alpha/(xbins*ybins); % Uncomment to use Bonferonni Correction

% Calculate windowed time-series PAC:
miPACWin = zeros(size(idx,2), 1);
for i = 1:size(idx, 2)
  
  z = pacAmp(idx(:,i)) .* exp(1i * modPhase(idx(:,i))); % Create composite signal
  miRaw = mean(z);  % Compute the mean length of composite signal
  miPACWin(i) = abs(miRaw);
  
end

% Add to data structure - TEMP FIX THIS
data.xFreq = struct;
data.xFreq.miPAC = miPAC;
data.xFreq.miPACWin = miPACWin;
data.xFreq.morlFreqP = morlFreqP;
data.xFreq.morlFreqA = morlFreqA;
data.xFreq.modPhase = modPhase;
data.xFreq.pacAmp = pacAmp;
data.xFreq.idx = idx;






