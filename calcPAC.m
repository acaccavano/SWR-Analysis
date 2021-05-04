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
nSample   = length(data.LFP.tSeries);
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

end


function y = morletPhase(f, s, Fs, width)
% function y = morletPhase(f, s, Fs, width)
%
% Returns a vector 'y' containing the instantaneous phase values of signal 
% 's' filtered for frequency 'f' (via convolution with a complex Morlet 
% wavelet described by 'width')
%
% Create a Morlet wavelet 'y' with frequency resolution 'f' and temporal 
% resolution 't'. The wavelet will be normalized so the total energy is 1.
% The 'width' defines the temporal and frequency resolution for the given
% centre frequency 'f' by determining the number of cycles of the wavelet
% itself (see Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997) or 
% Event-Related Potentials: A Methods Handbook, Handy (editor), MIT Press, 
% (2005))
%
% INPUTS:
% f - frequency to filter at 
% s - signal to be filtered
% Fs - sampling frequency of s
% width - parameter which defines the mother wavelet (this is then 
% scaled and translated in order to filter for different frequencies, 
% >= 5 is suggested, see Tallon-Baudry et al., J. Neurosci. 15, 722-734 
% (1997))
% 
% Code History:
% Author: Angela Onslow, May 2010, see PAC
% Modified from Ole Jensen, August 1998, see ENERGYVEC
% Adapted by Adam Caccavano, May 2021

dt = 1 / Fs;
sf = f / width;
st = 1 / (2*pi*sf); 
t  = -(width/2)*st : dt : (width/2)*st; 
A  = 1 / sqrt((st*sqrt(pi)));
m  = A * exp(-t.^2/(2*st^2)) .* exp(1i*2*pi*f.*t);
y = conv(s, m');
y = angle(y);
y = y(ceil(length(m)/2) : length(y) - floor(length(m)/2));

end



function y = morletAmp(f, s, Fs, width)
% function y = morletPhase(f, s, Fs, width)
%
% Returns a vector 'y' containing the instantaneous amplitude values of signal 
% 's' filtered for frequency 'f' (via convolution with a complex Morlet 
% wavelet described by 'width')
%
% Create a Morlet wavelet 'y' with frequency resolution 'f' and temporal 
% resolution 't'. The wavelet will be normalized so the total energy is 1.
% The 'width' defines the temporal and frequency resolution for the given
% centre frequency 'f' by determining the number of cycles of the wavelet
% itself (see Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997) or 
% Event-Related Potentials: A Methods Handbook, Handy (editor), MIT Press, 
% (2005))
%
% INPUTS:
% f - frequency to filter at 
% s - signal to be filtered
% Fs - sampling frequency of s
% width - parameter which defines the mother wavelet (this is then 
% scaled and translated in order to filter for different frequencies, 
% >= 5 is suggested, see Tallon-Baudry et al., J. Neurosci. 15, 722-734 
% (1997))
% 
% Code History:
% Author: Angela Onslow, May 2010, see PAC
% Modified from Ole Jensen, August 1998, see ENERGYVEC
% Adapted by Adam Caccavano, May 2021

dt = 1 / Fs;
sf = f / width;
st = 1 / (2*pi*sf); 
t  = -(width/2)*st : dt : (width/2)*st; 
A  = 1 / sqrt((st*sqrt(pi)));
m  = A * exp(-t.^2/(2*st^2)) .* exp(1i*2*pi*f.*t);
y = conv(s, m');
y = abs(y);
y = y(ceil(length(m)/2) : length(y) - floor(length(m)/2));

end


