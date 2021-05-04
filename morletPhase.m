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