function xFreq = calcPACMI(xFreq, Fs, nShuffle)
% Computes normalized phase-amplitude coupled modulation index
% Adapted from Canolty et al 2006 and Onslow et al 2011

nSample = size(xFreq.phsPAC, 1);
nPhs    = size(xFreq.phsPAC, 2);
nAmp    = size(xFreq.ampPAC, 2);

% Initialize output arrays:
xFreq.pacMI       = zeros(nAmp, nPhs);
xFreq.pacMI_Len   = zeros(nAmp, nPhs);
xFreq.pacMI_Phase = zeros(nAmp, nPhs);

if nShuffle > 1
  xFreq.pacMIShf_LenAve = zeros(nAmp, nPhs);
  xFreq.pacMIShf_LenSTD = zeros(nAmp, nPhs);
  xFreq.pacMI_Z         = zeros(nAmp, nPhs);
  xFreq.pacMI_LenZ      = zeros(nAmp, nPhs);
end

for i = 1 : nAmp
  for j = 1 : nPhs
    z = xFreq.ampPAC(:,i) .* exp(1i * xFreq.phsPAC(:,j)); % Create composite signal
    xFreq.pacMI(i,j)       = mean(z);  % Compute the mean length of composite signal
    xFreq.pacMI_Len(i,j)   = abs(xFreq.pacMI(i,j));
    xFreq.pacMI_Phase(i,j) = angle(xFreq.pacMI(i,j));
    
    if nShuffle > 1
      pacMIShf_Len = zeros(nShuffle, 1);
      
      % Compute shuffled values:
      for k = 1 : nShuffle
        skip = shuffleInd(nSample, Fs, nShuffle);
        ampPacShf = vertcat(xFreq.ampPAC(skip(k):end,j), xFreq.ampPAC(1:skip(k)-1,j));
        pacMIShf_Len(k) = abs(mean(ampPacShf .* exp(1i * xFreq.phsPAC(:,i))));
      end
      
      % Fit gaussian to shuffled data, uses normfit.m from MATLAB Statistics toolbox
      [xFreq.pacMIShf_LenAve(i,j), xFreq.pacMIShf_LenSTD(i,j)] = normfit(pacMIShf_Len);
      
      % Normalize length using shuffled data (z-score)
      xFreq.pacMI_LenZ(i,j) = (xFreq.pacMI_Len(i,j) - xFreq.pacMIShf_LenAve(i,j)) / xFreq.pacMIShf_LenSTD(i,j);
      xFreq.pacMI_Z(i,j) = xFreq.pacMI_LenZ(i,j) * exp(1i * xFreq.pacMI_Phase(i,j));
    end
  end
end