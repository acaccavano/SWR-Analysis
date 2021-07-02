function hand = plotSpecAnalysis(data, hand, param)
%% hand = plotSpecAnalysis(data, hand, param, dsPlot)
%
%  Function to plot spectral output of analyzeLFPFile (FFT and PAC)

%% Handle input arguments - if not entered
if (nargin < 3); param  = struct;  end
if (nargin < 2); hand   = struct;  end
if (nargin < 1); data   = struct;  end

% Handle case in which empty variables are supplied:
if isempty(param); param = struct; end
if isempty(hand);  hand  = struct; end
if isempty(data);  data  = struct; end

% Set default parameters if not specified
if ~isfield(param,'fftOption');   param.fftOption   = 1; end
if ~isfield(param,'xFreqOption'); param.xFreqOption = 1; end
if ~isfield(param,'spectLim1');   param.spectLim1   = 1; end
if ~isfield(param,'spectLim2');   param.spectLim2   = 100; end

% Plotting parameters
nData      = length(data);
lnWidth    = 1.5;
marginSz   = 0.1;
fontSz     = 12;
fitColor   = [1 0 0];
plotHeight = (1 - 3*marginSz) / 2;
plotWidth  = (1 - (nData + 1)*marginSz) / nData;
xPos       = marginSz;
minC       =  999999;
maxC       = -999999;

hand.axFFT = gobjects(nData, 1);
hand.axPAC = gobjects(nData, 1);

for i = 1:nData
  
  % FFT
  xLFP = data(i).LFP.FFT.fftRange(data(i).LFP.FFT.subRange);
  yLFP = data(i).LFP.FFT.F1(data(i).LFP.FFT.subRange);
  
  hand.axFFT(i) = subplot('Position',[xPos 2*marginSz+plotHeight plotWidth plotHeight]);
  bar(hand.axFFT(i), xLFP, yLFP, 1)
  hold on
  
  % Overlay gaussian fits (if present)
  if isfield(data(i), 'theta')
    xTheta     = data(i).theta.FFT.fftRange(data(i).theta.FFT.subRange);
    yThetaFit  = data(i).theta.FFT.fitGauss(xTheta);
    plot(hand.axFFT(i), xTheta, yThetaFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  if isfield(data(i), 'alpha')
    xAlpha     = data(i).alpha.FFT.fftRange(data(i).alpha.FFT.subRange);
    yAlphaFit  = data(i).alpha.FFT.fitGauss(xAlpha);
    plot(hand.axFFT(i), xAlpha, yAlphaFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  if isfield(data(i), 'beta')
    xBeta     = data(i).beta.FFT.fftRange(data(i).beta.FFT.subRange);
    yBetaFit  = data(i).beta.FFT.fitGauss(xGamma);
    plot(hand.axFFT(i), xBeta, yBetaFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  if isfield(data(i), 'gamma')
    xGamma     = data(i).gamma.FFT.fftRange(data(i).gamma.FFT.subRange);
    yGammaFit  = data(i).gamma.FFT.fitGauss(xGamma);
    plot(hand.axFFT(i), xGamma, yGammaFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  if isfield(data(i), 'hgamma')
    xHGamma    = data(i).hgamma.FFT.fftRange(data(i).hgamma.FFT.subRange);
    yHGammaFit = data(i).hgamma.FFT.fitGauss(xHGamma);
    plot(hand.axFFT(i), xHGamma, yHGammaFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  if isfield(data(i), 'R')
    xRipple    = data(i).R.FFT.fftRange(data(i).R.FFT.subRange);
    yRippleFit = data(i).R.FFT.fitGauss(xRipple);
    plot(hand.axFFT(i), xRipple, yRippleFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  if isfield(data(i), 'fR')
    xFRipple    = data(i).fR.FFT.fftRange(data(i).fR.FFT.subRange);
    yFRippleFit = data(i).fR.FFT.fitGauss(xFRipple);
    plot(hand.axFFT(i), xFRipple, yFRippleFit, '-', 'Color', fitColor, 'LineWidth', lnWidth);
  end
  
  axis(hand.axFFT(i), [param.spectLim1 param.spectLim2 0 inf]);
  set(hand.axFFT(i), 'FontSize', fontSz);
  title(hand.axFFT(i), 'FFT')
  %   xlabel(hand.axFFT(i), 'Frequency (Hz)')
  
  % Phase-Amplitude Coupling
  hand.axPAC(i) = subplot('Position',[xPos marginSz plotWidth plotHeight]);

  if isfield(data(i).LFP.xFreq, 'pacMI_LenZ')
    pacMI  = data(i).LFP.xFreq.pacMI_LenZ;
    pacLbl = 'Phase-Amplitude Modulation Index Z-Score';
  else
    pacMI = data(i).LFP.xFreq.pacMI_Len;
    pacLbl = 'Phase-Amplitude Modulation Index';
  end
  imagesc(hand.axPAC(i), data(i).LFP.xFreq.morlFreq, data(i).LFP.xFreq.morlFreq, pacMI)
  
  axis xy
  title(hand.axPAC(i), pacLbl)
  ylabel(hand.axPAC(i), 'Amplitude Freq. (Hz)')
  xlabel(hand.axPAC(i), 'Phase Freq. (Hz)')
  set(hand.axPAC(i), 'FontSize', fontSz);
  
  minC = min(minC, min(min(pacMI)));
  maxC = max(maxC, max(max(pacMI)));
  xPos = xPos + plotWidth + marginSz;
  
  linkaxes([hand.axFFT(i), hand.axPAC(i)], 'x');
end

% Colorbar only on furthest right plot:
colorbar(hand.axPAC(nData), 'Position', [1.01-marginSz marginSz 0.02 plotHeight]);
colormap hot
for i = 1:nData
  caxis(hand.axPAC(i), [minC maxC]);
end


