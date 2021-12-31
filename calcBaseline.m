function [baseMean, baseSD] = calcBaseline(tSeries, param)
%% [baseMean, baseSD] = calcBaseline(tSeries, param)
%
%  Function to determine baseline noise of signal with true upward events.
%  Assumes signal is normalized (RMS or dFoF)
%
%  Inputs: (will use defaults if not supplied)
%   tSeries    = dFoF time series for single cell
%   param      = structure containing all parameters including:
%     param.baseDetectMethod = Method for baseline stats detection (0: none, 1: lower quantile, 2: iterative gaussian fitting (default))
%     param.baseQuant        = Lower quantile for baseline cutoff (default = 0.80)
%     param.pkDiffMin        = min fractional distance between double gaussian peaks to consider them equivalent = abs(B1-B2) (default = 0.1)
%     param.pkSimLim         = Peak amplitude similarity metric = (A1^2 + A2^2)/(A1*A2) (default = 2)
%     param.kurtosisMin      = Min kurtosis limit to fit with 2 gaussians (otherwise skip 1st fit) (default = 0)
%     param.kurtosisMax      = Max kurtosis limit until exclude high points (otherwise fit can fail) (default = 5)
%     param.excludeQuant     = quantile above which to exclude if max kurtosis limit reached (default = 0.95)
%     param.plotFitHisto     = boolean option to plot histograms and fits for the signal
%     param.minBins          = quantile above which to exclude if max kurtosis limit reached (default = 0.98)
%     param.skewedBL         = boolean option to indicate skewed BL distribution, and use both gaussians just for BL 
%
%  Outputs:
%   baseMean     = baseline mean
%   baseSD       = baseline SD

if (nargin < 2); param = struct; end
if ~isfield(param,'baseDetectMethod'); param.baseDetectMethod = 2;    end
if ~isfield(param,'baseQuant');        param.baseQuant        = 0.80; end
if ~isfield(param,'pkDiffMin');        param.pkDiffMin        = 0.1;  end
if ~isfield(param,'pkSimLim');         param.pkSimLim         = 2;    end
if ~isfield(param,'kurtosisMin');      param.kurtosisMin      = 0;    end
if ~isfield(param,'kurtosisMax');      param.kurtosisMax      = 5;    end
if ~isfield(param,'excludeQuant');     param.excludeQuant     = 0.95; end
if ~isfield(param,'plotFitHisto');     param.plotFitHisto     = 1;    end
if ~isfield(param,'minBins');          param.minBins          = 500;  end
if ~isfield(param,'skewedBL');         param.skewedBL         = 0;    end
if ~isfield(param,'weightFact');       param.weightFact       = 5;    end
if ~isfield(param,'nSDLim');           param.nSDLim           = 2;    end

% Calc simple baseline for both methods:
baseline = tSeries;
baseline(baseline > quantile(baseline, param.baseQuant)) = [];
  
if (param.baseDetectMethod == 1)
  baseMean = mean(baseline);
  baseSD   = std(baseline);
  
elseif (param.baseDetectMethod == 2)

  % Bin all data across tSeries:
  nBins = floor(param.minBins * range(tSeries) / range(baseline));
  [binValues, binEdges] = histcounts(tSeries, nBins);
  binWidth = binEdges(2) - binEdges(1);
    
  xVal = binEdges(1:length(binEdges) - 1)' + 0.5*binWidth;
  yVal = binValues';
  exRange = zeros(length(xVal), 1);
  
  if param.plotFitHisto
    figure
    hold on
  end
  
  % Rough estimate of peaks:
  [~, xPeak1] = maxk(yVal, 10); % Finds indeces of 10-highest yVals
  xPeak1 = xVal(round(mean(xPeak1))); % To get approximate index of 1st peak
  xPeak2 = xPeak1 + param.pkDiffMin * range(xVal); % Estimate location of 2nd peak by param.pkDiffMin (if 0, same as xPeak1)
  
  % 1st iteration: 2-term gaussian to estimate baseline range (only use if min kurtosis limit reached)
  if (kurtosis(tSeries) > param.kurtosisMin)

    % Double gaussian model:
    ft = fittype('gauss2');

    loParams = [0 max(0.5*(xPeak1 + min(xVal)), min(xVal)) 0 0 max(0.5*(xPeak2 + min(xVal)), min(xVal)) 0];
    hiParams = [Inf min(0.5*(3*xPeak1 - min(xVal)), max(xVal)) range(xVal) Inf min(0.5*(3*xPeak2 - min(xVal)), max(xVal)) range(xVal)];

    % Exclude extreme kurtosis values to help with fitting:
    if kurtosis(tSeries) > param.kurtosisMax
      hiLimX = quantile(tSeries, param.excludeQuant);
      exRange = xVal > hiLimX;
    end
    
    % Weight lower values more within inclusion range:
    weightEq = exp(-param.weightFact*linspace(1,nBins,nBins)'/(nBins-sum(exRange)));
    
    f = fit(xVal, yVal, ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare', 'Weights', weightEq, 'Exclude', exRange);
    
    if param.plotFitHisto 
      yyaxis left
      plot(f, xVal, yVal, exRange) 
      yyaxis right
      plot(xVal, weightEq, 'DisplayName','weight')
    end

    pkSimMetric  = sqrt((f.a1^2 + f.a2^2)/(f.a1*f.a2));
    pkDiffMetric = abs(f.b1 - f.b2)/range(xVal);
    
    % If skewed BL distribution, use fits of right-side of distribution to estimate SD
    if param.skewedBL
      baseMean = mean([f.b1 f.b2]);
      baseSD   = sqrt((f.c1*f.c1 + f.c2*f.c2)/2);
      
    % Special case where peaks are on top of eachother (often when noise >> signal):
    elseif (pkDiffMetric < param.pkDiffMin)
      % Use average of fits to define BL range:
      loLimX = mean([f.b1 f.b2]) - param.nSDLim * sqrt((f.c1*f.c1 + f.c2*f.c2)/2);
      hiLimX = mean([f.b1 f.b2]) + param.nSDLim * sqrt((f.c1*f.c1 + f.c2*f.c2)/2);
      
    % Special case when peaks are similar amplitude (often when noise similar or < signal):
    elseif pkSimMetric < param.pkSimLim
      % Use more negative peak to define BL range:
      if (f.b1 < f.b2)
        loLimX = f.b1 - param.nSDLim * f.c1/sqrt(2);
        hiLimX = f.b1 + param.nSDLim * f.c1/sqrt(2);
      else
        loLimX = f.b2 - param.nSDLim * f.c2/sqrt(2);
        hiLimX = f.b2 + param.nSDLim * f.c2/sqrt(2);
      end
      
    % Standard case: use peak with greater amplitude (noise > signal):
    elseif (f.a1 > f.a2)
      loLimX = f.b1 - param.nSDLim * f.c1/sqrt(2);
      hiLimX = f.b1 + param.nSDLim * f.c1/sqrt(2);
    elseif (f.a1 < f.a2)
      loLimX = f.b2 - param.nSDLim * f.c2/sqrt(2);
      hiLimX = f.b2 + param.nSDLim * f.c2/sqrt(2);
    end
    
  % Single Gaussian Model:
  else
    ft = fittype('gauss1');
    loParams = [0 max(0.5*(xPeak1 + min(xVal)), min(xVal)) 0];
    hiParams = [Inf min(0.5*(3*xPeak1 - min(xVal)), max(xVal)) range(xVal)];
    f = fit(xVal, yVal, ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare');
    baseMean = f.b1;
    baseSD   = f.c1/sqrt(2);
    if param.plotFitHisto
      yyaxis left
      plot(f, xVal, yVal, exRange)
    end
  end
  
  % 2nd iteration (for non-skewed BL distributions): single gaussian fit restricted to identified hiLimX:
  if ~param.skewedBL
    ft = fittype('gauss1');
    loParams = [0 min(tSeries) 0];
    hiParams = [Inf max(tSeries) range(tSeries)];
    exRange  = (xVal < loLimX) | (xVal > hiLimX);
    f = fit(xVal, yVal, ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare', 'Exclude', exRange);
    baseMean = f.b1;
    baseSD   = f.c1/sqrt(2);
    if param.plotFitHisto
      yyaxis left
      plot(f, xVal, yVal, exRange)
    end
  end
end
