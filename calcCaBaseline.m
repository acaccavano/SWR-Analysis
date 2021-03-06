function [baseMean, baseSD] = calcCaBaseline(tSeries, param)
%% [baseMean, baseSD] = calcCaBaseline(tSeries, param)
%
%  Function to determine baseline of individual cell tSeries
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   tSeries    = dFoF time series for single cell
%   param      = structure containing all parameters including:
%     param.baseDetectMethod = Method for baseline stats detection (0: none, 1: lower quantile, 2: iterative gaussian fitting (default))
%     param.baseQuant        = Lower quantile for baseline cutoff (default = 0.8)
%     param.pkDiffMin        = min distance between double gaussian peaks to consider them equivalent = abs(B1-B2) (default = 0.1)
%     param.pkSimLim         = Peak amplitude similarity metric = (A1^2 + A2^2)/(A1*A2) (default = 2)
%     param.kurtosisMin      = Min kurtosis limit to fit with 2 gaussians (otherwise skip 1st fit) (default = 0)
%     param.kurtosisMax      = Max kurtosis limit until exclude high points (otherwise fit can fail) (default = 5)
%     param.excludeQuant     = quantile above which to exclude if max kurtosis limit reached (default = 0.98)
%     param.plotFitHisto     = boolean option to plot histograms and fits for each cell
%
%  Outputs:
%   baseMean     = baseline mean
%   baseSD       = baseline SD

if (nargin < 2); param = struct; end
if ~isfield(param,'baseDetectMethod'); param.baseDetectMethod = 2;    end
if ~isfield(param,'baseQuant');        param.baseQuant        = 0.8;  end
if ~isfield(param,'pkDiffMin');        param.pkDiffMin        = 0.1;  end
if ~isfield(param,'pkSimLim');         param.pkSimLim         = 2;    end
if ~isfield(param,'kurtosisMin');      param.kurtosisMin      = 0;    end
if ~isfield(param,'kurtosisMax');      param.kurtosisMax      = 5;    end
if ~isfield(param,'excludeQuant');     param.excludeQuant     = 0.98; end
if ~isfield(param,'plotFitHisto');     param.plotFitHisto     = 0;    end

if (param.baseDetectMethod == 1)
  tSeries(tSeries > quantile(tSeries, param.baseQuant)) = [];
  baseMean = mean(tSeries);
  baseSD   = std(tSeries);
  
elseif (param.baseDetectMethod == 2)

  % Bin all data across tSeries:
  % Ensure  there are at least 200 bins per 1 dF/F:
  nBins = max(floor(200*range(tSeries)), 200);
  [binValues, binEdges] = histcounts(tSeries, nBins);
  binWidth = binEdges(2) - binEdges(1);
  
  % 1st pass: two-term gaussian to estimate baseline range (only if skewed dist)
  xVal = binEdges(1:length(binEdges) - 1)' + 0.5*binWidth;
  yVal = binValues';
  exRange = zeros(length(xVal), 1);
  
  if param.plotFitHisto
    figure
    hold on
  end

  if kurtosis(tSeries) > param.kurtosisMin
      
    ft = fittype('gauss2');
    loParams = [0 min(tSeries) 0 0 min(tSeries) 0];
    hiParams = [Inf max(tSeries) range(tSeries) Inf max(tSeries) range(tSeries)];
    weightEq = exp(-10*linspace(1,nBins,nBins)'/nBins);
    
    if kurtosis(tSeries) > param.kurtosisMax
      hiLimX = quantile(tSeries, param.excludeQuant);
      exRange = xVal > hiLimX;
    end
    
    f = fit(xVal, yVal, ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare', 'Weights', weightEq, 'Exclude', exRange);
    
    if param.plotFitHisto 
      plot(f, xVal, yVal, exRange) 
    end

    pkSimMetric   = sqrt((f.a1^2 + f.a2^2)/(f.a1*f.a2));
    
    % Special case where peaks are on top of eachother (often when noise >> signal):
    if abs(f.b1 - f.b2) < param.pkDiffMin
      % Use average of fits to define BL range:
      hiLimX = mean([f.b1 f.b2]) + mean([f.c1 f.c2])/sqrt(2);
      
      % Special case when peaks are similar amplitude (often when noise similar or < signal):
    elseif pkSimMetric < param.pkSimLim
      % Use more negative peak to define BL range:
      if (f.b1 < f.b2)
        hiLimX = f.b1 + f.c1/sqrt(2);
      else
        hiLimX = f.b2 + f.c2/sqrt(2);
      end
      
      % Standard case: use peak with greater amplitude (noise > signal):
    elseif (f.a1 > f.a2)
      hiLimX = f.b1 + f.c1/sqrt(2);
      
    elseif (f.a1 < f.a2)
      hiLimX = f.b2 + f.c2/sqrt(2);
      
    end
  end

  ft = fittype('gauss1');
  loParams = [0 min(tSeries) 0];
  hiParams = [Inf max(tSeries) range(tSeries)];
  exRange  = xVal > hiLimX;
  
  f = fit(xVal, yVal, ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare', 'Exclude', exRange);

  baseMean = f.b1;
  baseSD   = f.c1/sqrt(2);
  
  if param.plotFitHisto 
    plot(f, xVal, yVal, exRange) 
  end 
  
end

end