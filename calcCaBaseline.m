function [baseMean, baseSD] = calcCaBaseline(tSeries, param) % [baseMean, baseSD]

if (nargin < 2) param = struct; end
if ~isfield(param,'baseDetectMethod') param.baseDetectMethod  = 2; end
if ~isfield(param,'baseQuant')   param.baseQuant   = 0.8; end
if ~isfield(param,'minBaseline') param.minBaseline = 0.1; end
if ~isfield(param,'plotHisto')   param.plotHisto   = 0;   end

if (param.baseDetectMethod == 1)
  tSeries(tSeries > quantile(tSeries, param.baseQuant)) = [];
  baseMean = mean(tSeries);
  baseSD   = std(tSeries);
  
elseif (param.baseDetectMethod == 2)
  
  % Fitting parameter metrics
  pkDiffMin    = 0.1; % How close peaks can be before considered eqivalent = abs(B1-B2)
  pkSimLim     = 2; % Peak similarity metric = (A1^2 + A2^2)/(A1*A2)
  kurtosisMax  = 5; % Max kurtosis limit until exclude high points (otherwise fit can fail)
  kurtosisMin  = 0; % Min kurtosis limit to fit with 2 gaussians (otherwise skip 1st fit)
  excludeQuant = 0.98; % quantile to limit 1st fit to if kurtosis limit exceeded
  
  % Bin all data across tSeries:
  
  % Ensure  there are at least 200 bins per 1 dF/F:
  nBins = max(floor(200*range(tSeries)), 200);
  [binValues, binEdges] = histcounts(tSeries, nBins);
  binWidth = binEdges(2) - binEdges(1);
  
  % 1st pass: two-term gaussian to estimate baseline range (only if skewed dist)
  xVal = binEdges(1:length(binEdges) - 1)' + 0.5*binWidth;
  yVal = binValues';
  exRange = zeros(length(xVal), 1);
  
  if param.plotHisto
    figure
    hold on
  end

  if kurtosis(tSeries) > kurtosisMin
      
    ft = fittype('gauss2');
    loParams = [0 min(tSeries) 0 0 min(tSeries) 0];
    hiParams = [Inf max(tSeries) range(tSeries) Inf max(tSeries) range(tSeries)];
    weightEq = exp(-10*linspace(1,nBins,nBins)'/nBins);
    
    if kurtosis(tSeries) > kurtosisMax
      hiLimX = quantile(tSeries, excludeQuant);
      exRange = xVal > hiLimX;
    end
    
    f = fit(xVal, yVal, ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare', 'Weights', weightEq, 'Exclude', exRange);
    
    if param.plotHisto 
      plot(f, xVal, yVal, exRange) 
    end

    pkSimMetric   = sqrt((f.a1^2 + f.a2^2)/(f.a1*f.a2));
    
    % Special case where peaks are on top of eachother (often when noise >> signal):
    if abs(f.b1 - f.b2) < pkDiffMin
      % Use average of fits to define BL range:
      hiLimX = mean([f.b1 f.b2]) + mean([f.c1 f.c2])/sqrt(2);
      
      % Special case when peaks are similar amplitude (often when noise similar of < signal):
    elseif pkSimMetric < pkSimLim
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
  
  if param.plotHisto 
    plot(f, xVal, yVal, exRange) 
  end 
  
end

end