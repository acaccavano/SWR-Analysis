function [baseMean, baseSD, baseThresh, peakThresh] = calcCaThresh(tSeries, timing, param)

if (nargin < 3) param = struct; end

if ~isfield(param,'sdMult')        param.sdMult        = 5;    end
if ~isfield(param,'sdBaseFactor')  param.sdBaseFactor  = 0.6; end
if ~isfield(param,'skipDetectLim') param.skipDetectLim = 2;    end

% Test if tSeries is cell array, if so it is array of files with consistent baselines
if iscell(tSeries)
  [~, nCols] = cellfun(@size,tSeries);
  
  if (all(nCols == nCols(1)))
    % Remove from threshold detection periods before param.skipDetectLim
    for i = 1:length(tSeries)
      tSeries{i} = tSeries{i}(timing{i} >= 1000 * param.skipDetectLim, :);
    end
    
    % concatenate into one continuous file:
    tSeries = vertcat(tSeries{:});
  else
    error('Array mismatch. Ensure all files have same number of cells to use consistent baseline');
  end
else
  % Remove from threshold detection periods before param.skipDetectLim
  tSeries = tSeries(timing >= 1000 * param.skipDetectLim, :);
end

% tSeries diminsion > 1 signifies multiple cell recording (Ca imaging)
nChannels  = size(tSeries, 2);
baseMean   = zeros(1, nChannels);
baseSD     = zeros(1, nChannels);
baseThresh = zeros(1, nChannels);
peakThresh = zeros(1, nChannels);

for ch = 1:nChannels
  [baseMean(ch), baseSD(ch)] = calcCaBaseline(tSeries(:,ch), param);
  baseThresh(ch) = baseMean(ch) + param.sdBaseFactor * baseSD(ch) * param.sdMult;
  peakThresh(ch) = baseMean(ch) + baseSD(ch) * param.sdMult;
end

end