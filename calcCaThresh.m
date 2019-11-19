function [baseMean, baseSD, baseThresh, peakThresh] = calcCaThresh(tSeries, param)

if (nargin < 2) param = struct; end

if ~isfield(param,'sdMult')        param.sdMult        = 4;    end
if ~isfield(param,'sdBaseFactor')  param.sdBaseFactor  = 0.8;  end

% Test if tSeries is cell array, if so it is array of files with consistent baselines
if iscell(tSeries)
  [~, nCols] = cellfun(@size,tSeries);
  
  % Concatenate into one continuous file if same number of cells
  if (all(nCols == nCols(1)))
    tSeries = vertcat(tSeries{:});
  else
    error('Array mismatch. Ensure all files have same number of cells to use consistent baseline');
  end
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