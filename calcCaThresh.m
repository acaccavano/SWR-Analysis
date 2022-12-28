function [baseMean, baseSD, baseThresh, peakThresh, hand] = calcCaThresh(tSeries, hand, param)
%% [baseMean, baseSD, baseThresh, peakThresh, hand] = calcCaThresh(tSeries, hand, param)
%
%  Script to detect thresholds for peak detection for array of cells for one or more files
%
%  Inputs: (all optional - will be prompted for or use defaults)
%   tSeries    = dFoF time series, either array of cells, or cell array of arrays for multiple files
%   hand       = handle structure to specify where figure should be drawn
%   param      = structure containing all parameters including:
%     param.sdMult        = SD of baseline for threshold detection (default = 4)
%     param.sdBaseFactor  = Factor of sdMult to consider for event start/end times (default = 0.75 eg 3SD)
%
%  Outputs:
%   baseMean     = baseline mean array
%   baseSD       = baseline SD array
%   baseThresh   = threshold array for event start/end times
%   peakThresh   = threshold array for event peak detection
%   hand         = handle structure to specify where figure should be drawn

%% Handle input arguments - if not entered
if (nargin < 3); param       = struct; end
if (nargin < 2); hand        = struct; end

% Handle case in which empty variables are supplied:
if isempty(param); param     = struct; end
if isempty(hand);  hand      = struct; end

if ~isfield(param,'sdMult');           param.sdMult           = 4;    end
if ~isfield(param,'sdBaseFactor');     param.sdBaseFactor     = 0.5; end

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

if param.plotFitHisto
  hand.blFig = figure('Name','Baseline Analysis','units','normalized','outerposition',[0.25 0.25 0.5 0.75]);
  hold on
end

for ch = 1:nChannels
  [baseMean(ch), baseSD(ch), hand] = calcBaseline(tSeries(:,ch), hand, param);
  baseThresh(ch) = baseMean(ch) + param.sdBaseFactor * baseSD(ch) * param.sdMult;
  peakThresh(ch) = baseMean(ch) + baseSD(ch) * param.sdMult;
end

end