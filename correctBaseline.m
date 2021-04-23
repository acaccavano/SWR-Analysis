function tSeries = correctBaseline(timing, tSeries, param)
%% tSeries = correctBaseline(timing, tSeries, param)
%
%  Function to subtract slowly drifting baseline for individual dFoF time
%  series for one cell
%
%  Inputs: (param optional)
%   timing     = array of timing
%   tSeries    = dFoF time series for single cell
%   param      = structure containing all parameters including:
%     param.baseCorrectMethod    = Method for baseline correction (0: none, 1: gassuian filter, 2: smoothed average (default))
%     param.CaFiltLim1           = Lower limit for gaussian filter (default = 0.03Hz)
%     param.CaFiltLim2           = Upper limit for gaussian filter (default = 4Hz)
%     param.CaFiltOrder          = Gaussian filter order (default = 80)
%     param.CaFiltAlpha          = Gaussian filter alpha (default = 2.5)
%     param.smoothFactor         = Proportion of file duration for moving linear average (default = 0.25)
%
%  Output:
%   tSeries    = dFoF time series for single cell with corrected baseline

if (nargin < 3); param = struct; end

if ~isfield(param,'baseCorrectMethod');  param.baseCorrectMethod = 2;    end
if ~isfield(param,'CaFiltLim1');         param.CaFiltLim1        = 0.03; end
if ~isfield(param,'CaFiltLim2');         param.CaFiltLim2        = 4;    end
if ~isfield(param,'CaFiltOrder');        param.CaFiltOrder       = 80;   end
if ~isfield(param,'CaFiltAlpha');        param.CaFiltAlpha       = 2.5;  end
if ~isfield(param,'smoothFactor');       param.smoothFactor      = 0.25; end

if (param.baseCorrectMethod == 1)
  samplingInt = timing(2) - timing(1);
  tSeries = gaussianFilt(tSeries, param.CaFiltLim1, param.CaFiltLim2, samplingInt, param.CaFiltOrder, param.CaFiltAlpha);
  
elseif (param.baseCorrectMethod == 2)
  tSmooth = smooth(timing, tSeries, param.smoothFactor, 'rloess');
  tSeries = tSeries - tSmooth;
end

% % re-adjust baseline to ensure no negative dFoF values
% tSeries = tSeries - min(tSeries);

end