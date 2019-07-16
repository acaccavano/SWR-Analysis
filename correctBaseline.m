function tSeries = correctBaseline(timing, tSeries, param)

if (nargin < 3) param = struct; end

if ~isfield(param,'baseCorrectMethod')  param.baseCorrectMethod = 2;    end
if ~isfield(param,'CaFiltLim1')         param.CaFiltLim1        = 0.03; end
if ~isfield(param,'CaFiltLim2')         param.CaFiltLim2        = 4;    end
if ~isfield(param,'CaFiltOrder')        param.CaFiltOrder       = 80;   end
if ~isfield(param,'CaFiltAlpha')        param.CaFiltAlpha       = 2.5;  end
if ~isfield(param,'smoothFactor')       param.smoothFactor      = 0.25; end

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