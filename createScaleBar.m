function h_scale = createScaleBar(h_ax, h_scale, XOption, YOption, fontSz, XUnit, YUnit)
%% h_scale = createScaleBar(h_ax, h_scale, XOption, YOption, fontSz, XUnit, YUnit)
%
%  Function to automate scale bar creation. Requires scalebar.m @Chenxinfeng, 2016-9-10

% Set default values of parameters
if (nargin < 7) YUnit = '\muV'; end
if (nargin < 6) XUnit = 's'; end
if (nargin < 5) fontSz = 16; end
if (nargin < 4) YOption = 1; end
if (nargin < 3) XOption = 1; end

xScaleBarFactor = 15;
yScaleBarFactor = 10;

if isempty(h_scale)
  h_scale = scalebar(h_ax);
end

% x (timing)
if (XOption == 1)
  h_scale.XLen = round5sd((h_ax.XLim(2) - h_ax.XLim(1)) / xScaleBarFactor, 1, 'ceil');
  h_scale.XUnit = XUnit;
  h_scale.hTextX.FontSize = fontSz;
  h_scale.hLineX(1).Color = [0.1 0.1 0.1];
  h_scale.hLineX(2).Color = [0.1 0.1 0.1];
  h_scale.hTextX.Color = [0.1 0.1 0.1];
else
  h_scale.XLen = 0;
  h_scale.hTextX.String = '';
end

% y (data)
if(YOption == 1)
  h_scale.YLen = round5sd((h_ax.YLim(2) - h_ax.YLim(1)) / yScaleBarFactor, 1, 'ceil');
  h_scale.YUnit = YUnit;
  h_scale.hTextY.FontSize = fontSz;
  h_scale.hTextY_Rot = 0;
  h_scale.hTextY.HorizontalAlignment = 'right';
  h_scale.hLineY(1).Color = [0.1 0.1 0.1];
  h_scale.hLineY(2).Color = [0.1 0.1 0.1];
  h_scale.hTextY.Color = [0.1 0.1 0.1];
else
  h_scale.YLen = 0;
  h_scale.hTextY.String = '';
end

% Adjust Positions:
h_scale.Position = [h_ax.XLim(2) - h_scale.XLen, h_ax.YLim(2) - 2 * h_scale.YLen];
h_scale.hTextX_Pos = [0.3 * h_scale.XLen, -0.3 * h_scale.YLen];
h_scale.hTextY_Pos = [-0.05 * h_scale.XLen, 0.5 * h_scale.YLen];

end