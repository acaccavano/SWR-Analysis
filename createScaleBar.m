function h_scale = createScaleBar(h_ax, h_scale, XOption, YOption, fontSz, XUnit, YUnit)
%% h_scale = createScaleBar(h_ax, h_scale, XOption, YOption, fontSz, XUnit, YUnit)
%
%  Function to automate scale bar creation. Requires scalebar.m @Chenxinfeng, 2016-9-10

% Set default values of parameters
if (nargin < 7); YUnit = '\muV'; end
if (nargin < 6); XUnit = 's'; end
if (nargin < 5); fontSz = 16; end
if (nargin < 4); YOption = 1; end
if (nargin < 3); XOption = 1; end

xSBFactor = 15;
ySBFactor = 10;
xRange = h_ax.XLim(2) - h_ax.XLim(1);
yRange = h_ax.YLim(2) - h_ax.YLim(1);
xSBLength = round5sd(xRange / xSBFactor, 1, 'ceil');
ySBLength = round5sd(yRange / ySBFactor, 1, 'ceil');
lnWidth = 2;

if isempty(h_scale)
  h_scale = scalebar(h_ax);
end

h_scale.Border = 'LR';

% x (timing)
if (XOption == 1)
  h_scale.XLen = xSBLength;
  h_scale.XUnit = XUnit;
  h_scale.hTextX.FontSize = fontSz;
  h_scale.hLineX(1).LineStyle = '-';
  h_scale.hLineX(2).LineStyle = '-';
  h_scale.hLineX(1).LineWidth = lnWidth;
  h_scale.hLineX(2).LineWidth = lnWidth;
  h_scale.hLineX(1).Color = [0 0 0];
  h_scale.hLineX(2).Color = [0 0 0];
  h_scale.hTextX.Color = [0 0 0];
else
  h_scale.XLen = 0;
  h_scale.hTextX.String = '';
end

% y (data)
if(YOption == 1)
  h_scale.YLen = ySBLength;
  h_scale.YUnit = YUnit;
  h_scale.hTextY.FontSize = fontSz;
  h_scale.hTextY_Rot = 0;
  h_scale.hTextY.HorizontalAlignment = 'right';
  h_scale.hLineY(1).LineStyle = '-';
  h_scale.hLineY(2).LineStyle = '-';
  h_scale.hLineY(1).LineWidth = lnWidth;
  h_scale.hLineY(2).LineWidth = lnWidth;
  h_scale.hLineY(1).Color = [0 0 0];
  h_scale.hLineY(2).Color = [0 0 0];
  h_scale.hTextY.Color = [0 0 0];
else
  h_scale.YLen = 0;
  h_scale.hTextY.String = '';
end

% Adjust Positions:
h_scale.Position = [h_ax.XLim(2) - h_scale.XLen - round(0.01*xRange,1), h_ax.YLim(2) - 2*ySBLength];
h_scale.hTextX_Pos = [0.4 * xSBLength, -0.15 * ySBLength];
h_scale.hTextY_Pos = [- round(0.01*xRange,1), 0.5 * ySBLength];

end