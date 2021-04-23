function hand = plotPSCAnalysis(data, hand, param, dsPlot)
%% hand = plotPSCAnalysis(data, hand, param, dsPlot)
% 
%  Function to plot output of analyzePSCFile

nData = length(data);

% Set default parameters if not specified
if ~isfield(param,'useSWRDurationOption'); param.useSWRDurationOption = 0; end
if ~isfield(param,'useSWRWindowOption');   param.useSWRWindowOption = 1; end

% Plot dimension parameters
convFact = 1000;
marginSz = 0.04;
rasterSz = 0.04;
spacerSz = 0.02;
markerSz = 2;
fontSz   = 5 * (5 - nData);
tFact    = 0.1;
lnWidth  = 0.5;
axWidth  = (1 - (nData + 1) * marginSz) / nData;
axSz     = (1 - 2*marginSz - 2*rasterSz - 4*spacerSz)/2;

% Plot colors
% lfpCol  = [0.2 0.2 0.5];
% cellCol = [0.5 0.2 0.2];
% swrCol  = [0.5 0.5 0.8];
% swrCCol = [0.3 0.0 0.0];
% pscCol  = [0.8 0.5 0.5];
% pscCCol = [0.0 0.0 0.3];

lfpCol{1}  = [ 48  70 160]/255;
lfpCol{2}  = [ 50  50  50]/255;
% cellCol    = [  0  90   0]/255; % PC
cellCol    = [128   0   0]/255; % PVBC
% cellCol    = [ 65   2  87]/255; % PVBSC
% cellCol    = [192  96   0]/255; % PVAAC
swrCol     = [180 180 180]/255;
swrCCol{1} = lfpCol{1};
swrCCol{2} = lfpCol{2};
pscCol     = swrCol;
pscCCol{1} = lfpCol{1};
pscCCol{2} = lfpCol{2};

% Initialize graphical structures
hand.axTrLFP   = gobjects(nData, 1);
hand.axRsLFP   = gobjects(nData, 1);
hand.axTrCell  = gobjects(nData, 1);
hand.axRsCell  = gobjects(nData, 1);
hand.lblLFP    = gobjects(nData, 1);
hand.lblCell   = gobjects(nData, 1);
hand.scaleLFP  = struct;
hand.scaleCell = struct;

% Initialize plot data
for i = 1:nData
  timing{i}  = downsampleMean(data(i).LFP.timing/1000, dsPlot);
  trLFP{i}   = downsampleMean(convFact * data(i).LFP.tSeries, dsPlot);
  trCell{i}  = downsampleMean(data(i).C.tSeries, dsPlot);
  
  if param.useSWRDurationOption
    rsSWR{i} = downsampleMax(data(i).SWR.evStatus, dsPlot);
  elseif param.useSWRWindowOption
    rsSWR{i} = downsampleMax(data(i).SWR.evStatusStand, dsPlot);
  end
  
  rsSWRC{i}  = downsampleMax(data(i).SWR.PSC.evStatusC, dsPlot);
  rsPSC{i}   = downsampleMax(data(i).C.PSC.evStatusPeak, dsPlot);
  rsPSCC{i}  = downsampleMax(data(i).C.PSC.evStatusPeakC, dsPlot);
end

%% In reverse y order, initialize trace/raster plots
yPos = marginSz;

% In forward x order, initialize cell raster(s)
xPos = marginSz;
for i = 1:nData
  hand.axRsCell(i) = subplot('Position',[xPos yPos axWidth rasterSz]);
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + rasterSz + spacerSz;

% In forward x order, initialize cell trace(s)
xPos = marginSz;
for i = 1:nData
  hand.axTrCell(i) = subplot('Position',[xPos yPos axWidth axSz]);
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + axSz + spacerSz;

% In forward x order, initialize SWR raster(s)
xPos = marginSz;
for i = 1:nData
  hand.axRsLFP(i) = subplot('Position',[xPos yPos axWidth rasterSz]);
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + rasterSz + spacerSz;

% In forward x order, initialize LFP trace(s)
xPos = marginSz;
for i = 1:nData
  hand.axTrLFP(i) = subplot('Position',[xPos yPos axWidth axSz]);
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + axSz + spacerSz;

%% LFP Plot
minY =  999999;
maxY = -999999;

for i = 1:nData
  plot(hand.axTrLFP(i), timing{i}, trLFP{i}, 'LineWidth', lnWidth, 'Color', lfpCol{i});
  minY = min(minY, min(trLFP{i}));
  maxY = max(maxY, max(trLFP{i}));
end

% Set uniform axes ranges
for i = 1:nData
  axis(hand.axTrLFP(i), [timing{i}(1) timing{i}(length(timing{i})) minY maxY]);
  axis(hand.axTrLFP(i), 'off');
  hand.lblLFP(i) = text(hand.axTrLFP(i), hand.axTrLFP(i).XLim(1), hand.axTrLFP(i).YLim(2) - tFact * (hand.axTrLFP(i).YLim(2) - hand.axTrLFP(i).YLim(1)), 'Local Field Potential and SWR Events', 'FontSize', fontSz, 'Color', [0 0 0]);
  hand.scaleLFP.(['scale' int2str(i)]) = createScaleBar(hand.axTrLFP(i), [], 1, 1, fontSz);
end

%% SWR raster plot
for i = 1:nData
  
  hold(hand.axRsLFP(i), 'on');
  
  % Assign SWR points
  if (sum(rsSWR{i}) == 0)
    swrPoints  = timing{i};
    swrChannel = NaN * ones(1,20);
  else
    swrPoints  = timing{i}(find(rsSWR{i} .* timing{i}));
    swrChannel = (1:20);
  end
  swrChannel = swrChannel(ones(1,length(swrPoints)),:);
  plot(hand.axRsLFP(i), swrPoints, swrChannel, '.', 'MarkerSize', markerSz, 'Color', swrCol);
  
  % Assign coincident SWR points
  if (sum(rsSWRC{i}) == 0)
    swrCPoints  = timing{i};
    swrCChannel = NaN * ones(1,20);
  else
    swrCPoints  = timing{i}(find(rsSWRC{i} .* timing{i}));
    swrCChannel = (1:20);
  end
  swrCChannel = swrCChannel(ones(1,size(swrCPoints,1)),:);
  plot(hand.axRsLFP(i), swrCPoints, swrCChannel,'.','MarkerSize', markerSz, 'Color', swrCCol{i});
  
  axis(hand.axRsLFP(i), [timing{i}(1) timing{i}(length(timing{i})) 1 20]);
  axis(hand.axRsLFP(i), 'off');
  hold(hand.axRsLFP(i), 'off');
  
end

%% Cell Plot
minY =  999999;
maxY = -999999;

for i = 1:nData
  plot(hand.axTrCell(i), timing{i}, trCell{i}, 'LineWidth', lnWidth, 'Color', cellCol);
  minY = min(minY, min(trCell{i}));
  maxY = max(maxY, max(trCell{i}));
end

% Set uniform axes ranges:
for i = 1:nData
  axis(hand.axTrCell(i), [timing{i}(1) timing{i}(length(timing{i})) minY maxY]);
  axis(hand.axTrCell(i), 'off');
  hand.lblCell(i) = text(hand.axTrCell(i), hand.axTrCell(i).XLim(1), hand.axTrCell(i).YLim(2) - tFact * (hand.axTrCell(i).YLim(2) - hand.axTrCell(i).YLim(1)), 'Post-Synaptic Currents', 'FontSize', fontSz, 'Color', [0 0 0]);
  hand.scaleCell.(['scale' int2str(i)]) = createScaleBar(hand.axTrCell(i), [], 1, 1, fontSz, 's', 'pA');
end

%% PSC raster plot
for i = 1:nData
  
  hold(hand.axRsCell(i), 'on');
  
  % Assign PSC points
  if (sum(rsPSC{i}) == 0)
    pscPoints  = timing{i};
    pscChannel = NaN * ones(1,20);
  else
    pscPoints  = timing{i}(find(rsPSC{i} .* timing{i}));
    pscChannel = (1:20);
  end
  pscChannel = pscChannel(ones(1,length(pscPoints)),:);
  plot(hand.axRsCell(i), pscPoints, pscChannel, '.', 'MarkerSize', markerSz, 'Color', pscCol);
  
  % Assign coincident PSC points
  if (sum(rsPSCC{i}) == 0)
    pscCPoints  = timing{i};
    pscCChannel = NaN * ones(1,20);
  else
    pscCPoints  = timing{i}(find(rsPSCC{i} .* timing{i}));
    pscCChannel = (1:20);
  end
  pscCChannel = pscCChannel(ones(1,size(pscCPoints,1)),:);
  plot(hand.axRsCell(i), pscCPoints, pscCChannel,'.','MarkerSize', markerSz, 'Color', pscCCol{i});
  
  axis(hand.axRsCell(i), [timing{i}(1) timing{i}(length(timing{i})) 1 20]);
  axis(hand.axRsCell(i), 'off');
  hold(hand.axRsCell(i), 'off');
  
end

%% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
set(pan(hand.pscFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});
set(zoom(hand.pscFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});

for i = 1:nData
  linkaxes([hand.axTrLFP(i), hand.axRsLFP(i), hand.axTrCell(i), hand.axRsCell(i)], 'x');
end

end

% Callback function to redraw scalebar and labels after zooming/panning
function redrawPZCallback(obj, evd, hand, fontSz, tFact)
scaleNames = fieldnames(hand.scaleLFP);
for i = 1:length(hand.lblLFP)
  hand.lblLFP(i).Position  = [hand.axTrLFP(i).XLim(1), hand.axTrLFP(i).YLim(2) - tFact * (hand.axTrLFP(i).YLim(2) - hand.axTrLFP(i).YLim(1))];
  hand.lblCell(i).Position = [hand.axTrCell(i).XLim(1), hand.axTrCell(i).YLim(2) - tFact * (hand.axTrCell(i).YLim(2) - hand.axTrCell(i).YLim(1))];
  createScaleBar(hand.axTrLFP(i), hand.scaleLFP.(scaleNames{i}), 1, 1, fontSz);
  createScaleBar(hand.axTrCell(i), hand.scaleCell.(scaleNames{i}), 1, 1, fontSz, 's', 'pA');
end
end

