function hand = plotSpkAnalysis(data, hand, param, dsPlot)

nData = length(data);

% Set default parameters if not specified
if ~isfield(param,'procBstOption') param.procBstOption = 1; end

% Plot dimension parameters
convFact = 1000;
marginSz = 0.04;
rasterSz = 0.04;
spacerSz = 0.01;
markerSz = 1;
fontSz   = 5 * (5 - nData);
tFact    = 0.1;
lnWidth  = 1.0;
axWidth  = (1 - (nData + 1) * marginSz) / nData;
axSz     = (1 - 2*marginSz - 2*rasterSz - 5*spacerSz)/2;

% Plot colors
% lfpCol     = [0.2 0.2 0.5];
% cellCol    = [0.5 0.2 0.2];
% swrCol     = [0.5 0.5 0.8];
% swrSpkCCol = [0.6 0.3 0.3];
% swrBstCCol = [0.3 0.0 0.0];
% spkCol     = [0.8 0.5 0.5];
% spkCCol    = [0.0 0.0 0.3];
% bstCol     = [0.8 0.5 0.5];
% bstCCol    = [0.0 0.0 0.3];

lfpCol{1}  = [ 48  70 160]/255;
lfpCol{2}  = [ 50  50  50]/255;
cellCol    = [  0  90   0]/255;
swrCol     = [180 180 180]/255;
swrSpkCCol = cellCol;
swrBstCCol = cellCol/2;
spkCol     = swrCol;
spkCCol{1} = lfpCol{1};
spkCCol{2} = lfpCol{2};
bstCol     = swrCol;
bstCCol{1} = lfpCol{1};
bstCCol{2} = lfpCol{2};

% Initialize graphical structures
hand.axTrLFP   = gobjects(nData, 1);
hand.axRsLFP   = gobjects(nData, 1);
hand.axTrCell  = gobjects(nData, 1);
hand.axRsSpk   = gobjects(nData, 1);
hand.axRsBst   = gobjects(nData, 1);
hand.lblLFP    = gobjects(nData, 1);
hand.lblCell   = gobjects(nData, 1);
hand.scaleLFP  = struct;
hand.scaleCell = struct;

% Initialize plot data
for i = 1:nData
  maxIndex     = length(data(i).SWR.spike.timingA);
  timing{i}    = downsampleMean(data(i).SWR.spike.timingA/1000, dsPlot);
  trLFP{i}     = downsampleMean(convFact * data(i).LFP.tSeries(1:maxIndex), dsPlot);
  trCell{i}    = downsampleMean(data(i).C.tSeries(1:maxIndex), dsPlot);
  rsSWR{i}     = downsampleMax(data(i).SWR.spike.evStatusA, dsPlot);
  rsSWRSpkC{i} = downsampleMax(data(i).SWR.spike.evStatusC, dsPlot);
  rsSpk{i}     = downsampleMax(data(i).C.spike.evStatusA, dsPlot);
  rsSpkC{i}    = downsampleMax(data(i).C.spike.evStatusC, dsPlot);
  
  if param.procBstOption
    rsSWRBstC{i} = downsampleMax(data(i).SWR.burst.evStatusC, dsPlot);
    rsBst{i}     = downsampleMax(data(i).C.burst.evStatusA, dsPlot);
    rsBstC{i}    = downsampleMax(data(i).C.burst.evStatusC, dsPlot);
  end
end

%% In reverse y order, initialize trace/raster plots
yPos = marginSz;

% In forward x order, initialize burst raster(s)
xPos = marginSz;
for i = 1:nData
  hand.axRsBst(i) = subplot('Position',[xPos yPos axWidth spacerSz]);
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + spacerSz + spacerSz;

% In forward x order, initialize spike raster(s)
xPos = marginSz;
for i = 1:nData
  hand.axRsSpk(i) = subplot('Position',[xPos yPos axWidth rasterSz]);
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
  
  % Assign coincident SWR points (with spikes)
  if (sum(rsSWRSpkC{i}) == 0)
    swrSpkCPoints  = timing{i};
    swrSpkCChannel = NaN * ones(1,20);
  else
    swrSpkCPoints  = timing{i}(find(rsSWRSpkC{i} .* timing{i}));
    swrSpkCChannel = (1:20);
  end
  swrSpkCChannel = swrSpkCChannel(ones(1,size(swrSpkCPoints,1)),:);
  plot(hand.axRsLFP(i), swrSpkCPoints, swrSpkCChannel,'.','MarkerSize', markerSz, 'Color', swrSpkCCol);
  
  % Assign coincident SWR points (with bursts - if selected)
  if param.procBstOption
    if (sum(rsSWRBstC{i}) == 0)
      swrBstCPoints  = timing{i};
      swrBstCChannel = NaN * ones(1,20);
    else
      swrBstCPoints  = timing{i}(find(rsSWRBstC{i} .* timing{i}));
      swrBstCChannel = (1:20);
    end
    swrBstCChannel = swrBstCChannel(ones(1,size(swrBstCPoints,1)),:);
    plot(hand.axRsLFP(i), swrBstCPoints, swrBstCChannel,'.','MarkerSize', markerSz, 'Color', swrBstCCol);
  end
  
  axis(hand.axRsLFP(i), [timing{i}(1) timing{i}(length(timing{i})) 1 20]);
  axis(hand.axRsLFP(i), 'off');
  hold(hand.axRsLFP(i), 'off');
  
end

%% Spike Plot
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
  hand.lblCell(i)= text(hand.axTrCell(i), hand.axTrCell(i).XLim(1), hand.axTrCell(i).YLim(2) - tFact * (hand.axTrCell(i).YLim(2) - hand.axTrCell(i).YLim(1)), 'Post-Synaptic Currents', 'FontSize', fontSz, 'Color', [0 0 0]);
  hand.scaleCell.(['scale' int2str(i)]) = createScaleBar(hand.axTrCell(i), [], 1, 1, fontSz, 's', 'pA');
end

%% Spike raster plot
for i = 1:nData
  
  hold(hand.axRsSpk(i), 'on');
  
  % Assign Spike points
  if (sum(rsSpk{i}) == 0)
    spkPoints  = timing{i};
    spkChannel = NaN * ones(1,20);
  else
    spkPoints  = timing{i}(find(rsSpk{i} .* timing{i}));
    spkChannel = (1:20);
  end
  spkChannel = spkChannel(ones(1,length(spkPoints)),:);
  plot(hand.axRsSpk(i), spkPoints, spkChannel, '.', 'MarkerSize', markerSz, 'Color', spkCol);
  
  % Assign coincident Spike points
  if (sum(rsSpkC{i}) == 0)
    spkCPoints  = timing{i};
    spkCChannel = NaN * ones(1,20);
  else
    spkCPoints  = timing{i}(find(rsSpkC{i} .* timing{i}));
    spkCChannel = (1:20);
  end
  spkCChannel = spkCChannel(ones(1,size(spkCPoints,1)),:);
  plot(hand.axRsSpk(i), spkCPoints, spkCChannel,'.','MarkerSize', markerSz, 'Color', spkCCol{i});
  
  axis(hand.axRsSpk(i), [timing{i}(1) timing{i}(length(timing{i})) 1 20]);
  axis(hand.axRsSpk(i), 'off');
  hold(hand.axRsSpk(i), 'off');

end

%% Burst raster plot (if selected)
if param.procBstOption
  for i = 1:nData
    
    hold(hand.axRsBst(i), 'on');
    
    % Assign Spike points
    if (sum(rsBst{i}) == 0)
      bstPoints  = timing{i};
      bstChannel = NaN * ones(1,20);
    else
      bstPoints  = timing{i}(find(rsBst{i} .* timing{i}));
      bstChannel = (1:20);
    end
    bstChannel = bstChannel(ones(1,length(bstPoints)),:);
    plot(hand.axRsBst(i), bstPoints, bstChannel, '.', 'MarkerSize', markerSz, 'Color', bstCol);
    
    % Assign coicident Burst points
    if (sum(rsBstC{i}) == 0)
      bstCPoints  = timing{i};
      bstCChannel = NaN * ones(1,20);
    else
      bstCPoints  = timing{i}(find(rsBstC{i} .* timing{i}));
      bstCChannel = (1:20);
    end
    bstCChannel = bstCChannel(ones(1,size(bstCPoints,1)),:);
    plot(hand.axRsBst(i), bstCPoints, bstCChannel,'.','MarkerSize', markerSz, 'Color', bstCCol{i});
    
    axis(hand.axRsBst(i), [timing{i}(1) timing{i}(length(timing{i})) 1 20]);
    axis(hand.axRsBst(i), 'off');
    hold(hand.axRsBst(i), 'off');
    
  end
end

%% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
set(pan(hand.spkFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});
set(zoom(hand.spkFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});

for i = 1:nData
  linkaxes([hand.axTrLFP(i), hand.axTrCell(i), hand.axRsLFP(i), hand.axRsSpk(i), hand.axRsBst(i)], 'x');
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
