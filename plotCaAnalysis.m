function hand = plotCaAnalysis(data, hand, param, dsPlot)

nData = length(data);

% Set default parameters if not specified
if ~isfield(param,'skipDetectLim')        param.skipDetectLim        = 2;   end
if ~isfield(param,'swrCaOption')          param.swrCaOption          = 0;   end
if ~isfield(param,'spkCaOption')          param.spkCaOption          = 0;   end

% Plot dimension parameters
convFact     = 1000;
marginSz     = 0.04;
spacerSz     = 0.02;
fontSz       = 5 * (5 - nData);
tFact        = 0.01;
lnWidth      = 1.5;
markerSz     = 3;
colInd       = 1;
axWidth      = (1 - (nData + 1) * marginSz) / nData;

% Initialize graphical structures
hand.scale   = struct;
hand.axCaTr  = gobjects(nData, 1);
hand.axCaRs  = gobjects(nData, 1);
hand.lblCaTr = gobjects(nData, 1);
nChannelsMax = 1;
axCaRsSz     = 0.22;
axCaTrSz     = 1 - axCaRsSz - spacerSz - 2*marginSz;

% Calcium timing arrays and dFoF data - always necessary
for i = 1:nData
  timingCa{i} = downsampleMean(data(i).Ca.timing/1000, dsPlot);
  CaRange{i}  = find(data(i).Ca.timing >= 1000 * param.skipDetectLim);
  
  for ch = 1:data(i).Ca.nChannels
    tSeriesCa{i}(:,ch)  = downsampleMean(data(i).Ca.tSeries(:,ch), dsPlot);
    tSeriesEv{i}(:,ch)  = data(i).Ca.tSeries(CaRange{i},ch);
    
    if param.swrCaOption
      timingCaEv{i}      = data(i).SWR.Ca.timingA/1000;
      timingCaRs{i}      = downsampleMean(timingCaEv{i}, dsPlot);
      rasterCa{i}(:,ch)  = downsampleMax(data(i).Ca.SWR.evStatusA(:,ch), dsPlot);
      rasterCaC{i}(:,ch) = downsampleMax(data(i).Ca.SWR.evStatusC(:,ch), dsPlot);
    else
      timingCaEv{i}      = data(i).Ca.timing(CaRange{i})/1000;
      timingCaRs{i}      = downsampleMean(timingCaEv{i}, dsPlot);
      rasterCa{i}(:,ch)  = downsampleMax(data(i).Ca.evStatus(:,ch), dsPlot);
    end
  end
  % Find maximum number to plot on consistent axes
  nChannelsMax = max(nChannelsMax, data(i).Ca.nChannels);
end

% %% Initialize Spike traces/raster data (if selected) - DATA STRUCTURES DON'T EXIST YET
% if param.spkCaOption
%   hand.axSpkTr = gobjects(nData, 1);
%   hand.axSpkRs = gobjects(nData, 1);
%   axSpkTrSz    = 0.04;
%   axSpkRsSz    = 0.02;
%   colInd       = colInd + 1;
%   axCaTrSz     = axCaTrSz - (axSpkTrSz + axSpkRsSz + 2*spacerSz);
%   
%   for i = 1:nData
%     timingSpk{i} = downsampleMean(data(i).C.Ca.timingA/1000, dsPlot); 
%     [tSeriesSpk{i}, ~, ~, ~, ~, ~, ~] = timeAlign(data(i).C.tSeries, data(i).C.Ca.evStatusA, data(i).C.timing, data(i).C.Ca.timingA);
%     tSeriesSpk{i} = downsampleMean(tSeriesSpk{i}, dsPlot); 
%     rasterSpk{i}  = downsampleMax(data(i).C.Ca.evStatusA, dsPlot);
%   end
% end

%% Initialize SWR traces/raster data (if selected)
if param.swrCaOption
  hand.axLFPTr = gobjects(nData, 1);
  hand.axSWRRs = gobjects(nData, 1);
  axLFPTrSz    = 0.04;
  axSWRRsSz    = 0.02;
  colInd       = colInd + 1;
  axCaTrSz     = axCaTrSz - (axLFPTrSz + axSWRRsSz + 2*spacerSz);

  for i = 1:nData
    timingLFP{i}   = downsampleMean(data(i).LFP.timing/1000, dsPlot);
    tSeriesLFP{i}  = downsampleMean(convFact * data(i).LFP.tSeries, dsPlot);
    timingLFPRs{i} = downsampleMean(data(i).SWR.Ca.timingA/1000, dsPlot);
    rasterSWR{i}   = downsampleMax(data(i).SWR.Ca.evStatusA, dsPlot);
    rasterSWRC{i}  = downsampleMax(data(i).SWR.Ca.evStatusSumC, dsPlot);
  end
end

%% In reverse y order, initialize trace and raster plots
yPos = marginSz;

% In forward x order, initialize Ca raster plot
xPos = marginSz;
for i = 1:nData
  hand.axCaRs(i) = subplot('Position',[xPos yPos axWidth axCaRsSz]);
  axis(hand.axCaRs(i), 'off');
  hold(hand.axCaRs(i), 'on');
  hand.axCaRs(i).ColorOrderIndex = colInd;
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + axCaRsSz + spacerSz;

% In forward x order, initialize Ca trace plot
xPos = marginSz;
for i = 1:nData
  hand.axCaTr(i) = subplot('Position',[xPos yPos axWidth axCaTrSz]);
  axis(hand.axCaTr(i), 'off');
  hold(hand.axCaTr(i), 'on');
  hand.axCaTr(i).ColorOrderIndex = colInd;
  xPos = xPos + axWidth + spacerSz;
end
yPos = yPos + axCaTrSz + spacerSz;

% % If selected, initialize in forward x-order Spike trace and raster
% if param.spkCaOption
%   colInd = colInd - 1;
%   xPos = marginSz;
%   for i = 1:nData
%     hand.axSpkRs(i) = subplot('Position',[xPos yPos axWidth axSpkRsSz]);
%     axis(hand.axSpkRs(i), 'off');
%     hold(hand.axSpkRs(i), 'on');
%     hand.axSpkRs(i).ColorOrderIndex = colInd;
%     xPos = xPos + axWidth + spacerSz;
%   end
%   yPos = yPos + axSpkRsSz + spacerSz;
%   
%   xPos = marginSz;
%   for i = 1:nData
%     hand.axSpkTr(i) = subplot('Position',[xPos yPos axWidth axSpkTrSz]);
%     axis(hand.axSpkTr(i), 'off');
%     hold(hand.axSpkTr(i), 'on');
%     hand.axSpkTr(i).ColorOrderIndex = colInd;
%     xPos = xPos + axWidth + spacerSz;
%   end
%   yPos = yPos + axCaTrSz + spacerSz;
% end

% If selected, initialize in forward x-order SWR trace and raster
if param.swrCaOption
  colInd = colInd - 1;
  xPos = marginSz;
  for i = 1:nData
    hand.axSWRRs(i) = subplot('Position',[xPos yPos axWidth axSWRRsSz]);
    axis(hand.axSWRRs(i), 'off');
    hold(hand.axSWRRs(i), 'on');
    hand.axSWRRs(i).ColorOrderIndex = colInd;
    xPos = xPos + axWidth + spacerSz;
  end
  yPos = yPos + axSWRRsSz + spacerSz;
  
  xPos = marginSz;
  for i = 1:nData
    hand.axLFPTr(i) = subplot('Position',[xPos yPos axWidth axLFPTrSz]);
    axis(hand.axLFPTr(i), 'off');
    hold(hand.axLFPTr(i), 'on');
    hand.axLFPTr(i).ColorOrderIndex = colInd;
    xPos = xPos + axWidth + spacerSz;
  end
  yPos = yPos + axLFPTrSz + spacerSz;
end


%% Plot Ca traces/raster
minY =  999999;
maxY = -999999;
for i = 1:nData
  
  % Create offset arrays to distribute traces
  offset = linspace(1, data(i).Ca.nChannels, data(i).Ca.nChannels);
  offsetArray = offset(ones(1,length(timingCa{i})), :);
  colInd = hand.axCaRs(i).ColorOrderIndex;

  % Plot dFoF as lines
  plot(hand.axCaTr(i), timingCa{i}, offsetArray + tSeriesCa{i}, 'LineWidth', lnWidth);
  hand.axCaTr(i).ColorOrderIndex = colInd;
  
  % Add thresholds channel by channel
  for ch = 1:data(i).Ca.nChannels
    offsetThresh = offsetArray(1,ch) + data(i).Ca.peakThresh(ch);
    plot(hand.axCaTr(i), [timingCaRs{i}(1) timingCaRs{i}(length(timingCaRs{i}))],[offsetThresh offsetThresh], 'LineWidth', 0.25, 'LineStyle', '--');
  end
  hand.axCaTr(i).ColorOrderIndex = colInd;
  
  % Add peaks channel by channel, and plot raster data
  for ch = 1:data(i).Ca.nChannels
    rasterCol = hand.axCaRs(i).ColorOrder(hand.axCaRs(i).ColorOrderIndex,:);
    
    % Assign Ca events to raster:
    if (sum(rasterCa{i}(:,ch)) == 0)
      plot(hand.axCaTr(i), timingCaRs{i}, NaN * ones(length(timingCaRs{i}), 1));
      evPoints  = timingCaRs{i};
      evChannel = NaN * ones(length(timingCaRs{i}), 1);
    else
      offsetPeaks = tSeriesEv{i}(data(i).Ca.evPeak{ch}, ch) + offsetArray(1,ch) + 0.2;
      plot(hand.axCaTr(i), timingCaEv{i}(data(i).Ca.evPeak{ch}), offsetPeaks, 'v', 'MarkerSize', markerSz);
      evPoints  = timingCaRs{i}(find(rasterCa{i}(:,ch) .* timingCaRs{i}));
      evChannel = ch * ones(size(evPoints,1),1);
    end
    plot(hand.axCaRs(i), evPoints, evChannel, 's', 'MarkerSize', markerSz, 'MarkerEdgeColor', rasterCol, 'MarkerFaceColor', rasterCol);
    
    % Assign Calcium events coicident with SWRs to an additional raster:
    if param.swrCaOption
      rasterCol = [0 0 0];
      if (sum(rasterCaC{i}(:,ch)) == 0)
        evPoints  = timingCaRs{i};
        evChannel = NaN * ones(length(timingCaRs{i}), 1);
      else
        evPoints  = timingCaRs{i}(find(rasterCaC{i}(:,ch) .* timingCaRs{i}));
        evChannel = ch * ones(size(evPoints,1),1);
      end
      plot(hand.axCaRs(i), evPoints, evChannel, 's', 'MarkerSize', markerSz, 'MarkerEdgeColor', rasterCol, 'MarkerFaceColor', rasterCol);
    end
    
  end
  
  minY = min(minY, min(min(offsetArray + tSeriesCa{i})));
  maxY = max(maxY, max(max(offsetArray + tSeriesCa{i})));
  
  hold(hand.axCaTr(i), 'off');
  hold(hand.axCaRs(i), 'off');
end

% Set uniform axis ranges:
for i = 1:nData
  axis(hand.axCaTr(i), [timingCa{i}(1) timingCa{i}(length(timingCa{i})) minY maxY]);
  axis(hand.axCaRs(i), [timingCaRs{i}(1) timingCaRs{i}(length(timingCaRs{i})) 0 nChannelsMax]);
  hand.lblCaTr(i) = text(hand.axCaTr(i), hand.axCaTr(i).XLim(1), hand.axCaTr(i).YLim(2) - tFact * (hand.axCaTr(i).YLim(2) - hand.axCaTr(i).YLim(1)), data(i).saveName, 'FontSize', fontSz, 'Interpreter', 'none');
  hand.scale.(['scale' int2str(i)]) = createScaleBar(hand.axCaTr(i), [], 1, 1, fontSz, 's', '\DeltaF/F');
end

%% Plot LFP trace and SWR raster (if selected)
if param.swrCaOption
  for i = 1:nData

    % Assign SWR events to raster:
    rasterCol = hand.axSWRRs(i).ColorOrder(hand.axSWRRs(i).ColorOrderIndex,:);
    if (sum(rasterSWR{i}) == 0)
      evPoints  = timingLFPRs{i};
      evChannel = NaN * ones(1,20);
    else
      evPoints  = timingLFPRs{i}(find(rasterSWR{i} .* timingLFPRs{i}));
      evChannel = (1:20);
    end
    evChannel = evChannel(ones(1,length(evPoints)),:);
    plot(hand.axSWRRs(i), evPoints, evChannel, '.', 'MarkerSize', 1, 'Color', rasterCol);
    
    % Assign SWR events coincident with at least one Ca transient to an additional raster:
    rasterCol = [0 0 0];
    if (sum(rasterSWRC{i}) == 0)
      evPoints  = timingLFPRs{i};
      evChannel = NaN * ones(1,20);
    else
      evPoints  = timingLFPRs{i}(find(rasterSWRC{i} .* timingLFPRs{i}));
      evChannel = (1:20);
    end
    evChannel = evChannel(ones(1,length(evPoints)),:);
    plot(hand.axSWRRs(i), evPoints, evChannel, '.', 'MarkerSize', 1, 'Color', rasterCol);
    
    axis(hand.axSWRRs(i), [timingLFPRs{i}(1) timingLFPRs{i}(length(timingLFPRs{i})) -inf inf]);
  end
  
  % Plot LFP trace:
  minY =  999999;
  maxY = -999999;
  for i = 1:nData
    plot(hand.axLFPTr(i), timingLFP{i}, tSeriesLFP{i}, 'LineWidth', lnWidth);
    minY = min(minY, min(tSeriesLFP{i}));
    maxY = max(maxY, max(tSeriesLFP{i}));
  end
  
  % Set uniform axis ranges:
  for i = 1:nData
    axis(hand.axLFPTr(i), [timingLFP{i}(1) timingLFP{i}(length(timingLFP{i})) minY maxY]);
  end
  
end

% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
set(pan(hand.CaFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});
set(zoom(hand.CaFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});

for i = 1:nData
  if param.swrCaOption
    linkaxes([hand.axCaTr(i), hand.axCaRs(i), hand.axLFPTr(i), hand.axSWRRs(i)], 'x');
  else
    linkaxes([hand.axCaTr(i), hand.axCaRs(i)], 'x');
  end
end
% 
% linkaxes([hand.axCaTr(1), hand.axCaTr(2)], 'y');
% linkaxes([hand.axCaRs(1), hand.axCaRs(2)], 'y');

end


% Callback function to redraw scalebar and labels after zooming/panning
function redrawPZCallback(obj, evd, hand, fontSz, tFact)
scaleNames = fieldnames(hand.scale);
for i = 1:length(hand.lblCaTr)
  hand.lblCaTr(i).Position = [hand.axCaTr(i).XLim(1), hand.axCaTr(i).YLim(2) - tFact * (hand.axCaTr(i).YLim(2) - hand.axCaTr(i).YLim(1))];
  createScaleBar(hand.axCaTr(i), hand.scale.(scaleNames{i}), 1, 1, fontSz, 's', '\DeltaF/F');
end
end

