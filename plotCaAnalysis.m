function hand = plotCaAnalysis(data, hand, param, dsPlot)

nData = length(data);

% Set default parameters if not specified
if ~isfield(param,'skipDetectLim') param.skipDetectLim = 2; end
if ~isfield(param,'swrCaOption')   param.swrCaOption   = 0; end
if ~isfield(param,'spkCaOption')   param.spkCaOption   = 0; end
if ~isfield(param,'colOption')     param.colOption     = 1; end % If enabled plots keeps all traces same colors, reserves below defined colors for different datasets
if ~isfield(param,'threshOption')  param.threshOption  = 1; end
if ~isfield(param,'peakOption')    param.peakOption    = 1; end
if ~isfield(param,'stimCaOption')  param.swrCaOption   = 0; end

% Plot dimension parameters
convFact     = 1000;
marginSz     = 0.04;
spacerSz     = 0.02;
fontSz       = 5 * (5 - nData);
tFact        = 0.01;
lnWidth      = 1;
markerSz     = 3;
colInd       = 1;
axWidth      = (1 - (nData + 1) * marginSz) / nData;

% Plot colors (for colOption = 1)
lfpCol{1}  = [ 48  70 160]/255;
lfpCol{2}  = [ 50  50  50]/255;
cellCol    = [  0  90   0]/255;
swrCol     = [180 180 180]/255;
swrCCol    = cellCol;
CaCol      = swrCol;
CaCCol{1}  = lfpCol{1};
CaCCol{2}  = lfpCol{2};

% Initialize graphical structures
hand.scale   = struct;
hand.axCaTr  = gobjects(nData, 1);
hand.axCaRs  = gobjects(nData, 1);
hand.lblCaTr = gobjects(nData, 1);
nChannelsMax = 1;
axCaRsSz     = 0.20;
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
      if param.peakOption
        peakCaVal{i}{ch}  = tSeriesEv{i}(data(i).Ca.SWR.evPeakA{ch}, ch);
        peakCaTime{i}{ch} = timingCaEv{i}(data(i).Ca.SWR.evPeakA{ch});
      end
      
    elseif param.stimCaOption
      timingCaEv{i}      = data(i).stim.Ca.timingA/1000;
      timingCaRs{i}      = downsampleMean(timingCaEv{i}, dsPlot);
      rasterCa{i}(:,ch)  = downsampleMax(data(i).Ca.stim.evStatusA(:,ch), dsPlot);
      rasterCaC{i}(:,ch) = downsampleMax(data(i).Ca.stim.evStatusC(:,ch), dsPlot);
      if param.peakOption
        peakCaVal{i}{ch}  = tSeriesEv{i}(data(i).Ca.stim.evPeakA{ch}, ch);
        peakCaTime{i}{ch} = timingCaEv{i}(data(i).Ca.stim.evPeakA{ch});
      end
      
    else
      timingCaEv{i}      = data(i).Ca.timing(CaRange{i})/1000;
      timingCaRs{i}      = downsampleMean(timingCaEv{i}, dsPlot);
      rasterCa{i}(:,ch)  = downsampleMax(data(i).Ca.evStatus(:,ch), dsPlot);
      if param.peakOption
        peakCaVal{i}{ch}  = tSeriesEv{i}(data(i).Ca.evPeak{ch}, ch);
        peakCaTime{i}{ch} = timingCaEv{i}(data(i).Ca.evPeak{ch});
      end
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
if param.swrCaOption || param.stimCaOption
  hand.axLFPTr = gobjects(nData, 1);
  hand.axLFPRs = gobjects(nData, 1);
  axLFPTrSz    = 0.22;
  axLFPRsSz    = 0.02;
  colInd       = colInd + 1;
  axCaTrSz     = axCaTrSz - (axLFPTrSz + axLFPRsSz + 2*spacerSz);
  
  for i = 1:nData
    timingLFP{i}   = downsampleMean(data(i).LFP.timing/1000, dsPlot);
    tSeriesLFP{i}  = downsampleMean(convFact * data(i).LFP.tSeries, dsPlot);
    
    if param.swrCaOption
      timingLFPRs{i} = downsampleMean(data(i).SWR.Ca.timingA/1000, dsPlot);
      rasterLFP{i}   = downsampleMax(data(i).SWR.Ca.evStatusA, dsPlot);
      rasterLFPC{i}  = downsampleMax(data(i).SWR.Ca.evStatusSumC, dsPlot);
      
    elseif param.stimCaOption
      timingLFPRs{i} = downsampleMean(data(i).stim.Ca.timingA/1000, dsPlot);
      rasterLFP{i}   = downsampleMax(data(i).stim.Ca.evStatusA, dsPlot);
      rasterLFPC{i}  = downsampleMax(data(i).stim.Ca.evStatusSumC, dsPlot);
      
    end
    
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

% If selected, initialize in forward x-order LFP trace and raster
if param.swrCaOption || param.stimCaOption
  colInd = colInd - 1;
  xPos = marginSz;
  for i = 1:nData
    hand.axLFPRs(i) = subplot('Position',[xPos yPos axWidth axLFPRsSz]);
    axis(hand.axLFPRs(i), 'off');
    hold(hand.axLFPRs(i), 'on');
    hand.axLFPRs(i).ColorOrderIndex = colInd;
    xPos = xPos + axWidth + spacerSz;
  end
  yPos = yPos + axLFPRsSz + spacerSz;
  
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
  hand.axCaTr(i).ColorOrderIndex = colInd;
  if param.colOption
    plot(hand.axCaTr(i), timingCa{i}, offsetArray + tSeriesCa{i}, 'LineWidth', lnWidth, 'Color', cellCol);
  else
    plot(hand.axCaTr(i), timingCa{i}, offsetArray + tSeriesCa{i}, 'LineWidth', lnWidth);
  end
  
  % Add thresholds channel by channel
  if param.threshOption
    for ch = 1:data(i).Ca.nChannels
      offsetThresh = offsetArray(1,ch) + data(i).Ca.peakThresh(ch);
      if param.colOption
        plot(hand.axCaTr(i), [timingCaRs{i}(1) timingCaRs{i}(length(timingCaRs{i}))],[offsetThresh offsetThresh], 'LineWidth', 0.25, 'LineStyle', '--', 'Color', cellCol);
      else
        plot(hand.axCaTr(i), [timingCaRs{i}(1) timingCaRs{i}(length(timingCaRs{i}))],[offsetThresh offsetThresh], 'LineWidth', 0.25, 'LineStyle', '--');
      end
    end
    hand.axCaTr(i).ColorOrderIndex = colInd;
  end
  
  % Add peaks channel by channel, and plot raster data
  for ch = 1:data(i).Ca.nChannels
    
    if param.colOption
      rasterCol = CaCol;
    else
      rasterCol = hand.axCaRs(i).ColorOrder(hand.axCaRs(i).ColorOrderIndex,:);
    end
    
    % Assign Ca events to raster:
    if (sum(rasterCa{i}(:,ch)) == 0)
      
      plot(hand.axCaTr(i), timingCaRs{i}, NaN * ones(length(timingCaRs{i}), 1));
      evPoints  = timingCaRs{i};
      evChannel = NaN * ones(length(timingCaRs{i}), 1);
      
    else
      
      % Plot peaks:
      if param.peakOption
        offsetPeaks = peakCaVal{i}{ch} + offsetArray(1,ch) + 0.2;
        if param.colOption
          plot(hand.axCaTr(i), peakCaTime{i}{ch}, offsetPeaks, 'v', 'MarkerSize', markerSz - 1, 'Color', cellCol);
        else
          plot(hand.axCaTr(i), peakCaTime{i}{ch}, offsetPeaks, 'v', 'MarkerSize', markerSz - 1);
        end
      end
      
      evPoints  = timingCaRs{i}(find(rasterCa{i}(:,ch) .* timingCaRs{i}));
      evChannel = ch * ones(size(evPoints,1),1);
    end
    plot(hand.axCaRs(i), evPoints, evChannel, 's', 'MarkerSize', markerSz, 'MarkerEdgeColor', rasterCol, 'MarkerFaceColor', rasterCol);
    
    % Assign Calcium events coicident with SWRs to an additional raster:
    if param.swrCaOption || param.stimCaOption
      
      if param.colOption
        rasterCol = CaCCol{i};
      else
        rasterCol = [0 0 0];
      end
      
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
  axis(hand.axCaTr(i), [0 timingCa{i}(length(timingCa{i})) minY maxY]); % timingCa{i}(1)
  axis(hand.axCaRs(i), [0 timingCaRs{i}(length(timingCaRs{i})) 0 nChannelsMax]); % timingCaRs{i}(1)
  hand.lblCaTr(i) = text(hand.axCaTr(i), hand.axCaTr(i).XLim(1), hand.axCaTr(i).YLim(2) - tFact * (hand.axCaTr(i).YLim(2) - hand.axCaTr(i).YLim(1)), data(i).saveName, 'FontSize', fontSz, 'Interpreter', 'none');
  hand.scale.(['scale' int2str(i)]) = createScaleBar(hand.axCaTr(i), [], 1, 1, fontSz, 's', '\DeltaF/F');
end

%% Plot LFP trace and SWR raster (if selected)
if param.swrCaOption || param.stimCaOption
  for i = 1:nData
    
    % Assign SWR events to raster:
    if param.colOption
      rasterCol = swrCol;
    else
      rasterCol = hand.axLFPRs(i).ColorOrder(hand.axLFPRs(i).ColorOrderIndex,:);
    end

    if (sum(rasterLFP{i}) == 0)
      evPoints  = timingLFPRs{i};
      evChannel = NaN * ones(1,20);
    else
      evPoints  = timingLFPRs{i}(find(rasterLFP{i} .* timingLFPRs{i}));
      evChannel = (1:20);
    end
    evChannel = evChannel(ones(1,length(evPoints)),:);
    plot(hand.axLFPRs(i), evPoints, evChannel, '.', 'MarkerSize', markerSz, 'Color', rasterCol);
    
    % Assign SWR events coincident with at least one Ca transient to an additional raster:
    if param.colOption
      rasterCol = swrCCol;
    else
      rasterCol = [0 0 0];
    end

    if (sum(rasterLFPC{i}) == 0)
      evPoints  = timingLFPRs{i};
      evChannel = NaN * ones(1,20);
    else
      evPoints  = timingLFPRs{i}(find(rasterLFPC{i} .* timingLFPRs{i}));
      evChannel = (1:20);
    end
    evChannel = evChannel(ones(1,length(evPoints)),:);
    plot(hand.axLFPRs(i), evPoints, evChannel, '.', 'MarkerSize', markerSz, 'Color', rasterCol);
    
    axis(hand.axLFPRs(i), [0 timingLFPRs{i}(length(timingLFPRs{i})) -inf inf]); % timingLFPRs{i}(1)
  end
  
  % Plot LFP trace:
  minY =  999999;
  maxY = -999999;
  for i = 1:nData
    if param.swrCaOption
      plot(hand.axLFPTr(i), timingLFP{i}, tSeriesLFP{i}, 'LineWidth', lnWidth, 'Color', lfpCol{i});
    else
      plot(hand.axLFPTr(i), timingLFP{i}, tSeriesLFP{i}, 'LineWidth', lnWidth);
    end
    minY = min(minY, min(tSeriesLFP{i}));
    maxY = max(maxY, max(tSeriesLFP{i}));
  end
  
  % Set uniform axis ranges:
  for i = 1:nData
    axis(hand.axLFPTr(i), [0 timingLFP{i}(length(timingLFP{i})) minY maxY]); % timingLFP{i}(1)
  end
  
end

% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
set(pan(hand.CaFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});
set(zoom(hand.CaFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact});

for i = 1:nData
  if param.swrCaOption || param.stimCaOption
    linkaxes([hand.axCaTr(i), hand.axCaRs(i), hand.axLFPTr(i), hand.axLFPRs(i)], 'x');
  else
    linkaxes([hand.axCaTr(i), hand.axCaRs(i)], 'x');
  end
end

% Usually comment below out, only for making figures
% linkaxes([hand.axCaTr, hand.axCaRs, hand.axLFPTr, hand.axLFPRs], 'x');
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

