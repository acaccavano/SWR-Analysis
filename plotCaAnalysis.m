function hand = plotCaAnalysis(data, hand, param, dsPlot)
%% hand = plotCaAnalysis(data, hand, param, dsPlot)
%  Function to plot output of analyzeCaFile

if (nargin < 4); dsPlot    = 1; end
if (nargin < 3); param     = struct; end
if (nargin < 2); hand      = struct; end
if (nargin < 1); data      = struct; end

% Set default parameters if not specified
if ~isfield(param,'skipDetectLim');  param.skipDetectLim  =  2; end
if ~isfield(param,'swrCaOption');    param.swrCaOption    =  0; end
if ~isfield(param,'spkCaOption');    param.spkCaOption    =  0; end
if ~isfield(param,'threshOption');   param.threshOption   =  1; end
if ~isfield(param,'peakOption');     param.peakOption     =  0; end
if ~isfield(param,'stimCaOption');   param.swrCaOption    =  0; end
if ~isfield(param,'cellLimOption');  param.cellLimOption  =  0; end % Limits number of Ca cells plotted - useful for comparisons where equal number of traces preferable
if ~isfield(param,'cellLim');        param.cellLim        = 80; end 
if ~isfield(param,'cellTypeOption'); param.cellTypeOption =  1; end
if ~isfield(param,'nCellTypes');     param.nCellTypes     =  2; end
if ~isfield(param,'cellTypeName');   param.cellTypeName{param.nCellTypes} = []; end
if ~isfield(param,'alignEndOption'); param.alignEndOption =  0; end

nData = length(data);

% Plot dimension parameters
convFact     = 1000;
marginSz     = 0.04;
spacerSz     = 0.02;
fontSz       = 5 * (5 - nData);
tFact        = 0.01;
lnWidth      = 1;
markerSz     = 1;
colInd       = 1;
axWidth      = (1 - (nData + 1) * marginSz) / nData;
axCaRsSz     = 0.20;
axCaTrSz     = 1 - axCaRsSz - spacerSz - 2*marginSz;
if param.swrCaOption || param.stimCaOption
  axLFPTrSz  = 0.12;
  axLFPRsSz  = 0.01;
  axCaTrSz   = axCaTrSz - (axLFPTrSz + axLFPRsSz + 2*spacerSz);
end
offsetFactor = 1;

lfpCol{1}    = [  0   0   0]/255; % Baseline = black
% lfpCol{2}    = [ 66  93 154]/255; % Control  = blue
lfpCol{2}    = [128 128 128]/255; % Sema3a   = grey

cellCol{1,1} = lfpCol{1}; % Set same as LFP
cellCol{2,1} = lfpCol{2}; % Set same as LFP

swrCol{1}    = lfpCol{1}; % Non-coincident
swrCol{2}    = lfpCol{2}; % Non-coincident

swrCCol{1}   = lfpCol{1}; % Coincident
swrCCol{2}   = lfpCol{2}; % Coincident

% LFP colors for data sets 1-3:
% lfpCol{1}  = [ 48  70 160]/255;
% lfpCol{2}  = [ 50  50  50]/255;
% lfpCol{3}  = [160  70  48]/255;

% % SWR raster event colors:
% swrCol     = [180 180 180]/255; % Non-coincident
% swrCCol    = lfpCol;            % Coincident

% % For Cell Types 1-6:
% cellCol{1} = [ 68 114 196]/255; % 1 = blue
% cellCol{2} = [192   0   0]/255; % 2 = red
% cellCol{3} = [255 192   0]/255; % 3 = yellow
% cellCol{4} = [ 84 130  53]/255; % 4 = green
% cellCol{5} = [204   0 153]/255; % 5 = purple
% cellCol{6} = [237 125  49]/255; % 6 = orange

% colMap = lines;
% cellCol{1} = colMap(2,:); 
% cellCol{2} = colMap(3,:);
% cellCol{3} = colMap(4,:);
% cellCol{4} = colMap(5,:);
% cellCol{5} = colMap(6,:);
% cellCol{6} = colMap(7,:);

% Initialize graphical structures
hand.scale   = struct;
hand.axCaTr  = gobjects(nData, 1);
hand.axCaRs  = gobjects(nData, 1);
hand.lblCaTr = gobjects(nData, 1);
nChannelsMax = 1;

% Calcium timing arrays and dFoF data - always necessary
for i = 1:nData

  % Downsample if selected (otherwise pass-through) 
  timingCa{i} = downsampleMean(data(i).Ca.timing/1000, dsPlot);
  CaRange{i}  = find(data(i).Ca.timing >= 1000 * param.skipDetectLim);
  
  if param.cellLimOption
    nCells(i) = min(param.cellLim, data(i).Ca.nChannels);
    if param.cellTypeOption
      cellType(i,:) = data(i).Ca.cellType(1:param.cellLim);
      nCellTypes(i) = unique(cellType(i));
    end
  else
    nCells(i) = data(i).Ca.nChannels;
    if param.cellTypeOption
      cellType(i,:) = data(i).Ca.cellType;
      nCellTypes(i) = param.nCellTypes;
    end
  end
  
  for ch = 1:nCells(i)
    tSeriesCa{i}(:,ch) = downsampleMean(data(i).Ca.tSeries(:,ch), dsPlot);
    tSeriesEv{i}(:,ch) = data(i).Ca.tSeries(CaRange{i},ch);
    
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
      timingCaEv{i}     = data(i).Ca.timing(CaRange{i})/1000;
      timingCaRs{i}     = downsampleMean(timingCaEv{i}, dsPlot);
      rasterCa{i}(:,ch) = downsampleMax(data(i).Ca.evStatus(:,ch), dsPlot);
      if param.peakOption
        peakCaVal{i}{ch}  = tSeriesEv{i}(data(i).Ca.evPeak{ch}, ch);
        peakCaTime{i}{ch} = timingCaEv{i}(data(i).Ca.evPeak{ch});
      end
    end
  end
  % Find maximum number to plot on consistent axes
  nChannelsMax = max(nChannelsMax, nCells(i));
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
  colInd       = colInd + 1;
  
  for i = 1:nData
    
    if param.swrCaOption
      if param.alignEndOption
        timeArray = data(i).LFP.timing + (data(i).SWR.Ca.timingA(end) - data(i).LFP.timing(end));
      else
        timeArray = data(i).LFP.timing;
      end
      timingLFP{i}   = downsampleMean(timeArray/1000, dsPlot);
      tSeriesLFP{i}  = downsampleMean(convFact * data(i).LFP.tSeries, dsPlot);
      timingLFPRs{i} = downsampleMean(data(i).SWR.Ca.timingA/1000, dsPlot);
      rasterLFP{i}   = downsampleMax(data(i).SWR.Ca.evStatusA, dsPlot);
      rasterLFPC{i}  = downsampleMax(data(i).SWR.Ca.evStatusSumC, dsPlot);
      
    elseif param.stimCaOption
      if param.alignEndOption
        timeArray = data(i).LFP.timing + (data(i).stim.Ca.timingA(end) - data(i).LFP.timing(end));
      else
        timeArray = data(i).LFP.timing;
      end
      timingLFP{i}   = downsampleMean(timeArray/1000, dsPlot);
      tSeriesLFP{i}  = downsampleMean(convFact * data(i).LFP.tSeries, dsPlot);
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
  offset = offsetFactor * linspace(1, nCells(i), nCells(i));
  offsetArray = offset(ones(1,length(timingCa{i})), :);
  colInd = hand.axCaRs(i).ColorOrderIndex;

  % Plot dFoF as lines
  hand.axCaTr(i).ColorOrderIndex = colInd;
  if param.cellTypeOption
    for j = 1:nCellTypes(i)
      cellInd = (cellType(i,:) == j);
      plot(hand.axCaTr(i), timingCa{i}, offsetArray(:, cellInd) + tSeriesCa{i}(:, cellInd), 'LineWidth', lnWidth, 'Color', cellCol{i,j});
      axis(hand.axCaTr(i), [0 timingCa{i}(length(timingCa{i})) min(min(offsetArray + tSeriesCa{i})) max(max(offsetArray(:, cellInd) + tSeriesCa{i}(:, cellInd)))]);
      text(hand.axCaTr(i), hand.axCaTr(i).XLim(1), min(offsetArray(:, find(cellInd, 1, 'last')) + tSeriesCa{i}(:, find(cellInd, 1, 'last'))) - (fontSz - 5), param.cellTypeName{j}, 'FontSize', fontSz - 5, 'Color', [0 0 0]);
    end
  else
    plot(hand.axCaTr(i), timingCa{i}, offsetArray + tSeriesCa{i}, 'LineWidth', lnWidth, 'Color', cellCol{i,1});
  end
  
  % Add thresholds channel by channel
  if param.threshOption
    for ch = 1:nCells(i)
      offsetThresh = offsetArray(1,ch) + data(i).Ca.peakThresh(ch);
      if param.cellTypeOption
        plot(hand.axCaTr(i), [timingCaRs{i}(1) timingCaRs{i}(length(timingCaRs{i}))],[offsetThresh offsetThresh], 'LineWidth', 0.25, 'LineStyle', '--', 'Color', cellCol{i,cellType(i,ch)});
      else
        plot(hand.axCaTr(i), [timingCaRs{i}(1) timingCaRs{i}(length(timingCaRs{i}))],[offsetThresh offsetThresh], 'LineWidth', 0.25, 'LineStyle', '--', 'Color', cellCol{i,1});
      end
    end
    hand.axCaTr(i).ColorOrderIndex = colInd;
  end
  
  % Add peaks channel by channel, and plot raster data
  for ch = 1:nCells(i)
    
    % Assign Ca events to raster:
    if (sum(rasterCa{i}(:,ch)) == 0)
      plot(hand.axCaTr(i), timingCaRs{i}, NaN * ones(length(timingCaRs{i}), 1));
      evPoints  = timingCaRs{i};
      evChannel = NaN * ones(length(timingCaRs{i}), 1);
    else
      
      % Plot peaks:
      if param.peakOption
        offsetPeaks = peakCaVal{i}{ch} + offsetArray(1,ch) + 0.2;
        plot(hand.axCaTr(i), peakCaTime{i}{ch}, offsetPeaks, 'v', 'MarkerSize', markerSz - 1, 'Color', [1 1 1]);
      end
      
      evPoints  = timingCaRs{i}(find(rasterCa{i}(:,ch) .* timingCaRs{i}));
      evChannel = ch * ones(size(evPoints,1),1);
    end
    plot(hand.axCaRs(i), evPoints, evChannel, 's', 'MarkerSize', markerSz, 'MarkerEdgeColor', swrCol{i}, 'MarkerFaceColor', swrCol{i});
    
    % Assign Calcium events coicident with SWRs to an additional raster:
    if param.swrCaOption || param.stimCaOption
      
      if (sum(rasterCaC{i}(:,ch)) == 0)
        evPoints  = timingCaRs{i};
        evChannel = NaN * ones(length(timingCaRs{i}), 1);
      else
        evPoints  = timingCaRs{i}(find(rasterCaC{i}(:,ch) .* timingCaRs{i}));
        evChannel = ch * ones(size(evPoints,1),1);
      end
      
      if param.cellTypeOption
        plot(hand.axCaRs(i), evPoints, evChannel, 's', 'MarkerSize', markerSz, 'MarkerEdgeColor', cellCol{i,cellType(i,ch)}, 'MarkerFaceColor', cellCol{i,cellType(i,ch)});
      else
        plot(hand.axCaRs(i), evPoints, evChannel, 's', 'MarkerSize', markerSz, 'MarkerEdgeColor', cellCol{i,1}, 'MarkerFaceColor', cellCol{i,1});
      end
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
%   hand.lblCaTr(i) = text(hand.axCaTr(i), hand.axCaTr(i).XLim(1), hand.axCaTr(i).YLim(2) - tFact * (hand.axCaTr(i).YLim(2) - hand.axCaTr(i).YLim(1)), data(i).saveName, 'FontSize', fontSz, 'Interpreter', 'none');
  hand.scale.(['scale' int2str(i)]) = createScaleBar(hand.axCaTr(i), [], 1, 1, fontSz, 's', '\DeltaF/F', 10, 30);
end

%% Plot LFP trace and SWR raster (if selected)
if param.swrCaOption || param.stimCaOption
  for i = 1:nData
    
    % Assign SWR events to raster:
    if (sum(rasterLFP{i}) == 0)
      evPoints  = timingLFPRs{i};
      evChannel = NaN * ones(1,20);
    else
      evPoints  = timingLFPRs{i}(find(rasterLFP{i} .* timingLFPRs{i}));
      evChannel = (1:20);
    end
    evChannel = evChannel(ones(1,length(evPoints)),:);
    plot(hand.axLFPRs(i), evPoints, evChannel, '.', 'MarkerSize', markerSz, 'Color', swrCol{i});
    
    % Assign SWR events coincident with at least one Ca transient to an additional raster:
    if (sum(rasterLFPC{i}) == 0)
      evPoints  = timingLFPRs{i};
      evChannel = NaN * ones(1,20);
    else
      evPoints  = timingLFPRs{i}(find(rasterLFPC{i} .* timingLFPRs{i}));
      evChannel = (1:20);
    end
    evChannel = evChannel(ones(1,length(evPoints)),:);
    plot(hand.axLFPRs(i), evPoints, evChannel, '.', 'MarkerSize', markerSz, 'Color', swrCCol{i});
    
    axis(hand.axLFPRs(i), [0 timingLFPRs{i}(length(timingLFPRs{i})) -inf inf]);
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
  createScaleBar(hand.axCaTr(i), hand.scale.(scaleNames{i}), 1, 1, fontSz, 's', '\DeltaF/F', 30, 50);
end
end

