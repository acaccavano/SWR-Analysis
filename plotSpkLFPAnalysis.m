function hand = plotSpkLFPAnalysis(data, hand, param, dsPlot)
%% hand = plotSpkLFPAnalysis(data, hand, param, dsPlot)
% 
%  Function to plot output of analyzeSpkFile for Spk-LFP phase analysis

nData = length(data);

% Initialization
convFact = 1000; % Convert from mV to uV
nTrace   = 1;
nRaster  = 1;
colOption  = false; % If true will plot all traces of one data structure the same below defined colors, otherwise uses default ColorOrder
dataCol{1} = [48 70 160]/255;
dataCol{2} = [50 50  50]/255;

% Timing and Cell arrays - always necessary
for i = 1:nData
  timingPlot{i}          = downsampleMean(data(i).C.timing/1000, dsPlot);
  dataPlot{i, nTrace}    = downsampleMean(convFact * data(i).C.tSeries, dsPlot);
  dataPhase{i, nTrace}   = NaN * ones(length(data(i).C.timing),1); % Empty Placeholder
  dataName{i, nTrace}    = 'Cell-Attached';
  dataRaster{i, nRaster} = downsampleMax(data(i).C.spike.evStatus, dsPlot);
end

% Theta plots
if isfield(data(1).C.spike, 'theta')
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).theta.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).theta.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Theta (' int2str(data(i).theta.lim1) '-' int2str(data(i).theta.lim2) 'Hz)'];
  end
end

% Beta plots
if isfield(data(1).C.spike, 'beta')
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).beta.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).beta.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Beta (' int2str(data(i).beta.lim1) '-' int2str(data(i).beta.lim2) 'Hz)'];
  end
end

% Gamma plots
if isfield(data(1).C.spike, 'gamma')
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).gamma.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).gamma.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Gamma (' int2str(data(i).gamma.lim1) '-' int2str(data(i).gamma.lim2) 'Hz)'];
  end
end

% High Gamma plots
if isfield(data(1).C.spike, 'hgamma')
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).hgamma.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).hgamma.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['High Gamma (' int2str(data(i).hgamma.lim1) '-' int2str(data(i).hgamma.lim2) 'Hz)'];
  end
end

% Plotting parameters
marginSz  = 0.04;
rasterSz  = 0.01;
spacerSz  = 0.01;
fontSz    = 16;
tFact     = 0.1;
lnWidth   = 1.5;
yPos      = marginSz;
colInd    = nTrace;
colMatrix = colororder;
plotWidth = (1 - (nData + 1) * marginSz) / nData;
plotSz    = (1 - marginSz - spacerSz - yPos - nRaster*rasterSz - nTrace*spacerSz)/nTrace;

% Plotting in order from lowest y to highest
hand.axTr  = gobjects(nData, nTrace);
hand.axRs  = gobjects(nData, nRaster);
hand.lblTr = gobjects(1, nTrace);
hand.scale = struct;

%% In reverse y order, initialize and plot filtered LFP traces:
if nTrace > 1
  for tr = nTrace : -1 : 2
    
    % In forward x order, plot LFP traces:
    xPos = marginSz;
    minY =   999999;
    maxY =  -999999;
    
    for i = 1:nData
      
      hand.axTr(i, tr) = subplot('Position',[xPos yPos plotWidth plotSz]);
      hold(hand.axTr(i, tr), 'on');
      yyaxis(hand.axTr(i, tr), 'left')
      
      if colOption
        traceColor = dataCol{i};
      else
        traceColor = colMatrix(colInd, :);
      end
      
      hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPlot{i, tr}, 'LineWidth', lnWidth, 'Color', traceColor);
      yyaxis(hand.axTr(i, tr), 'right')
      hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPhase{i, tr}, 'LineWidth', 0.5*lnWidth, 'Color', traceColor);
      
      ylim(hand.axTr(i, tr), [0 20*pi])
      yyaxis(hand.axTr(i, tr), 'left')
      hold(hand.axTr(i, tr), 'off');
      
      minY = min(minY, min(dataPlot{i, tr}));
      maxY = max(maxY, max(dataPlot{i, tr}));
      xPos = xPos + plotWidth + marginSz;
    end
    
    % Set uniform axis ranges:
    for i = 1:nData
      axis(hand.axTr(i, tr), [timingPlot{i}(1) timingPlot{i}(length(timingPlot{i})) minY maxY]);
      axis(hand.axTr(i, tr), 'off');
    end
    
    hand.lblTr(tr) = text(hand.axTr(1, tr), hand.axTr(1, tr).XLim(1), hand.axTr(1, tr).YLim(2) - tFact * (hand.axTr(1, tr).YLim(2) - hand.axTr(1, tr).YLim(1)), dataName{1, tr}, 'FontSize', fontSz);
    
    % Scale bar:
    hand.scale.(['scale' int2str(tr)]) = createScaleBar(hand.axTr(nData, tr), [], 0, 1, fontSz);
    
    % Update running y position and color index
    yPos = yPos + plotSz + spacerSz;
    colInd = colInd - 1;
  end
end


%% Initialize and plot cell-attached tSeries and raster data:
tr = 1;
rs = 1;

% In forward x order, plot raster data:
if strcmp(dataName{1, tr}(1:4),'Cell')
  xPos = marginSz;
  for i = 1:nData
    
    hand.axRs(i, rs) = subplot('Position',[xPos yPos plotWidth rasterSz]);
    hand.axRs(i, rs).ColorOrderIndex = colInd;
    if colOption
      rasterCol = dataCol{i};
    else
      rasterCol = hand.axRs(i, rs).ColorOrder(hand.axRs(i, rs).ColorOrderIndex,:);
    end
    
    % Assign points to raster
    if (sum(dataRaster{i, rs}) == 0)
      rsPoints  = timingPlot{i};
      rsChannel = NaN * ones(1,20);
    else
      rsPoints  = timingPlot{i}(find(dataRaster{i, rs} .* timingPlot{i}));
      rsChannel = (1:20);
    end
    rsChannel = rsChannel(ones(1,length(rsPoints)),:);
    plot(hand.axRs(i, rs), rsPoints, rsChannel, '.', 'MarkerSize', 1, 'Color', rasterCol);
    axis(hand.axRs(i, rs), [timingPlot{i}(1) timingPlot{i}(length(timingPlot{i})) -inf inf]);
    axis(hand.axRs(i, rs), 'off');
    
    xPos = xPos + plotWidth + marginSz;
  end
  
  % Update running y position
  yPos = yPos + rasterSz + spacerSz;
end

%% In forward x order, plot traces:
xPos = marginSz;
minY =   999999;
maxY =  -999999;

for i = 1:nData
  
  hand.axTr(i, tr) = subplot('Position',[xPos yPos plotWidth plotSz]);
  hold(hand.axTr(i, tr), 'on');
  hand.axTr(i, tr).ColorOrderIndex = colInd;
  if colOption
    hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPlot{i, tr}, 'LineWidth', lnWidth, 'Color', dataCol{i});
  else
    hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPlot{i, tr}, 'LineWidth', lnWidth);
  end
  hold(hand.axTr(i, tr), 'off');
  
  minY = min(minY, min(dataPlot{i, tr}));
  maxY = max(maxY, max(dataPlot{i, tr}));
  xPos = xPos + plotWidth + marginSz;
end

% Set uniform axis ranges:
for i = 1:nData
  axis(hand.axTr(i, tr), [timingPlot{i}(1) timingPlot{i}(length(timingPlot{i})) minY maxY]);
  axis(hand.axTr(i, tr), 'off');
end

hand.lblTr(tr) = text(hand.axTr(1, tr), hand.axTr(1, tr).XLim(1), hand.axTr(1, tr).YLim(2) - tFact * (hand.axTr(1, tr).YLim(2) - hand.axTr(1, tr).YLim(1)), dataName{1, tr}, 'FontSize', fontSz);

% Scale bar with x-scale
hand.scale.(['scale' int2str(tr)]) = createScaleBar(hand.axTr(nData, tr), [], 1, 1, fontSz);


%% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
hand.scale = orderfields(hand.scale);
set(pan(hand.spkFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact, nData});
set(zoom(hand.spkFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact, nData});

% Link data axes
for i = 1:nData
  linkaxes([hand.axTr(i,:), hand.axRs(i,:)], 'x');
end
% linkaxes(hand.axTr(:,:), 'y'); % Comment out - usually don't want this behavior (helpful for making figures)

end


% Callback function to redraw scalebar and labels after zooming/panning
function redrawPZCallback(obj, evd, hand, fontSz, tFact, nData)
scaleNames = fieldnames(hand.scale);
for i = 1:length(hand.lblTr)
  hand.lblTr(i).Position = [hand.axTr(1, i).XLim(1), hand.axTr(1, i).YLim(2) - tFact * (hand.axTr(1, i).YLim(2) - hand.axTr(1, i).YLim(1))];
  
  if i == 1
    createScaleBar(hand.axTr(nData, i), hand.scale.(scaleNames{i}), 1, 1, fontSz);
  else
    createScaleBar(hand.axTr(nData, i), hand.scale.(scaleNames{i}), 0, 1, fontSz);
  end
end

end