function hand = plotLFPAnalysis(data, hand, param, dsPlot)
%% Initialize plotting data based on options entered in param structure
% Will downsample if selected, otherwise just a pass-through and name change
% If data2 and/or data3 entered will plot side-by-side on same scale

% Initialization 
convFact = 1000; % Convert from mV to uV
nData    = length(data);
nTrace   = 1;
nRaster  = 0;
limSpectCol = true;
maxPZScore = 10;

% Timing and LFP array - always necessary
for i = 1:nData
  timingPlot{i} = downsampleMean(data(i).LFP.timing/1000, dsPlot);
  dataPlot{i, nTrace} = downsampleMean(convFact * data(i).LFP.tSeries, dsPlot);
  dataName{i, nTrace} = 'LFP';
  if isfield(data(i).LFP, 'lim1')
    dataName{i, nTrace} = [dataName{i, nTrace} ' (' int2str(data(i).LFP.lim1) '-' int2str(data(i).LFP.lim2) 'Hz)'];
  end
end

% SWR raster events
if param.swrOption
  nRaster = 1;
  for i = 1:nData
    dataRaster{i, nRaster} = downsampleMax(data(i).SWR.evStatus, dsPlot);
  end
end

% SW plots and raster events
if param.swOption
  nTrace  = nTrace + 1;
  nRaster = nRaster + 1;
  for i = 1:nData
    dataPlot{i, nTrace}    = downsampleMean(convFact * data(i).SW.tSeries, dsPlot);
    dataRaster{i, nRaster} = downsampleMax(data(i).SW.evStatus, dsPlot);
    dataName{i, nTrace}    = ['SW (' int2str(data(i).SW.lim1) '-' int2str(data(i).SW.lim2) 'Hz)'];
  end
  
  if param.rmsOption
    nTrace = nTrace + 1;
    for i = 1:nData
      dataPlot{i, nTrace} = downsampleMean(convFact * data(i).SW.RMS, dsPlot);
      dataName{i, nTrace} = 'SW RMS';
    end
  end
end

% Ripple plots and raster events
if param.rOption
  nTrace  = nTrace + 1;
  nRaster = nRaster + 1;
  for i = 1:nData
    dataPlot{i, nTrace}    = downsampleMean(convFact * data(i).R.tSeries, dsPlot);
    dataRaster{i, nRaster} = downsampleMax(data(i).R.evStatus, dsPlot);
    dataName{i, nTrace}    = ['R (' int2str(data(i).R.lim1) '-' int2str(data(i).R.lim2) 'Hz)'];
  end
  
  if param.rmsOption
    nTrace = nTrace + 1;
    for i = 1:nData
      dataPlot{i, nTrace} = downsampleMean(convFact * data(i).R.RMS, dsPlot);
      dataName{i, nTrace} = 'R RMS';
    end
  end
end

% Theta plot
if param.thetaOption
  nTrace = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace} = downsampleMean(convFact * data(i).theta.tSeries, dsPlot);
    dataName{i, nTrace} = ['Theta (' int2str(data(i).theta.lim1) '-' int2str(data(i).theta.lim2) 'Hz)'];
  end
end

% Beta plot
if param.betaOption
  nTrace = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace} = downsampleMean(convFact * data(i).beta.tSeries, dsPlot);
    dataName{i, nTrace} = ['Beta (' int2str(data(i).beta.lim1) '-' int2str(data(i).beta.lim2) 'Hz)'];
  end
end

% Gamma plot
if param.gammaOption
  nTrace = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace} = downsampleMean(convFact * data(i).gamma.tSeries, dsPlot);
    dataName{i, nTrace} = ['Low Gamma (' int2str(data(i).gamma.lim1) '-' int2str(data(i).gamma.lim2) 'Hz)'];
  end
end

% High gamma plot
if param.hgammaOption
  nTrace = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace} = downsampleMean(convFact * data(i).hgamma.tSeries, dsPlot);
    dataName{i, nTrace} = ['High Gamma (' int2str(data(i).hgamma.lim1) '-' int2str(data(i).hgamma.lim2) 'Hz)'];
  end
end

% Plotting parameters
marginSz  = 0.04;
rasterSz  = 0.01;
spacerSz  = 0.01;
spectSz   = 0.20;
fontSz    = 16;
tFact     = 0.1;
lnWidth   = 1.5;
plotWidth = (1 - (nData + 1) * marginSz) / nData;

%% Plotting in order from lowest y to highest
hand.axTr  = gobjects(nData, nTrace);
hand.axRs  = gobjects(nData, nRaster);
hand.axSp  = gobjects(nData, 1);
hand.lblTr = gobjects(1, nTrace);
hand.scale = struct;

xPos   = marginSz;
yPos   = marginSz;
colInd = nTrace;

%% Spectrogram
if param.spectOption
  
  minC =   999999;
  maxC =  -999999;
  
  for i = 1:nData
    if limSpectCol % Option to remove outliers to get better dynamic range
      data(i).LFP.spect.pZScore(data(i).LFP.spect.pZScore > maxPZScore) = maxPZScore;
    end
    
    hand.axSp(i) = subplot('Position',[xPos yPos plotWidth spectSz]);
    imagesc(hand.axSp(i),'XData', data(i).LFP.spect.tRange,'YData', param.spectLim1:param.spectLim2, 'CData', data(i).LFP.spect.pZScore(param.spectLim1:param.spectLim2,:));
    axis(hand.axSp(i), [data(i).LFP.spect.tRange(1) data(i).LFP.spect.tRange(length(data(i).LFP.spect.tRange)) param.spectLim1 param.spectLim2]);
    set(hand.axSp(i),'FontSize',fontSz);
    
    minC = min(minC, min(min(data(i).LFP.spect.pZScore)));
    maxC = max(maxC, max(max(data(i).LFP.spect.pZScore)));
    xPos = xPos + plotWidth + marginSz;
  end
  
  colormap jet
  hand.lblSp  = text(hand.axSp(1), hand.axSp(1).XLim(1), hand.axSp(1).YLim(2) - tFact * (hand.axSp(1).YLim(2) - hand.axSp(1).YLim(1)), 'Spectrogram Z-Score', 'FontSize', fontSz, 'Color', [1 1 1]);
  hand.colbar = colorbar(hand.axSp(nData), 'Position', [1.01-marginSz marginSz 0.01 spectSz]);
  
  % Adjust color scale for all spectrogams to same:
  for i = 1:nData
    caxis(hand.axSp(i), [minC maxC]);
  end
  
  yPos = yPos + spectSz + spacerSz;
end
plotSz = (1 - marginSz - spacerSz - yPos - nRaster*rasterSz - nTrace*spacerSz)/nTrace;

%% In reverse y order, initialize and plot trace and raster data:
rs = nRaster;
for tr = nTrace : -1 : 1
  
  %% In forward x order, plot raster data (if applicable):
  if strcmp(dataName{1, tr}(1:4),'SW (') || strcmp(dataName{1, tr}(1:3),'R (') || (strcmp(dataName{1, tr}(1:3),'LFP') && param.swrOption)
    xPos = marginSz;
    for i = 1:nData
      
      hand.axRs(i, rs) = subplot('Position',[xPos yPos plotWidth rasterSz]);
      hand.axRs(i, rs).ColorOrderIndex = colInd;
      rasterCol = hand.axRs(i, rs).ColorOrder(hand.axRs(i, rs).ColorOrderIndex,:);
      
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
    
    % Increment raster index and update running y position
    rs = rs - 1;
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
    hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPlot{i, tr}, 'LineWidth', lnWidth);
    
    % If RMS selected, plot thresholds for peak detection:
    if strcmp(dataName{i, tr},'SW RMS')
      plot(hand.axTr(i, tr), [data(i).LFP.timing(1) data(i).LFP.timing(length(data(i).LFP.timing))/1000], convFact * [data(i).SW.peakThresh data(i).SW.peakThresh], 'Color', hand.plot.Color, 'LineWidth', 1.0, 'LineStyle', '--');
      plot(hand.axTr(i, tr), [data(i).LFP.timing(1) data(i).LFP.timing(length(data(i).LFP.timing))/1000], convFact * [data(i).SW.baseThresh data(i).SW.baseThresh], 'Color', hand.plot.Color, 'LineWidth', 0.6, 'LineStyle', '--');
    elseif strcmp(dataName{tr},'R RMS')
      plot(hand.axTr(i, tr), [data(i).LFP.timing(1) data(i).LFP.timing(length(data(i).LFP.timing))/1000], convFact * [data(i).R.peakThresh data(i).R.peakThresh], 'Color', hand.plot.Color, 'LineWidth', 1.0, 'LineStyle', '--');
      plot(hand.axTr(i, tr), [data(i).LFP.timing(1) data(i).LFP.timing(length(data(i).LFP.timing))/1000], convFact * [data(i).R.baseThresh data(i).R.baseThresh], 'Color', hand.plot.Color, 'LineWidth', 0.6, 'LineStyle', '--');
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
  hand.scale.(['scale' int2str(tr)]) = createScaleBar(hand.axTr(nData, tr), [], 0, 1, fontSz);

  % Update running y position and color index
  yPos = yPos + plotSz + spacerSz;
  colInd = colInd - 1;
end

% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
hand.scale = orderfields(hand.scale);
set(pan(hand.lfpFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact, nData});
set(zoom(hand.lfpFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact, nData});

% Link data axes
for i = 1:nData
  if (nRaster > 0) && param.spectOption
    linkaxes([hand.axTr(i,:), hand.axRs(i,:), hand.axSp(i)'], 'x');
  elseif (nRaster > 0)
    linkaxes([hand.axTr(i,:), hand.axRs(i,:)], 'x');
  elseif param.spectOption
    linkaxes([hand.axTr(i,:), hand.axSp(i)'], 'x');
  else
    linkaxes(hand.axTr(i,:), 'x');
  end
end
linkaxes(hand.axTr(:,:), 'y'); % Comment out - usually don't want this behavior (helpful for making figures)

end

% Callback function to redraw scalebar and labels after zooming/panning
function redrawPZCallback(obj, evd, hand, fontSz, tFact, nData)
scaleNames = fieldnames(hand.scale);
for i = 1:length(hand.lblTr)
  hand.lblTr(i).Position = [hand.axTr(1, i).XLim(1), hand.axTr(1, i).YLim(2) - tFact * (hand.axTr(1, i).YLim(2) - hand.axTr(1, i).YLim(1))];
  createScaleBar(hand.axTr(nData, i), hand.scale.(scaleNames{i}), 0, 1, fontSz);
end

if isfield(hand, 'lblSp')
  hand.lblSp.Position = [hand.axSp(1).XLim(1), hand.axSp(1).YLim(2) - tFact * (hand.axSp(1).YLim(2) - hand.axSp(1).YLim(1))];
end

end