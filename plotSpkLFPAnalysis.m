function hand = plotSpkLFPAnalysis(data, hand, param, dsPlot)
%% hand = plotSpkLFPAnalysis(data, hand, param, dsPlot)
% 
%  Function to plot output of analyzeSpkFile for Spk-LFP phase analysis

nData = length(data);

% Input parameters
if isempty(param); param = struct; end
if ~isfield(param,'thetaOption');      param.thetaOption       = 1;    end
if ~isfield(param,'alphaOption');      param.thetaOption       = 0;    end
if ~isfield(param,'betaOption');       param.betaOption        = 0;    end
if ~isfield(param,'gammaOption');      param.gammaOption       = 1;    end
if ~isfield(param,'hgammaOption');     param.hgammaOption      = 1;    end

% Initialization
convFact = 1000; % Convert from mV to uV
nTrace   = 1;
colOption  = false; % If true will plot all traces of one data structure the same below defined colors, otherwise uses default ColorOrder
% dataCol{1} = [48 70 160]/255;
% dataCol{2} = [50 50  50]/255;

% Timing and Cell arrays - always necessary
for i = 1:nData
  timingPlot{i}          = downsampleMean(data(i).C.timing/1000, dsPlot);
%   peakTime{i}            = downsampleMean(data(i).C.timing(data.C.spike.evPeak)/1000, dsPlot);
  peakTime{i}            = data(i).C.timing(data.C.spike.evPeak)/1000;
  dataPlot{i, nTrace}    = downsampleMean(data(i).C.tSeries, dsPlot);
  dataPhase{i, nTrace}   = NaN * ones(length(data(i).C.timing),1); % Empty Placeholder
  dataName{i, nTrace}    = 'Cell-Attached';
end

% Theta plots
if isfield(data(1).C.spike, 'theta') && param.thetaOption
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).theta.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).theta.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Theta (' int2str(data(i).theta.lim1) '-' int2str(data(i).theta.lim2) 'Hz)'];
  end
end

% Alpha plots
if isfield(data(1).C.spike, 'alpha') && param.alphaOption
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).alpha.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).alpha.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Alpha (' int2str(data(i).alpha.lim1) '-' int2str(data(i).alpha.lim2) 'Hz)'];
  end
end

% Beta plots
if isfield(data(1).C.spike, 'beta') && param.betaOption
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).beta.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).beta.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Beta (' int2str(data(i).beta.lim1) '-' int2str(data(i).beta.lim2) 'Hz)'];
  end
end

% Gamma plots
if isfield(data(1).C.spike, 'gamma') && param.gammaOption
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).gamma.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).gamma.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['Gamma (' int2str(data(i).gamma.lim1) '-' int2str(data(i).gamma.lim2) 'Hz)'];
  end
end

% High Gamma plots
if isfield(data(1).C.spike, 'hgamma') && param.hgammaOption
  nTrace  = nTrace + 1;
  for i = 1:nData
    dataPlot{i, nTrace}  = downsampleMean(convFact * data(i).hgamma.tSeries, dsPlot);
    dataPhase{i, nTrace} = downsampleMean(data(i).hgamma.phase.tPhase, dsPlot);
    dataName{i, nTrace}  = ['High Gamma (' int2str(data(i).hgamma.lim1) '-' int2str(data(i).hgamma.lim2) 'Hz)'];
  end
end

% Plotting parameters
marginSz  = 0.04;
spacerSz  = 0.00;
fontSz    = 16;
tFact     = 0.1;
lnWidth   = 1.5;
yPos      = marginSz;
colInd    = nTrace;
colMatrix = colororder;
plotWidth = (1 - (nData + 1) * marginSz) / nData;
plotSz    = (1 - marginSz - spacerSz - yPos - nTrace*spacerSz)/nTrace;

% Plotting in order from lowest y to highest
hand.axTr  = gobjects(nData, nTrace);
hand.lblTr = gobjects(1, nTrace);
hand.scale = struct;

%% In reverse y order, initialize and plot traces:
for tr = nTrace : -1 : 1
  
  % In forward x order, plot LFP traces:
  xPos = marginSz;
  minY =   999999;
  maxY =  -999999;
  
  for i = 1:nData
    
    hand.axTr(i, tr) = subplot('Position',[xPos yPos plotWidth plotSz]);
    hold(hand.axTr(i, tr), 'on');
    
    if colOption
      traceColor = dataCol{i};
    else
      traceColor = colMatrix(colInd, :);
    end
    
    if tr > 1
      yyaxis(hand.axTr(i, tr), 'left')
      hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPlot{i, tr}, 'LineWidth', lnWidth, 'Color', traceColor);
      yyaxis(hand.axTr(i, tr), 'right')
      hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPhase{i, tr}, 'LineWidth', 0.5*lnWidth, 'Color', traceColor);
      ylim(hand.axTr(i, tr), [0 20*pi])
      yyaxis(hand.axTr(i, tr), 'left')
    else
      hand.plot = plot(hand.axTr(i, tr), timingPlot{i}, dataPlot{i, tr}, 'LineWidth', lnWidth, 'Color', traceColor);
    end
    
%     Plot viertical lines for spiked commented out - too slow   
%     for spk = 1:length(peakTime{i})
%       hand.plot = xline(hand.axTr(i, tr), peakTime{i}(spk), 'LineStyle', ':', 'LineWidth', 0.5*lnWidth, 'Color', colMatrix(1, :));
%     end

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
  if tr == 1
    hand.scale.(['scale' int2str(tr)]) = createScaleBar(hand.axTr(nData, tr), [], 1, 1, fontSz, 's', 'pA');
  else
    hand.scale.(['scale' int2str(tr)]) = createScaleBar(hand.axTr(nData, tr), [], 0, 1, fontSz);
  end
  
  % Update running y position and color index
  yPos = yPos + plotSz + spacerSz;
  colInd = colInd - 1;
end


%% Limit zoom and pan and set callbacks for scale bars
pan('xon');
zoom('xon');
hand.scale = orderfields(hand.scale);
set(pan(hand.spkFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact, nData});
set(zoom(hand.spkFig), 'ActionPostCallback', {@redrawPZCallback, hand, fontSz, tFact, nData});

% Link data axes
for i = 1:nData
  linkaxes([hand.axTr(i,:)], 'x');
end
% linkaxes(hand.axTr(:,:), 'y'); % Comment out - usually don't want this behavior (helpful for making figures)

end


% Callback function to redraw scalebar and labels after zooming/panning
function redrawPZCallback(obj, evd, hand, fontSz, tFact, nData)
scaleNames = fieldnames(hand.scale);
for i = 1:length(hand.lblTr)
  hand.lblTr(i).Position = [hand.axTr(1, i).XLim(1), hand.axTr(1, i).YLim(2) - tFact * (hand.axTr(1, i).YLim(2) - hand.axTr(1, i).YLim(1))];
  
  if i == 1
    createScaleBar(hand.axTr(nData, i), hand.scale.(scaleNames{i}), 1, 1, fontSz, 's', 'pA');
  else
    createScaleBar(hand.axTr(nData, i), hand.scale.(scaleNames{i}), 0, 1, fontSz);
  end
end

end