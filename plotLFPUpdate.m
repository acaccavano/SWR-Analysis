function hand = plotLFPUpdate(data, data2, hand, param, ev, dsPlot)
%% Resize LFPPlot centered around specific SWR event
% Will downsample if selected, otherwise just a pass-through and name change
% If data2 entered will plot both side-by-side on same scale

tMin = max(round(data.LFP.timing(data.SWR.evPeak(ev)) - 2*param.swrWindow), 0)/1000;
tMax = min(round(data.LFP.timing(data.SWR.evPeak(ev)) + 2*param.swrWindow), max(data.LFP.timing))/1000;

hand.axTr(1).XLim = [tMin tMax];

end