function hand = plotLFPUpdate(data, hand, param, ev)
%% hand = plotLFPUpdate(data, data2, hand, param, ev, dsPlot)
%
% Function to resize LFP plot centered around specific SWR event

tMin = max(round(data.LFP.timing(data.SWR.evPeak(ev)) - 2*param.swrWindow), 0)/1000;
tMax = min(round(data.LFP.timing(data.SWR.evPeak(ev)) + 2*param.swrWindow), max(data.LFP.timing))/1000;

hand.axTr(1).XLim = [tMin tMax];

end