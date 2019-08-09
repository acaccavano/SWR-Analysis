function [S, hand] = analyzePhase(T, param)

if (nargin < 2) param = struct; end

hand  = struct;

if ~isfield(param,'sz') param.sz = 20;  end
if ~isfield(param,'al') param.al = 0.8; end
if ~isfield(param,'lw') param.lw = 1.8; end

cMark_CT = [99  143 235]/256;
cMark_AD = [224 140  65]/256;
cLine_CT = [ 35  51 116]/256;
cLine_AD = [143  75   7]/256;

% Copy over imported table fields to new structure array:
S       = struct;
S.G     = struct;
S.R     = struct;
S.file  = T.file;
S.gen   = T.genotype;
S.sex   = T.sex;
S.G.Ph  = T.gPh;
S.G.Var = T.gVar;
S.G.Len = 1 - T.gVar;
S.R.Ph  = T.rPh;
S.R.Var = T.rVar;
S.R.Len = 1 - T.rVar;

% Calculate gamma stats
[S.G.avePhCT, S.G.avePhLoCT, S.G.avePhHiCT, S.G.aveVarCT, S.G.aveLenCT] = calcPhaseStats(S.G.Ph(S.gen=="CT"), S.G.Len(S.gen=="CT"));
[S.G.avePhAD, S.G.avePhLoAD, S.G.avePhHiAD, S.G.aveVarAD, S.G.aveLenAD] = calcPhaseStats(S.G.Ph(S.gen=="AD"), S.G.Len(S.gen=="AD"));

% Calculate ripple stats
[S.R.avePhCT, S.R.avePhLoCT, S.R.avePhHiCT, S.R.aveVarCT, S.R.aveLenCT] = calcPhaseStats(S.R.Ph(S.gen=="CT"), S.R.Len(S.gen=="CT"));
[S.R.avePhAD, S.R.avePhLoAD, S.R.avePhHiAD, S.R.aveVarAD, S.R.aveLenAD] = calcPhaseStats(S.R.Ph(S.gen=="AD"), S.R.Len(S.gen=="AD"));

%% Plot Gamma scatter
hand.gFig = figure('Name', 'Gamma Phase');
hand.gAx  = polaraxes(hand.gFig);

hand.gPl_CT_M = polarscatter(hand.gAx, S.G.Ph(S.gen=="CT" & S.sex=="M"), S.G.Len(S.gen=="CT" & S.sex=="M"), param.sz, cMark_CT, 'filled');
hand.gPl_CT_M = editPolarPlot(hand.gPl_CT_M, param);
hold on

hand.gPl_CT_F = polarscatter(hand.gAx, S.G.Ph(S.gen=="CT" & S.sex=="F"), S.G.Len(S.gen=="CT" & S.sex=="F"), param.sz, cMark_CT);
hand.gPl_CT_F = editPolarPlot(hand.gPl_CT_F, param);

hand.gPl_AD_M = polarscatter(hand.gAx, S.G.Ph(S.gen=="AD" & S.sex=="M"), S.G.Len(S.gen=="AD" & S.sex=="M"), param.sz, cMark_AD, 'filled');
hand.gPl_AD_M = editPolarPlot(hand.gPl_AD_M, param);

hand.gPl_AD_F = polarscatter(hand.gAx, S.G.Ph(S.gen=="AD" & S.sex=="F"), S.G.Len(S.gen=="AD" & S.sex=="F"), param.sz, cMark_AD);
hand.gPl_AD_F = editPolarPlot(hand.gPl_AD_F, param);

gAvePh_CT_Plot = S.G.avePhCT * ones(1,100);
gAvePhLL_CT_Plot = S.G.avePhLoCT * ones(1,100);
gAvePhUL_CT_Plot = S.G.avePhHiCT * ones(1,100);
gAveLen_CT_Plot = linspace(0, S.G.aveLenCT, 100);

gAvePh_AD_Plot = S.G.avePhAD * ones(1,100);
gAvePhLL_AD_Plot = S.G.avePhLoAD * ones(1,100);
gAvePhUL_AD_Plot = S.G.avePhHiAD * ones(1,100);
gAveLen_AD_Plot = linspace(0, S.G.aveLenAD, 100);

polarplot(hand.gAx, gAvePh_CT_Plot, gAveLen_CT_Plot, 'Color', cLine_CT, 'LineWidth', param.lw + 1);
polarplot(hand.gAx, gAvePhLL_CT_Plot, gAveLen_CT_Plot, '--', 'Color', cLine_CT, 'LineWidth', param.lw - 1);
polarplot(hand.gAx, gAvePhUL_CT_Plot, gAveLen_CT_Plot, '--', 'Color', cLine_CT, 'LineWidth', param.lw - 1);

polarplot(hand.gAx, gAvePh_AD_Plot, gAveLen_AD_Plot, 'Color', cLine_AD, 'LineWidth', param.lw + 1);
polarplot(hand.gAx, gAvePhLL_AD_Plot, gAveLen_AD_Plot, '--', 'Color', cLine_AD, 'LineWidth', param.lw - 1);
polarplot(hand.gAx, gAvePhUL_AD_Plot, gAveLen_AD_Plot, '--', 'Color', cLine_AD, 'LineWidth', param.lw - 1);

hand.gAx.ThetaTick = [0 90 180 270];
hand.gAx.RLim  = [0 1];


%% Plot Ripple scatter
hand.rFig = figure('Name', 'Ripple Phase');
hand.rAx  = polaraxes(hand.rFig);

hand.rPl_CT_M = polarscatter(hand.rAx, S.R.Ph(S.gen=="CT" & S.sex=="M"), S.R.Len(S.gen=="CT" & S.sex=="M"), param.sz, cMark_CT, 'filled');
hand.rPl_CT_M = editPolarPlot(hand.rPl_CT_M, param);
hold on

hand.rPl_CT_F = polarscatter(hand.rAx, S.R.Ph(S.gen=="CT" & S.sex=="F"), S.R.Len(S.gen=="CT" & S.sex=="F"), param.sz, cMark_CT);
hand.rPl_CT_F = editPolarPlot(hand.rPl_CT_F, param);

hand.rPl_AD_M = polarscatter(hand.rAx, S.R.Ph(S.gen=="AD" & S.sex=="M"), S.R.Len(S.gen=="AD" & S.sex=="M"), param.sz, cMark_AD, 'filled');
hand.rPl_AD_M = editPolarPlot(hand.rPl_AD_M, param);

hand.rPl_AD_F = polarscatter(hand.rAx, S.R.Ph(S.gen=="AD" & S.sex=="F"), S.R.Len(S.gen=="AD" & S.sex=="F"), param.sz, cMark_AD);
hand.rPl_AD_F = editPolarPlot(hand.rPl_AD_F, param);

rAvePh_CT_Plot = S.R.avePhCT * ones(1,100);
rAvePhLL_CT_Plot = S.R.avePhLoCT * ones(1,100);
rAvePhUL_CT_Plot = S.R.avePhHiCT * ones(1,100);
rAveLen_CT_Plot = linspace(0, S.R.aveLenCT, 100);

rAvePh_AD_Plot = S.R.avePhAD * ones(1,100);
rAvePhLL_AD_Plot = S.R.avePhLoAD * ones(1,100);
rAvePhUL_AD_Plot = S.R.avePhHiAD * ones(1,100);
rAveLen_AD_Plot = linspace(0, S.R.aveLenAD, 100);

polarplot(hand.rAx, rAvePh_CT_Plot, rAveLen_CT_Plot, 'Color', cLine_CT, 'LineWidth', param.lw + 1);
polarplot(hand.rAx, rAvePhLL_CT_Plot, rAveLen_CT_Plot, '--', 'Color', cLine_CT, 'LineWidth', param.lw - 1);
polarplot(hand.rAx, rAvePhUL_CT_Plot, rAveLen_CT_Plot, '--', 'Color', cLine_CT, 'LineWidth', param.lw - 1);

polarplot(hand.rAx, rAvePh_AD_Plot, rAveLen_AD_Plot, 'Color', cLine_AD, 'LineWidth', param.lw + 1);
polarplot(hand.rAx, rAvePhLL_AD_Plot, rAveLen_AD_Plot, '--', 'Color', cLine_AD, 'LineWidth', param.lw - 1);
polarplot(hand.rAx, rAvePhUL_AD_Plot, rAveLen_AD_Plot, '--', 'Color', cLine_AD, 'LineWidth', param.lw - 1);

hand.rAx.ThetaTick = [0 90 180 270];
hand.rAx.RLim  = [0 1];

end

% Commented out: don't use manual calc, instead use below
% function [avePh, avePhLL, avePhUL, aveVar, aveLen] = calcPhaseStats(ph, len)
% 
% phX = cos(ph);
% phY = sin(ph);
% 
% avePh = atan2(sum(phY.*len)/length(phX), sum(phX.*len)/length(phX));
% if (avePh < 0) avePh = avePh + 2*pi; end
% aveLen = sqrt(sum(phX.*len)*sum(phX.*len) + sum(phY.*len)*sum(phY.*len))/length(phX);
% aveVar = 1 - aveLen;
% 
% [~, sd] = circ_std(ph);
% 
% avePhLL = avePh - sd;
% avePhUL = avePh + sd;
% end

function [avePh, avePhLL, avePhUL, aveVar, aveLen] = calcPhaseStats(ph, len)

[avePh, avePhUL, avePhLL] = circ_mean(ph);
if (avePh < 0) 
  avePh   = avePh + 2*pi; 
  avePhUL = avePhUL + 2*pi; 
  avePhLL = avePhLL + 2*pi; 
end
aveLen = circ_r(ph);
aveVar = 1 - aveLen;

end


function h = editPolarPlot(h, param)

h.SizeData = param.sz;
h.MarkerFaceAlpha = param.al;
h.MarkerEdgeAlpha = param.al;
h.MarkerEdgeColor = 'flat';
h.LineWidth = param.lw;

end

