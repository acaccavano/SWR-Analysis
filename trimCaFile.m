function tSeriesOut = trimCaFile(tSeriesIn, samplingIntCa, samplingIntLFP, stimStart, param)

%% Handle input arguments
if (nargin < 5) param          = struct; end
if (nargin < 4) stimStart      = []; end
if (nargin < 3) samplingIntLFP = []; end
if (nargin < 2) samplingIntCa  = []; end
if (nargin < 1) tSeriesIn      = []; end

% Handle case in which empty variable is supplied:
if isempty(param)          param          = struct; end
if isempty(samplingIntLFP) samplingIntLFP = 0.05; end
if isempty(samplingIntCa)  samplingIntCa  = 0.5; end

% Set default parameters if not specified
if ~isfield(param,'skipDetectLim') param.skipDetectLim = 4;    end %  s
if ~isfield(param,'stimCaOption')  param.stimCaOption  = 1;    end
if ~isfield(param,'stimCaLim1')    param.stimCaLim1    = 0;    end % ms
if ~isfield(param,'stimCaLim2')    param.stimCaLim2    = 1000; end % ms

% Set spontaneous peak detection mask
detectMask = logical(ones(size(tSeriesIn, 1), 1));

% Limit initial detection based on param.skipDetectLim
hiWin = round(1000 * param.skipDetectLim / samplingIntCa);
if (hiWin > 1) detectMask(1 : hiWin) = 0; end

% Limit detection range to periods outside of stim window (if option selected)
if param.stimCaOption
  for ev = 1:length(stimStart)
    loWin = max(round(stimStart(ev) * samplingIntLFP / samplingIntCa) + round(param.stimCaLim1 / samplingIntCa), 1);
    hiWin = min(round(stimStart(ev) * samplingIntLFP / samplingIntCa) + round(param.stimCaLim2 / samplingIntCa), size(tSeriesIn, 1));
    detectMask(loWin : hiWin) = 0;
  end
end

CaRange = 1:size(tSeriesIn, 1);
CaRange = CaRange';
CaRange = CaRange(detectMask);
tSeriesOut = tSeriesIn(CaRange, :);

end

