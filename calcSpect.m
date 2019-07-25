function [S, BL] = calcSpect(S, BL, fRange, Fs, timeWindow, saveEvSpect)

% Handle input arguments - if not entered
if (nargin < 6) saveEvSpect = 0;     end % If true and > 1 tSeries will save all event spectra (Caution - large files)
if (nargin < 5) timeWindow  = 30;    end % ms
if (nargin < 4) Fs          = 20000; end % Hz
if (nargin < 3) fRange      = 1:500; end % Hz

% Reset spect structure
S.spect = struct;

% Initialize temporary tSeries cell array structure (to handle both full tSeries and series of events)
if isfield(S,'event') && ~isfield(S,'tSeries')
  tSeries = S.event;
elseif isfield(S,'tSeries')
  tSeries{1} = S.tSeries;
else
  error("no relevant tSeries found for spectogram")
end

% Set spectrogram parameters
nWin = round(3 * timeWindow * Fs / 1000);
nOv  = round(0.8 * nWin);
spectSamples = max(cellfun(@length, tSeries));

% Calculate baseline spectrogram (if supplied)
if ~isempty(BL)
  if ~isfield(BL,'spect') BL.spect = struct; end
  
  [~, BL.spect.fRange, BL.spect.tRange, BL.spect.power] = spectrogram(BL.tSeries, nWin, nOv, fRange, Fs, 'yaxis');
  
  pAve = mean(BL.spect.power, 2);
  pStd = std(BL.spect.power, 0, 2);
  pAve = pAve(:,ones(1, length(BL.spect.tRange)));
  pStd = pStd(:,ones(1, length(BL.spect.tRange)));
  
  BL.spect.pZScore = (BL.spect.power - pAve) ./ pStd;
  
  BL.spect = orderStruct(BL.spect);
end

% Count number of valid tSeries to determine size for array initialization:
fTemp = [];
tTemp = [];
ev    =  0;
for i = 1:length(tSeries)
  if ~isempty(tSeries{i}) && (length(tSeries{i}) == spectSamples)
    if isempty(fTemp) % only calculate fTemp and tTemp once at this stage
      [~, fTemp, tTemp, ~] = spectrogram(tSeries{i}, nWin, nOv, fRange, Fs, 'yaxis');
    end
    ev = ev + 1;
  end
end
power = zeros(length(fTemp), length(tTemp), ev);

% Spectral Analysis of valid tSeries:
ev = 1;
for i = 1:length(tSeries)
  if ~isempty(tSeries{i}) && (length(tSeries{i}) == spectSamples)
    [~, S.spect.fRange, S.spect.tRange, power(:,:,ev)] = spectrogram(tSeries{i}, nWin, nOv, fRange, Fs, 'yaxis');
    
    % Normalize to spectrogram unless baseline supplied:
    if isempty(BL)
      pAve = mean(power, 2);
      pStd = std(power, 0, 2);
      pAve = pAve(:,ones(1, length(S.spect.tRange)));
      pStd = pStd(:,ones(1, length(S.spect.tRange)));
    else
      pAve = mean(BL.spect.power, 2);
      pStd = std(BL.spect.power, 0, 2);
      pAve = pAve(:,ones(1, length(S.spect.tRange)));
      pStd = pStd(:,ones(1, length(S.spect.tRange)));
    end
    
    % Calculate Z-score for each tSeries
    pZScore = (power - pAve) ./ pStd;
    
    ev = ev + 1;
  end
end

% Calculate and save average spectrogram if more than one tSeries:
if length(tSeries) > 1
  
  % Save all event spectra if option selected:
  if saveEvSpect
    S.spect.power   = power;
    S.spect.pZScore = pZScore;
  end
  
  S.spect.powerAve = mean(power, 3);
  
  % Normalize to average spectrogram unless baseline supplied:
  if isempty(BL)
    pAve = mean(S.spect.powerAve, 2);
    pStd = std(S.spect.powerAve, 0, 2);
    pAve = pAve(:,ones(1, length(S.spect.tRange)));
    pStd = pStd(:,ones(1, length(S.spect.tRange)));
  else
    pAve = mean(BL.spect.power, 2);
    pStd = std(BL.spect.power, 0, 2);
    pAve = pAve(:,ones(1, length(S.spect.tRange)));
    pStd = pStd(:,ones(1, length(S.spect.tRange)));
  end
  
  S.spect.pZScoreAve = (S.spect.powerAve - pAve) ./ pStd;
  
else % only one tSeries - save spectra
  S.spect.power   = power;
  S.spect.pZScore = pZScore;
end

S.spect = orderStruct(S.spect);

end