function [S, BL] = calcSpect(S, BL, fRange, Fs)

% Handle input arguments - if not entered
if (nargin < 4) Fs  = 3000; end
if (nargin < 3) fRange  = 1:250; end

timeWindow = 50; % in ms

% Set spectrogram parameters
nWin = round(3 * timeWindow * Fs / 1000);
nOv  = round(0.8 * nWin);

if ~isfield(S,'spect') S.spect = struct; end

if isempty(BL) % No baseline, normalize to S
  
  [~, S.spect.fRange, S.spect.tRange, S.spect.power] = spectrogram(S.tSeries, nWin, nOv, fRange, Fs, 'yaxis');
  pAve = mean(S.spect.power')';
  pStd = std(S.spect.power')';
  pAve = pAve(:,ones(1, length(S.spect.tRange)));
  pStd = pStd(:,ones(1, length(S.spect.tRange)));
  S.spect.pZScore = (S.spect.power - pAve) ./ pStd;
  
else % Use BL to normalize S z-score
  
  if ~isfield(BL,'spect') BL.spect = struct; end

  [~, BL.spect.fRange, BL.spect.tRange, BL.spect.power] = spectrogram(BL.tSeries, nWin, nOv, fRange, Fs, 'yaxis');
  pAve = mean(BL.spect.power')';
  pStd = std(BL.spect.power')';
  pAve = pAve(:,ones(1, length(BL.spect.tRange)));
  pStd = pStd(:,ones(1, length(BL.spect.tRange)));
  BL.spect.pZScore = (BL.spect.power - pAve) ./ pStd;
  
  BL.spect = orderStruct(BL.spect);
  
  [~, S.spect.fRange, S.spect.tRange, S.spect.power] = spectrogram(S.tSeries, nWin, nOv, fRange, Fs, 'yaxis');
  pAve = mean(BL.spect.power')';
  pStd = std(BL.spect.power')';
  pAve = pAve(:,ones(1, length(S.spect.tRange)));
  pStd = pStd(:,ones(1, length(S.spect.tRange)));
  S.spect.pZScore = (S.spect.power - pAve) ./ pStd;
  
end

S.spect = orderStruct(S.spect);

end