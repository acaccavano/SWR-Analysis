function S = calcEvFFT(S, param, lim1, lim2)

if ~isfield(S,'FFT') S.FFT = struct; end

% (Re)Initialize data arrays:
S.FFT.F   = [];
S.FFT.F1  = [];
S.FFT.F2  = [];

S.FFT.F{length(S.event)}  = [];
S.FFT.F1{length(S.event)} = [];
S.FFT.F2{length(S.event)} = [];

S.FFT.F  = S.FFT.F';
S.FFT.F1 = S.FFT.F1';
S.FFT.F2 = S.FFT.F2';

S.FFT.pkFreq = NaN * zeros(length(S.event), 1);

% Compute FFT for valid events:
nSamples = max(cellfun(@length,S.event));
for ev = 1:length(S.event)
  if ~isempty(S.event{ev}) && (length(S.event{ev}) == nSamples)
    S.FFT.F{ev}  = fft(S.event{ev}, nSamples - 1);
    S.FFT.F2{ev} = abs(S.FFT.F{ev}/(nSamples - 1));
    S.FFT.F1{ev} = S.FFT.F2{ev}(1 : (nSamples - 1)/2 + 1);
    S.FFT.F1{ev}(2 : end - 1) = 2 * S.FFT.F1{ev}(2 : end - 1);
    S.FFT.fftRange = param.Fs * (0:((nSamples - 1)/2))/(nSamples - 1);
    S.FFT.subRange = find(S.FFT.fftRange>lim1 & S.FFT.fftRange<lim2);
    [~, pkFreqInd] = max(S.FFT.F1{ev}(S.FFT.subRange));
    S.FFT.pkFreq(ev) = S.FFT.fftRange(S.FFT.subRange(pkFreqInd));
  end
end

% Calculate average FFT and peak frequency:
S.FFT.fftAve    = mean(cat(2, S.FFT.F1{:}), 2);
S.FFT.fftAve    = S.FFT.fftAve';
[~, pkFreqInd]  = max(S.FFT.fftAve(S.FFT.subRange));
S.FFT.pkFreqAve = S.FFT.fftRange(S.FFT.subRange(pkFreqInd));
S.FFT = orderStruct(S.FFT);

end