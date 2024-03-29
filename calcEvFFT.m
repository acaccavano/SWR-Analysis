function S = calcEvFFT(S, param, lim1, lim2)
%% S = calcEvFFT(S, param, lim1, lim2)
%
%  Function to create new structure and calculate fast-fourier transform
%  (FFT) for an array of events on event-by-event basis and for the average
%  of all events. Unlike calcTotFFT, no windowing required because of short
%  duration and low frequency resolution. 

if ~isfield(S,'FFT');            S.FFT = struct;           end
if ~isfield(param,'fitFFT');     param.fitFFT = true;      end % fit average FFT with Gaussian
if ~isfield(param,'plotFitFFT'); param.plotFitFFT = false; end % plot gaussian fit

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

% Simple peak determination - unreliable, sensitive to noise
[~, pkFreqInd]  = max(S.FFT.fftAve(S.FFT.subRange));
S.FFT.pkFreqAve = S.FFT.fftRange(S.FFT.subRange(pkFreqInd));

% Gaussian fit peak determination - fitting could cause warnings
if param.fitFFT
  ft             = fittype('gauss1');
  loParams       = [0 min(S.FFT.fftRange(S.FFT.subRange)) 0];
  hiParams       = [Inf max(S.FFT.fftRange(S.FFT.subRange)) range(S.FFT.fftRange(S.FFT.subRange))];
  S.FFT.fitGauss = fit(S.FFT.fftRange(S.FFT.subRange)', S.FFT.fftAve(S.FFT.subRange)', ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare');
  S.FFT.fitMean  = S.FFT.fitGauss.b1;
  S.FFT.fitSD    = S.FFT.fitGauss.c1/sqrt(2);
  S.FFT.fitFWHM  = 2 * sqrt(2 * log(2)) * S.FFT.fitSD;
  
  if param.plotFitFFT
    plot(S.FFT.fitGauss, S.FFT.fftRange(S.FFT.subRange)', S.FFT.fftAve(S.FFT.subRange))
  end
end

S.FFT = orderStruct(S.FFT);

end