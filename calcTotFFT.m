function S = calcTotFFT(S, param)
%% S = calcTotFFT(S, param)
%
%  Function to create new structure and calculate fast-fourier transform
%  (FFT) over the total duration of a signal.

if ~isfield(S,'FFT');            S.FFT = struct;           end
if ~isfield(param,'fftRes');     param.fftRes = 0.1;       end % desired frequency resolution
if ~isfield(param,'fitFFT');     param.fitFFT = true;      end % fit FFT with Gaussian
if ~isfield(param,'plotFitFFT'); param.plotFitFFT = false; end % plot gaussian fit

% (Re)Initialize data arrays:
S.FFT.F   = [];
S.FFT.F1  = [];
S.FFT.F2  = [];

% Determine number of sections to result in desired frequency resolution
nSample    = length(S.tSeries);
nfft       = 2^nextpow2((param.Fs/param.fftRes) - 1); % number of spectral lines --> power of 2 because radix-2-algorithm                     
S.FFT.fRes = param.Fs/(nfft + 1); % real frequency resolution depending on nfft
nSection   = floor(nSample/nfft); % number of sections

% Frequency axis:
S.FFT.fftRange = linspace(0,param.Fs,nfft+1);
S.FFT.fftRange = S.FFT.fftRange(1:nfft/2+1);

% Divide the data into nSections with length of nfft
U = reshape(S.tSeries(1:nSection*nfft), [nfft, nSection]);

% Window function
win = hann(nfft,'Periodic')';   % create a hanning window with nfft values
W = repmat(win, nSection, 1)';  % create "window" matrix to multiply with each section
k = nfft/sum(win,2);            % correction coefficient of window function (not sure if correct)
Uw = U .* W;                    % values of the section multiplied by window function

% Calculate the frequency spectrum
S.FFT.F = fft(Uw,nfft,1);
S.FFT.F2 = k*abs(S.FFT.F); % Multiply with correction coefficient for window function
S.FFT.F1 = (1/(nSection*nfft))*sum(S.FFT.F2, 2);  % Average the sections
S.FFT.F1 = 2*S.FFT.F1(1:nfft/2+1); % Multiply the values by 2 --> two spectral lines make up one amplitude                       
S.FFT.F1(1) = S.FFT.F1(1)/2; % Only the first nfft/2+1 values needed because of multiplication by 2

S.FFT.subRange = find(S.FFT.fftRange>S.lim1 & S.FFT.fftRange<S.lim2);

% Simple peak determination - unreliable, sensitive to noise
[~, pkFreqInd] = max(S.FFT.F1(S.FFT.subRange));
S.FFT.pkFreq = S.FFT.fftRange(S.FFT.subRange(pkFreqInd));

% Gaussian fit peak determination - fitting could cause warnings
if param.fitFFT
  ft             = fittype('gauss1');
  loParams       = [0 min(S.FFT.fftRange(S.FFT.subRange)) 0];
  hiParams       = [Inf max(S.FFT.fftRange(S.FFT.subRange)) range(S.FFT.fftRange(S.FFT.subRange))];
  S.FFT.fitGauss = fit(S.FFT.fftRange(S.FFT.subRange)', S.FFT.F1(S.FFT.subRange), ft, 'Lower', loParams, 'Upper', hiParams, 'Robust', 'Bisquare');
  S.FFT.fitMean  = S.FFT.fitGauss.b1;
  S.FFT.fitSD    = S.FFT.fitGauss.c1/sqrt(2);
  S.FFT.fitFWHM  = 2 * sqrt(2 * log(2)) * S.FFT.fitSD;
  
  if param.plotFitFFT
    plot(S.FFT.fitGauss, S.FFT.fftRange(S.FFT.subRange)', S.FFT.F1(S.FFT.subRange))
  end
end

S.FFT = orderStruct(S.FFT);

end