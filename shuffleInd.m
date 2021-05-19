function skip = shuffleInd(nSample, Fs, nShuffle)
% Function to find shuffled indeces for tSeries

minSkip  = Fs;
maxSkip  = nSample - Fs;
skip     = ceil(nSample .* rand(nShuffle*2, 1));
skip(skip > maxSkip) = [];
skip(skip < minSkip) = [];
skip = skip(1:nShuffle, 1);
