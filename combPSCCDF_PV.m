function S = combPSCCDF_PV(param, dataFolders)
%% S = combPSCCDF(param, dataFolders)
%
%  Function to combine folder of *.mat files with PSC CDFs for amplitude, area, and kinetics

% Handle optional arguments
if (nargin < 2) dataFolders{6} = []; end
if (nargin < 1) param          = struct; end

% Handle case in which empty variable is supplied:
if isempty(param) param  = struct; end

% Set default parameters if not specified
if ~isfield(param,'nBins')    param.nBins    = 1000; end
if ~isfield(param,'minRise')  param.minRise  = 0.10; end 
if ~isfield(param,'maxRise')  param.maxRise  = 10;   end 
if ~isfield(param,'minDecay') param.minDecay = 3;    end 
if ~isfield(param,'maxDecay') param.maxDecay = 100;  end 

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

parentPath = [];

% If not supplied, prompt for folders to analyze
if isempty(dataFolders{1})
  dataFolders{1} = uigetdir(pwd, 'Select folder containing control PVAAC PSC recordings');
end
if (dataFolders{1} == 0) return; end

if isempty(dataFolders{2})
  dataFolders{2} = uigetdir(pwd, 'Select folder containing 5xFAD PVAAC PSC recordings');
end
if (dataFolders{2} == 0) return; end

if isempty(dataFolders{3})
  dataFolders{3} = uigetdir(pwd, 'Select folder containing control PVBC PSC recordings');
end
if (dataFolders{3} == 0) return; end

if isempty(dataFolders{4})
  dataFolders{4} = uigetdir(pwd, 'Select folder containing 5xFAD PVBC PSC recordings');
end
if (dataFolders{4} == 0) return; end

if isempty(dataFolders{5})
  dataFolders{5} = uigetdir(pwd, 'Select folder containing control PVBSC PSC recordings');
end
if (dataFolders{5} == 0) return; end

if isempty(dataFolders{6})
  dataFolders{6} = uigetdir(pwd, 'Select folder containing 5xFAD PVBSC PSC recordings');
end
if (dataFolders{6} == 0) return; end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Initialize output data structure:
S(6) = struct;
dataFile{6} = [];
nFiles = zeros(1,6);
for f = 1:6
  
  % Extract file names
  cd (dataFolders{f});
  dir_temp = dir('*.mat'); % Find only matlab files
  names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  nFiles(f) = length(file);
  cd (curPath);
  
  % Determine individual file names:
  dataFile{f}{nFiles(f)} = [];
  for i = 1:nFiles(f)
    dataFile{f}{i}  = [dataFolders{f} slash file{i}];
  end
  
  S(f).ampSpontX{nFiles(f)}   = [];
  S(f).ampSpontF{nFiles(f)}   = [];
  S(f).areaSpontX{nFiles(f)}  = [];
  S(f).areaSpontF{nFiles(f)}  = [];
  S(f).riseSpontX{nFiles(f)}  = [];
  S(f).riseSpontF{nFiles(f)}  = [];
  S(f).decaySpontX{nFiles(f)} = [];
  S(f).decaySpontF{nFiles(f)} = [];
  
  S(f).ampSWRX{nFiles(f)}    = [];
  S(f).ampSWRF{nFiles(f)}    = [];
  S(f).areaSWRX{nFiles(f)}   = [];
  S(f).areaSWRF{nFiles(f)}   = [];
  S(f).riseSWRX{nFiles(f)}   = [];
  S(f).riseSWRF{nFiles(f)}   = [];
  S(f).decaySWRX{nFiles(f)}  = [];
  S(f).decaySWRF{nFiles(f)}  = [];
  
  S(f).ampSpontAll    = [];
  S(f).areaSpontAll   = [];
  S(f).riseSpontAll   = [];
  S(f).decaySpontAll  = [];
  
  S(f).ampSWRAll      = [];
  S(f).areaSWRAll     = [];
  S(f).riseSWRAll     = [];
  S(f).decaySWRAll    = [];
  
  for i = 1:nFiles(f)
    data = load(dataFile{f}{i});
    
    S(f).ampSpontX{i}   = data.C.PSC.CDF.ampSpontX;
    S(f).ampSpontF{i}   = data.C.PSC.CDF.ampSpontF;
    S(f).areaSpontX{i}  = data.C.PSC.CDF.areaSpontX;
    S(f).areaSpontF{i}  = data.C.PSC.CDF.areaSpontF;
    S(f).riseSpontX{i}  = data.C.PSC.CDF.riseSpontX;
    S(f).riseSpontF{i}  = data.C.PSC.CDF.riseSpontF;
    S(f).decaySpontX{i} = data.C.PSC.CDF.decaySpontX;
    S(f).decaySpontF{i} = data.C.PSC.CDF.decaySpontF;
    
    S(f).ampSWRX{i}     = data.C.PSC.CDF.ampSWRX;
    S(f).ampSWRF{i}     = data.C.PSC.CDF.ampSWRF;
    S(f).areaSWRX{i}    = data.C.PSC.CDF.areaSWRX;
    S(f).areaSWRF{i}    = data.C.PSC.CDF.areaSWRF;
    S(f).riseSWRX{i}    = data.C.PSC.CDF.riseSWRX;
    S(f).riseSWRF{i}    = data.C.PSC.CDF.riseSWRF;
    S(f).decaySWRX{i}   = data.C.PSC.CDF.decaySWRX;
    S(f).decaySWRF{i}   = data.C.PSC.CDF.decaySWRF;
    
    % Exclude outliers for rise and decay (pClamp fitting often bad):
    exInd = data.C.PSC.riseTau < param.minRise | data.C.PSC.riseTau > param.maxRise;
    riseTau  = data.C.PSC.riseTau;
    riseTau(exInd) = [];
    swrMatRise = data.C.PSC.swrMatrix;
    swrMatRise(exInd) = [];
    
    exInd = data.C.PSC.decayTau < param.minDecay | data.C.PSC.decayTau > param.maxDecay;
    decayTau = data.C.PSC.decayTau;
    decayTau(exInd) = [];
    swrMatDecay = data.C.PSC.swrMatrix;
    swrMatDecay(exInd) = [];
    
    S(f).ampSpontAll   = vertcat(S(f).ampSpontAll, abs(data.C.PSC.amp(data.C.PSC.swrMatrix == 0)));
    S(f).areaSpontAll  = vertcat(S(f).areaSpontAll, abs(data.C.PSC.area(data.C.PSC.swrMatrix == 0)));
    S(f).riseSpontAll  = vertcat(S(f).riseSpontAll, abs(riseTau(swrMatRise == 0)));
    S(f).decaySpontAll = vertcat(S(f).decaySpontAll, abs(decayTau(swrMatDecay == 0)));
    
    S(f).ampSWRAll     = vertcat(S(f).ampSWRAll, abs(data.C.PSC.amp(data.C.PSC.swrMatrix == 1)));
    S(f).areaSWRAll    = vertcat(S(f).areaSWRAll, abs(data.C.PSC.area(data.C.PSC.swrMatrix == 1)));
    S(f).riseSWRAll    = vertcat(S(f).riseSWRAll, abs(riseTau(swrMatRise == 1)));
    S(f).decaySWRAll   = vertcat(S(f).decaySWRAll, abs(decayTau(swrMatDecay == 1)));
    
  end
  
  [S(f).ampSpontAllF, S(f).ampSpontAllX]     = ecdf(S(f).ampSpontAll);
  [S(f).areaSpontAllF, S(f).areaSpontAllX]   = ecdf(S(f).areaSpontAll);
  [S(f).riseSpontAllF, S(f).riseSpontAllX]   = ecdf(S(f).riseSpontAll);
  [S(f).decaySpontAllF, S(f).decaySpontAllX] = ecdf(S(f).decaySpontAll);
  
  [S(f).ampSWRAllF, S(f).ampSWRAllX]         = ecdf(S(f).ampSWRAll);
  [S(f).areaSWRAllF, S(f).areaSWRAllX]       = ecdf(S(f).areaSWRAll);
  [S(f).riseSWRAllF, S(f).riseSWRAllX]       = ecdf(S(f).riseSWRAll);
  [S(f).decaySWRAllF, S(f).decaySWRAllX]     = ecdf(S(f).decaySWRAll);
  
end

maxAmp = max([max(S(1).ampSpontAll) max(S(1).ampSWRAll) max(S(2).ampSpontAll) max(S(2).ampSWRAll) max(S(3).ampSpontAll) max(S(3).ampSWRAll) max(S(4).ampSpontAll) max(S(4).ampSWRAll) max(S(5).ampSpontAll) max(S(5).ampSWRAll) max(S(6).ampSpontAll) max(S(6).ampSWRAll)]);
maxArea = max([max(S(1).areaSpontAll) max(S(1).areaSWRAll) max(S(2).areaSpontAll) max(S(2).areaSWRAll) max(S(3).areaSpontAll) max(S(3).areaSWRAll) max(S(4).areaSpontAll) max(S(4).areaSWRAll) max(S(5).areaSpontAll) max(S(5).areaSWRAll) max(S(6).areaSpontAll) max(S(6).areaSWRAll)]);

for f = 1:6
  
  % Initialize standardized CDF arrays
  S(f).ampBins   = linspace(0, maxAmp, param.nBins)';
  S(f).areaBins  = linspace(0, maxArea, param.nBins)';
  S(f).riseBins  = linspace(param.minRise, param.maxRise, param.nBins)';
  S(f).decayBins = linspace(param.minDecay, param.maxDecay, param.nBins)';
  
  S(f).ampSpontCDF   = zeros(param.nBins, nFiles(f));
  S(f).areaSpontCDF  = zeros(param.nBins, nFiles(f));
  S(f).riseSpontCDF  = zeros(param.nBins, nFiles(f));
  S(f).decaySpontCDF = zeros(param.nBins, nFiles(f));
  
  S(f).ampSWRCDF     = zeros(param.nBins, nFiles(f));
  S(f).areaSWRCDF    = zeros(param.nBins, nFiles(f));
  S(f).riseSWRCDF    = zeros(param.nBins, nFiles(f));
  S(f).decaySWRCDF   = zeros(param.nBins, nFiles(f));
  
  S(f).ampSpontAllCDF   = zeros(param.nBins, 1);
  S(f).areaSpontAllCDF  = zeros(param.nBins, 1);
  S(f).riseSpontAllCDF  = zeros(param.nBins, 1);
  S(f).decaySpontAllCDF = zeros(param.nBins, 1);
  
  S(f).ampSWRAllCDF     = zeros(param.nBins, 1);
  S(f).areaSWRAllCDF    = zeros(param.nBins, 1);
  S(f).riseSWRAllCDF    = zeros(param.nBins, 1);
  S(f).decaySWRAllCDF   = zeros(param.nBins, 1);
  
  % Standardize combined CDFs:
  S(f).ampSpontAllCDF   = standardizeCDF(S(f).ampSpontAllX, S(f).ampSpontAllF, S(f).ampBins);
  S(f).areaSpontAllCDF  = standardizeCDF(S(f).areaSpontAllX, S(f).areaSpontAllF, S(f).areaBins);
  S(f).riseSpontAllCDF  = standardizeCDF(S(f).riseSpontAllX, S(f).riseSpontAllF, S(f).riseBins);
  S(f).decaySpontAllCDF = standardizeCDF(S(f).decaySpontAllX, S(f).decaySpontAllF, S(f).decayBins);
  
  S(f).ampSWRAllCDF   = standardizeCDF(S(f).ampSWRAllX, S(f).ampSWRAllF, S(f).ampBins);
  S(f).areaSWRAllCDF  = standardizeCDF(S(f).areaSWRAllX, S(f).areaSWRAllF, S(f).areaBins);
  S(f).riseSWRAllCDF  = standardizeCDF(S(f).riseSWRAllX, S(f).riseSWRAllF, S(f).riseBins);
  S(f).decaySWRAllCDF = standardizeCDF(S(f).decaySWRAllX, S(f).decaySWRAllF, S(f).decayBins);
  
  % Standardize individual file CDFs:
  for i = 1:nFiles(f)
    S(f).ampSpontCDF(:, i)    = standardizeCDF(S(f).ampSpontX{i}, S(f).ampSpontF{i}, S(f).ampBins);
    S(f).areaSpontCDF(:, i)   = standardizeCDF(S(f).areaSpontX{i}, S(f).areaSpontF{i}, S(f).areaBins);
    S(f).riseSpontCDF(:, i)   = standardizeCDF(S(f).riseSpontX{i}, S(f).riseSpontF{i}, S(f).riseBins);
    S(f).decaySpontCDF(:, i)  = standardizeCDF(S(f).decaySpontX{i}, S(f).decaySpontF{i}, S(f).decayBins);
    
    S(f).ampSWRCDF(:, i)    = standardizeCDF(S(f).ampSWRX{i}, S(f).ampSWRF{i}, S(f).ampBins);
    S(f).areaSWRCDF(:, i)   = standardizeCDF(S(f).areaSWRX{i}, S(f).areaSWRF{i}, S(f).areaBins);
    S(f).riseSWRCDF(:, i)   = standardizeCDF(S(f).riseSWRX{i}, S(f).riseSWRF{i}, S(f).riseBins);
    S(f).decaySWRCDF(:, i)  = standardizeCDF(S(f).decaySWRX{i}, S(f).decaySWRF{i}, S(f).decayBins);
  end
  
  % Average individual file CDFs:
  S(f).ampSpontAveCDF   = mean(S(f).ampSpontCDF, 2);
  S(f).areaSpontAveCDF  = mean(S(f).areaSpontCDF, 2);
  S(f).riseSpontAveCDF  = mean(S(f).riseSpontCDF, 2);
  S(f).decaySpontAveCDF = mean(S(f).decaySpontCDF, 2);
  
  S(f).ampSWRAveCDF     = mean(S(f).ampSWRCDF, 2);
  S(f).areaSWRAveCDF    = mean(S(f).areaSWRCDF, 2);
  S(f).riseSWRAveCDF    = mean(S(f).riseSWRCDF, 2);
  S(f).decaySWRAveCDF   = mean(S(f).decaySWRCDF, 2);
  
  % Calculate SEM:
  S(f).ampSpontSEMCDF   = std(S(f).ampSpontCDF, 0, 2) / sqrt(size(S(f).ampSpontCDF, 2));
  S(f).areaSpontSEMCDF  = std(S(f).areaSpontCDF, 0, 2) / sqrt(size(S(f).areaSpontCDF, 2));
  S(f).riseSpontSEMCDF  = std(S(f).riseSpontCDF, 0, 2) / sqrt(size(S(f).riseSpontCDF, 2));
  S(f).decaySpontSEMCDF = std(S(f).decaySpontCDF, 0, 2) / sqrt(size(S(f).decaySpontCDF, 2));
  
  S(f).ampSWRSEMCDF     = std(S(f).ampSWRCDF, 0, 2) / sqrt(size(S(f).ampSWRCDF, 2));
  S(f).areaSWRSEMCDF    = std(S(f).areaSWRCDF, 0, 2) / sqrt(size(S(f).areaSWRCDF, 2));
  S(f).riseSWRSEMCDF    = std(S(f).riseSWRCDF, 0, 2) / sqrt(size(S(f).riseSWRCDF, 2));
  S(f).decaySWRSEMCDF   = std(S(f).decaySWRCDF, 0, 2) / sqrt(size(S(f).decaySWRCDF, 2));
  
end

fprintf('complete\n');

end


function cdfOut = standardizeCDF(cdfXIn, cdfFIn, bins)
% function to interpolate CFD to standardized nBins

  % Trim duplicate X-values (would cause interpolation to crash)
  [cdfX, indX] = unique(cdfXIn, 'last');
  cdfF = cdfFIn(indX);
  
  % Cap CDF with X,Y = (0,0) and (xMax,1) if not present
  if cdfX(1) ~= bins(1)
    cdfX = vertcat(bins(1), cdfX);
    cdfF = vertcat(0, cdfF);
  end
  
  if cdfX(end) ~= bins(end)
    cdfX = vertcat(cdfX, bins(end));
    cdfF = vertcat(cdfF, 1);
  end
  
  cdfOut = interp1(cdfX, cdfF, bins, 'previous');
  
end