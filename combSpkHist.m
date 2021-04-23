function S = combSpkHist(param, dataFolder, saveFile)
%% S = combSpkHist(param, dataFolder, saveFile)
%
%  Function to combine folder of spike *.mat files to calculate average peri-SWR spike histogram

% Handle optional arguments
if (nargin < 3); saveFile   = []; end
if (nargin < 2); dataFolder = []; end
if (nargin < 1); param      = struct; end

% Handle case in which empty variable is supplied:
if isempty(param); param    = struct; end

% Set default parameters if not specified
if ~isfield(param,'nBins'); param.nBins = 100; end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% If not supplied, prompt for folders to analyze
if isempty(dataFolder)
  dataFolder = uigetdir(pwd, 'Select folder containing analyzed Spike recordings');
end
if (dataFolder == 0); return; end

% Parse dataFolder to determine default save name
[parentPath, dataFolderName, ~] = parsePath(dataFolder);

% Select folders to save analyzed matlab files
if isempty(saveFile)
  defaultPath = [parentPath dataFolderName '.mat'];
  [saveName, savePath] = uiputfile('.mat','Select file to save output matlab file', defaultPath);
  saveFile = [savePath saveName];
  if ~all(saveFile)
    warning('No file to be saved - no file selected');
  else
    [~, saveFileName, ~] = parsePath(saveFile);
  end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract file names
cd (dataFolder);
dir_temp = dir('*.mat'); % Find only matlab files
names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nFiles = length(file);
cd (curPath);

% Determine individual file names:
dataFile{nFiles}    = [];

for i = 1:nFiles
  dataFile{i}  = [dataFolder slash file{i}];
end

% Initialize output data structure:
S = struct;

S.sHist         = [];
S.sHist{nFiles} = [];
S.nSWRs         = zeros(nFiles, 1);
S.nSpks         = zeros(nFiles, 1);

for i = 1:nFiles
  data = load(dataFile{i});
  
  S.nSWRs(i) = length(data.C.SWR.spike.evPeakA);
  
  Bins = linspace(data.SWR.evTiming(1), data.SWR.evTiming(length(data.SWR.evTiming)), param.nBins)';
  binWidth = Bins(2) - Bins(1);
  
  S.sHist{i} = zeros(length(Bins), 1);
  
  for swr = 1:S.nSWRs(i)
    nSpkSWR = length(data.C.SWR.spike.evPeakA{swr});
    S.nSpks(i) = S.nSpks(i) + nSpkSWR;
    for spk = 1:nSpkSWR
      peakTime = data.SWR.evTiming(data.C.SWR.spike.evPeakA{swr}(spk));
      S.sHist{i} = S.sHist{i} + double(peakTime >= Bins & peakTime < Bins + binWidth);
    end
  end
  S.sHist{i} = S.sHist{i}/S.nSWRs(i);
end

S.Bins = Bins;

% Calculate average histogram
S.sHistAve = S.sHist{1};
for i = 2:nFiles
  S.sHistAve = S.sHistAve + S.sHist{i};
end
S.sHistAve = S.sHistAve/nFiles;

S.sHistSEM = (S.sHist{1} - S.sHistAve).*(S.sHist{1} - S.sHistAve);
for i = 2:nFiles
  S.sHistSEM = S.sHistSEM + (S.sHist{i} - S.sHistAve).*(S.sHist{i} - S.sHistAve);
end
S.sHistSEM = sqrt(S.sHistSEM/(nFiles*(nFiles - 1)));

% Save file
if all(saveFileName)
  fprintf(['saving file ' saveFileName '... ']);
  save(saveFile,'-struct','S');
  fprintf('done\n');
end

fprintf('complete\n');

end
