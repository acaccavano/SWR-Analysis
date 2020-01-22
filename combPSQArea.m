function S = combPSQArea(epscFolder, ipscFolder, epscIndex, ipscIndex, saveFile)
%% S = combPSQArea(epscFolder, ipscFolder, epscIndex, ipscIndex, saveFile)
%
%  Function to combine folder of EPSC and IPSC *.mat files to calculate average EPSQ and IPSQ curves

% Handle optional arguments
if (nargin < 5) saveFile   = []; end
if (nargin < 4) ipscIndex  = []; end
if (nargin < 3) epscIndex  = []; end
if (nargin < 2) ipscFolder = []; end
if (nargin < 1) epscFolder = []; end

% Optional additional place to specify file indeces - comment out if not needed

% PC-5xFAD
% epscIndex = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
% ipscIndex = [1 3 4 5 6 7 8 9 10 11 12 13 14 15];

% PC-Control
% epscIndex = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
% ipscIndex = [1 2 3 4 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

% PVAAC-5xFAD
% epscIndex = [5 6];
% ipscIndex = [1 2];

% PVAAC-Control
% epscIndex = [2 3 4 5];
% ipscIndex = [1 2 3 4];

% PVBC-5xFAD
% epscIndex = [2 3 9 11 12 13 14 15];
% ipscIndex = [1 2 3 4 5 6 7 8];

% PVBC-Control
% epscIndex = [1 7 8 9 10 11 12];
% ipscIndex = [1 3 4 5 6 7 8];

% PVBSC-5xFAD
epscIndex = [1 2 7 8];
ipscIndex = [1 2 3 4];

% PVBSC-Control
% epscIndex = [3 4];
% ipscIndex = [1 2];

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

parentPath = [];

% If not supplied, prompt for EPSC folder to analyze
if isempty(epscFolder)
  epscFolder = uigetdir(pwd, 'Select folder containing analyzed EPSC recordings');
end
if (epscFolder == 0) return; end

% Parse epscFolder to determine default save name
[parentPath, epscFolderName, ~] = parsePath(epscFolder);

% If not supplied, prompt for IPSC folder to analyze
if isempty(ipscFolder)
  ipscFolder = uigetdir(parentPath, 'Select folder containing analyzed IPSC recordings');
end
if (ipscFolder == 0) return; end

% Select folders to save analyzed matlab files
if isempty(saveFile)
  defaultPath = [parentPath epscFolderName '.mat'];
  [saveName, savePath] = uiputfile('.mat','Select file to save output matlab file', defaultPath);
  saveFile = [savePath saveName];
  if ~all(saveFile)
    warning('No file to be saved - no file selected');
  else
    [parentPath, saveFileName, ~] = parsePath(saveFile);
  end
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Initialize output data structure:
S = struct;

%% Import and average swrEPSQ:
% Extract file names
cd (epscFolder);
dir_temp = dir('*.mat'); % Find only matlab files
names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nEPSCFiles = length(file);
cd (curPath);

% Determine individual file names:
epscFile{nEPSCFiles}    = [];

for i = 1:nEPSCFiles
  epscFile{i}  = [epscFolder slash file{i}];
end

S.evEPSQ{nEPSCFiles} = [];

for i = 1:nEPSCFiles
  data = load(epscFile{i});
  S.evEPSQ{i} = data.C.SWR.evAreaAve;
end

% Calculate average evEPSQ
S.evEPSQAve = S.evEPSQ{1};
for i = 2:nEPSCFiles
  S.evEPSQAve = S.evEPSQAve + S.evEPSQ{i};
end
S.evEPSQAve = S.evEPSQAve/nEPSCFiles;

S.evEPSQSEM = (S.evEPSQ{1} - S.evEPSQAve).*(S.evEPSQ{1} - S.evEPSQAve);
for i = 2:nEPSCFiles
  S.evEPSQSEM = S.evEPSQSEM + (S.evEPSQ{i} - S.evEPSQAve).*(S.evEPSQ{i} - S.evEPSQAve);
end
S.evEPSQSEM = sqrt(S.evEPSQSEM/(nEPSCFiles*(nEPSCFiles - 1)));


%% Import and average swrIPSQ:
% Extract file names
cd (ipscFolder);
dir_temp = dir('*.mat'); % Find only matlab files
names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
file = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nIPSCFiles = length(file);
cd (curPath);

% Determine individual file names:
ipscFile{nIPSCFiles}    = [];

for i = 1:nIPSCFiles
  ipscFile{i}  = [ipscFolder slash file{i}];
end

S.evIPSQ{nIPSCFiles} = [];

for i = 1:nIPSCFiles
  data = load(ipscFile{i});
  S.evIPSQ{i} = data.C.SWR.evAreaAve;
end

% Calculate average evIPSQ
S.evIPSQAve = S.evIPSQ{1};
for i = 2:nIPSCFiles
  S.evIPSQAve = S.evIPSQAve + S.evIPSQ{i};
end
S.evIPSQAve = S.evIPSQAve/nIPSCFiles;

S.evIPSQSEM = (S.evIPSQ{1} - S.evIPSQAve).*(S.evIPSQ{1} - S.evIPSQAve);
for i = 2:nIPSCFiles
  S.evIPSQSEM = S.evIPSQSEM + (S.evIPSQ{i} - S.evIPSQAve).*(S.evIPSQ{i} - S.evIPSQAve);
end
S.evIPSQSEM = sqrt(S.evIPSQSEM/(nIPSCFiles*(nIPSCFiles - 1)));


%% Calculate E/I Ratios
% If not already specified set default IPSC and EPSC indeces to include all files
if isempty(epscIndex) && isempty(ipscIndex)
  if nEPSCFiles == nIPSCFiles
    epscIndex = 1:nEPSCFiles;
    ipscIndex = 1:nIPSCFiles;
  else % truncate to shortest
    warning('unequal number of files and no file indeces specified - likely mismatch of files for EI ratios');
    nFiles = min(nEPSCFiles, nIPSCFiles);
    epscIndex = 1:nFiles;
    ipscIndex = 1:nFiles;
  end
elseif length(epscIndex) ~= length(ipscIndex) % Index specified but with different amounts for each
  error('unequal number of files specified in index');
elseif max(epscIndex) > nEPSCFiles || max(ipscIndex) > nIPSCFiles % Index specified but exceeds nFiles
  error('specified index exceeds number of files');
end

nFiles = length(epscIndex);
S.EIRatio{nFiles} = [];
S.IERatio{nFiles} = [];

% Calculate EI and IE ratios:
for i = 1:nFiles
  S.EIRatio{i} = -S.evEPSQ{epscIndex(i)} ./ S.evIPSQ{ipscIndex(i)};
  S.IERatio{i} = -S.evIPSQ{ipscIndex(i)} ./ S.evEPSQ{epscIndex(i)};
end

% Calculate average EI and IE ratios:
S.EIRatioAve = S.EIRatio{1};
S.IERatioAve = S.IERatio{1};

for i = 2:nFiles
  S.EIRatioAve = S.EIRatioAve + S.EIRatio{i};
  S.IERatioAve = S.IERatioAve + S.IERatio{i};
end

S.EIRatioAve = S.EIRatioAve/nFiles;
S.IERatioAve = S.IERatioAve/nFiles;

S.EIRatioSEM = (S.EIRatio{1} - S.EIRatioAve).*(S.EIRatio{1} - S.EIRatioAve);
S.IERatioSEM = (S.IERatio{1} - S.IERatioAve).*(S.IERatio{1} - S.IERatioAve);

for i = 2:nFiles
  S.EIRatioSEM = S.EIRatioSEM + (S.EIRatio{i} - S.EIRatioAve).*(S.EIRatio{i} - S.EIRatioAve);
  S.IERatioSEM = S.IERatioSEM + (S.IERatio{i} - S.IERatioAve).*(S.IERatio{i} - S.IERatioAve);
end

S.EIRatioSEM = sqrt(S.EIRatioSEM/(nFiles*(nFiles - 1)));
S.IERatioSEM = sqrt(S.IERatioSEM/(nFiles*(nFiles - 1)));


%% Save file
if all(saveFileName)
  fprintf(['saving file ' saveFileName '... ']);
  save(saveFile,'-struct','S');
  fprintf('done\n');
end

fprintf('complete\n');

end
