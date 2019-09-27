function S = combPSQArea(dataFolder, saveFile)

% Handle optional arguments
if (nargin < 2) saveFile   = []; end
if (nargin < 1) dataFolder = []; end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

parentPath = [];

% If not supplied, prompt for folders to analyze
if isempty(dataFolder)
  dataFolder = uigetdir(pwd, 'Select folder containing analyzed PSC recordings');
end
if (dataFolder == 0) return; end

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
    [parentPath, saveFileName, ~] = parsePath(saveFile);
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

S.evPSQ{nFiles} = [];

for i = 1:nFiles
  data = load(dataFile{i});
  S.evPSQ{i} = data.C.SWR.evAreaAve;
end

% Calculate average evPSQ
S.evPSQAve = S.evPSQ{1};
for i = 2:nFiles
  S.evPSQAve = S.evPSQAve + S.evPSQ{i};
end
S.evPSQAve = S.evPSQAve/nFiles;

S.evPSQSEM = (S.evPSQ{1} - S.evPSQAve).*(S.evPSQ{1} - S.evPSQAve);
for i = 2:nFiles
  S.evPSQSEM = S.evPSQSEM + (S.evPSQ{i} - S.evPSQAve).*(S.evPSQ{i} - S.evPSQAve);
end
S.evPSQSEM = sqrt(S.evPSQSEM/(nFiles*(nFiles - 1)));

% Save file
if all(saveFileName)
  fprintf(['saving file ' saveFileName '... ']);
  save(saveFile,'-struct','S');
  fprintf('done\n');
end

fprintf('complete\n');

end
