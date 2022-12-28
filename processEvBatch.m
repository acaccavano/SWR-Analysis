function processEvBatch(param, evFolder, expFolder, statFolder, duration)
%% processPSCBatch(param, saveFolder, evFolder)
%
%  Function to clean and process batch of csv files of events from Clampfit.
%  Removes duplicate events (with same peak time), and calculates stats for
%  each trace. Not called from GUI, must be ran as standalone on csv folder.
%
%  Inputs: (all optional - will be prompted)
%   param      = structure containing all parameters - used for manual column selection (CAUTION - may change for different pClamp versions)
%   evFolder   = full path to pClamp psc event folder to import
%   expFolder  = full path to pClamp psc event folder to export (required, can be same to overwrite)
%   statFolder = full path to folder to save averaged trace stats (optional)
%   duration   = duration of trace in s

%% Handle input arguments - if not entered
if (nargin < 5); duration   = []; end
if (nargin < 4); statFolder = []; end
if (nargin < 3); expFolder  = []; end
if (nargin < 2); evFolder   = []; end
if (nargin < 1); param      = struct; end

% Handle case in which empty variable is supplied:
if isempty(param); param    = struct; end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% Select folder of PSC events exported from pClamp
if isempty(evFolder)
  evFolder = uigetdir(pwd, 'Select folder of event *.csv files exported from pClamp');
  if (evFolder == 0); error('No event folder selected'); end
  [parentPath, ~, ~] = parsePath(evFolder);
end

% Select folder to save processed PSC events
if isempty(expFolder)
  expFolder = uigetdir(parentPath, 'Select same or alternate folder to save processed *.csv files');
  if (expFolder == 0); error('No export folder selected'); end
  [parentPath, ~, ~] = parsePath(expFolder);
end

% Select folder to export processed PSC stats
if isempty(statFolder)
  statFolder = uigetdir(parentPath, 'Select folder to save processed event stats *.csv files');
  if (statFolder == 0); warning('No average stats will be exported - folder not selected'); end
end

% Define duration
if isempty(duration)
  answer = inputdlg({'Enter trace duration in s'}, 'Input', [1 35], {'60.000'});
  duration = str2double(answer{1});
end

% ensure current dir is in path so we can call helper funcs
curPath = pwd;
path(path, curPath);

% Extract input file names
cd (evFolder);
dir_temp  = dir('*.csv'); % Find only *.csv files
names     = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
evFiles  = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nEvFiles = length(evFiles);

cd (curPath);

% Determine individual file names:
evFile{nEvFiles}   = [];
expFile{nEvFiles}  = [];
statFile{nEvFiles} = [];

for i = 1:nEvFiles
  % Determine individual import files/folders
  evFile{i} = [evFolder slash evFiles{i}];
  [~, evFileName, ~] = parsePath(evFile{i});
  
  % Set export names
  expFile{i} = [expFolder slash evFileName '.csv'];
  
  % Determine individual exported PSC *.txt file names (if selected)
  if ~isempty(statFolder)
    if (statFolder ~= 0)
      statFile{i} = [statFolder slash evFileName '_Stats.csv'];
    end
  end 
end

for i = 1:nEvFiles
  processEvFile(param, evFile{i}, expFile{i}, statFile{i}, duration)
end
fprintf('complete\n');
end

