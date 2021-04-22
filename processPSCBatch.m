function processPSCBatch(param, pscFolder, expFolder, statFolder, duration)
%% processPSCBatch(param, saveFolder, pscFolder)
%
%  Function to clean and process batch of csv files of PSCs from Clampfit.
%  Removes duplicate PSCs (with same peak time), and calculates stats for
%  each trace. Not called from GUI, must be ran as standalone on csv folder.
%
%  Inputs: (all optional - will be prompted)
%   param      = structure containing all parameters including - PLACEHOLDER: CURRENTLY NOT USED
%     param.duration             = 1 = Single Recording, 2 = Multiple/Batch analysis (disables plotting)
%     param.pscEventPolarity     = boolean flag to indicate polarity of PSCs (0: neg/EPSCs, 1: pos/IPSCs) (default = 0)
%   pscFolder     = full path to pClamp psc event folder to import
%   expFolder     = full path to pClamp psc event folder to export (required, can be same to overwrite)
%   statFolder    = full path to folder to save averaged trace stats (optional)
%   duration      = duration of trace in s

%% Handle input arguments

if (nargin < 5) duration   = []; end
if (nargin < 4) statFolder = []; end
if (nargin < 3) expFolder  = []; end
if (nargin < 2) pscFolder  = []; end
if (nargin < 1) param      = struct; end

% Handle case in which empty variable is supplied:
if isempty(param) param      = struct; end

% % Set default parameters if not specified - PLACEHOLDER: CURRENTLY NO PARAMS
% if ~isfield(param,'fileNum')              param.fileNum              = 2;   end
% if ~isfield(param,'pscEventPolarity')     param.pscEventPolarity     = 0;   end

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% Select folder of PSC events exported from pClamp
if isempty(pscFolder)
  pscFolder = uigetdir(pwd, 'Select folder of PSC event *.csv files exported from pClamp');
  if (pscFolder == 0) error('No PSC event folder selected'); end
  [parentPath, ~, ~] = parsePath(pscFolder);
end

% Select folder to save processed PSC events
if isempty(expFolder)
  expFolder = uigetdir(parentPath, 'Select same or alternate folder to save processed *.csv files');
  if (expFolder == 0) error('No export folder selected'); end
  [parentPath, ~, ~] = parsePath(expFolder);
end

% Select folder to export processed PSC stats
if isempty(statFolder)
  statFolder = uigetdir(parentPath, 'Select folder to save processed PSC stats *.csv files');
  if (statFolder == 0) warning('No average stats will be exported - folder not selected'); end
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
cd (pscFolder);
dir_temp  = dir('*.csv'); % Find only *.csv files
names     = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
pscFiles  = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
nPSCFiles = length(pscFiles);

cd (curPath);

% Determine individual file names:
pscFile{nPSCFiles}  = [];
expFile{nPSCFiles}  = [];
statFile{nPSCFiles} = [];

for i = 1:nPSCFiles
  % Determine individual import files/folders
  pscFile{i} = [pscFolder slash pscFiles{i}];
  [~, pscFileName, ~] = parsePath(pscFile{i});
  
  % Set export names
  expFile{i} = [expFolder slash pscFileName '.csv'];
  
  % Determine individual exported PSC *.txt file names (if selected)
  if ~isempty(statFolder)
    if (statFolder ~= 0)
      statFile{i} = [statFolder slash pscFileName '_stats.csv'];
    end
  end 
end

for i = 1:nPSCFiles
  A = 5;
  processPSCFile(param, pscFile{i}, expFile{i}, statFile{i}, duration)
end
fprintf('complete\n');
end

