function dataOut = concatSweeps(inFile, outFile, param)

%% Handle input arguments - if not entered
if (nargin < 3) param   = struct; end
if (nargin < 2) outFile = []; end
if (nargin < 1) inFile  = []; end

% Handle case in which empty variables are supplied
if isempty(param) param = struct; end

% Set default parameters if not specified
if ~isfield(param,'fileType')      param.fileType     = 1;    end
if ~isfield(param,'Fs')            param.Fs           = 3000; end
if ~isfield(param,'dsFactor')      param.dsFactor     = 1;    end
if ~isfield(param,'channel')       param.channel      = [1 3 5]; end

nChannels = length(param.channel);

% If not supplied, prompt for files/folder to analyze
if isempty(inFile)
  if (param.fileType == 1)
    [fileName, filePath] = uigetfile('.abf', 'Select *.abf file to process');
    inFile = strcat(filePath,fileName);
  elseif (param.fileType == 2)
    inFile = uigetdir();
  end
  if ~all(inFile) return; end
end

% Parse inFile to determine default outFile name
[parentPath, inFileName, ~] = parsePath(inFile);

% If not supplied, prompt for outFile name
if isempty(outFile)
  defaultPath = [parentPath inFileName '.txt'];
  [outFileName, outFilePath] = uiputfile('.txt','Select file to save processed *.txt file', defaultPath);
  outFile = [outFilePath outFileName];
  if ~all(outFile) warning('No file to be saved - no file selected'); end
end

fprintf('importing data... ');
if (param.fileType == 1)
  [dataIn, samplingInt, headerInfo] = abfload(inFile);
  nSweeps = size(dataIn, 3);
  
  % Convert from um to ms
  samplingInt = samplingInt / 1000;
  
  % Take subset of channels specified:
  dataIn = dataIn(:, param.channel, :);
  
elseif (param.fileType == 2)
  
  samplingInt = 1000 / param.Fs; % (ms)
  
  % Ensure current directory is in path so helper functions work
  curPath = pwd;
  path(path, curPath);
  
  % Extract file names
  cd (inFile);
  dir_temp = dir('2*');
  names = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
  files = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."
  
  % Sort file names
  if iscell(files)
    for i = 1:length(files)
      temp = files{i};
      if ~contains(temp,'_')
        string = temp(1:strfind(temp,'.')-1);
      else
        if ~contains(temp,'.')
          string = temp(max(strfind(temp,'_'))+1:end);
        else
          string = temp(max(strfind(temp,'_'))+1:strfind(temp,'.')-1);
        end
      end
      % extract the number part
      str_ascii= str2num(int2str(string));
      index_nums = find(str_ascii <=57 & str_ascii >=48);
      filenum(i) = str2num(string(index_nums));
    end
    [filenum, index]=sort(filenum,'ascend');
    files = files(index);
  else
    temp_file = files;
    clear files;
    files = cell(1);
    files{1}=temp_file;
  end
  
  nSweeps = length(files);
  dataImport{nChannels, nSweeps} = []; % raw data cell array
  dataSize = zeros(nChannels, nSweeps); % file length for every channel and sweep
  
  %% Read Data file by file
  for sw = sw:nSweeps
    buf = fopen(files{sw},'r'); % open data file was files{i}
    firstline = fgetl(buf);
    while (isempty(firstline))
      firstline = fgetl(buf);
    end
    n = length(sscanf(firstline,'%g'));
    data_firstline = sscanf(firstline,'%g',[n inf]);
    
    dataTemp = fscanf(buf,'%g',[n inf]);
    dataTemp = [data_firstline dataTemp];
    dataTemp = dataTemp';
    
    fclose(buf);  % close data file after read
    
    % Import channel(s) of interest
    for ch = 1:nChannels
      dataImport{ch, sw} = dataTemp(:, param.channel(ch));
      dataSize(ch, sw)   = length(dataImport{ch, sw});
    end
  end
  
  % Check if last file is partial, and remove if so:
  if (dataSize(1, nSweeps) < max(dataSize(1, :)))
    dataImport{:, nSweeps} = [];
    dataSize(:, nSweeps)   = [];
    nSweeps = nSweeps - 1;
  end
  
  % Convert cell array to array, checking if filesize consistent
  if all(dataSize(:) == dataSize(1,1))
    dataIn = zeros(dataSize(1,1), nChannels, nSweeps);
    for ch = 1:nChannels
      for sw = 1:nSweeps
        dataIn(:, ch, sw) = dataImport{ch, sw};
      end
    end
  else
    error('Inconsistent file sizes, check for bad data file');
  end
  cd (curPath);
end
fprintf('done\n');

%% Concatenate sweeps
dataOut = dataIn(:, :, 1);
for sw = 2:nSweeps
  dataOut = vertcat(dataOut, dataIn(:, :, sw));
end

timing = 0 : samplingInt : samplingInt * (size(dataOut, 1) - 1);
timing = 0.5 + timing';
dataOut = horzcat(timing, dataOut);

% Output table names
tableVarNames = cell(1, nChannels);
nameIndex = 1;
tableVarNames{nameIndex} = 'time_ms';

% Sort cell array into output table
for ch = 1:nChannels
  nameIndex = nameIndex + 1;
  if (param.fileType == 1)
    tableVarNames{nameIndex} = [strrep(headerInfo.recChNames{param.channel(ch)},' ', '') '_' strrep(headerInfo.recChUnits{param.channel(ch)},' ', '')];
  elseif (param.fileType == 1)
    tableVarNames{nameIndex} = ['Ch' int2str(param.channel(ch))];
  end
end

%% Save output *.txt file
if all(outFile) 
  fprintf('saving combined file... ');
  dataTable = array2table(dataOut, 'VariableNames', tableVarNames);
  writetable(dataTable, outFile, 'Delimiter', '\t');
  fprintf('done\n');
end

end
