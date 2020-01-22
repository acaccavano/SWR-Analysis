function files = findFilesInDir(inFolder, fileType)
%% files = findFilesInDir(inFolder, fileType)
%
% function inputs drectory and filetype (entered as '*.mat' for example),
% and outputs cell array of file names

curPath = pwd; % Save curPath so it can be returned to

cd (inFolder);
dir_temp = dir(fileType); % Find only files of fileType
names    = {dir_temp.name}; % extract all the names in the struct returned by 'dir': ".", "..", file 1,2....
files    = names([dir_temp.isdir] == 0); % extract the name for all files, but no "." and ".."

cd (curPath);
end