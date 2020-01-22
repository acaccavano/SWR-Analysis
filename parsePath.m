function [parentPath, targetName, targetExt] = parsePath(targetPath)
%% [parentPath, targetName, targetExt] = parsePath(targetPath)
%
% Function parses an input file or folder path into three parts:
% the parent folder, the target file or folder, and the target extension,
% if target is a file. If target is a folder, targetExt = empty

% Assign OS specific variables:
if ispc
  slash = '\';
else
  slash = '/';
end

% Parse dataFile to determine default save name
slashPos   = strfind(targetPath, slash);
slashPos   = slashPos(size(slashPos, 2));
parentPath = targetPath(1 : slashPos);
targetName = targetPath(slashPos+1 : end);
perPos     = strfind(targetName,'.');
if isempty(perPos)
  targetExt = [];
else
  perPos = perPos(size(perPos, 2));
  targetExt  = targetName(perPos+1 : end);
  targetName = targetName(1 : perPos-1);
end

end

