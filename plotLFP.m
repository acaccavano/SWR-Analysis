% Script to make plotting a particular file faster
fprintf('plotting file... ');
dsPlot = 1;

% check if data file is open
if ~exist('data','var')
  [fileName, pathName] = uigetfile('.mat','Open output matlab file'); 
  folderPath = [pathName fileName];
  data = load(folderPath);
end

% Check if handle structure exists
if ~exist('hand','var')
  hand = struct;
end

% Plot figure
hand.lfpFig = figure('Name',['Raw and Filtered Data for ' data.LFP.dataFile],'units','normalized','outerposition',[0 0.05 1 0.95]);
hand = plotLFPAnalysis(data, [], hand, data.LFP.param, dsPlot);

fprintf('done\n');