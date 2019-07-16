function [IEI1, IEI2] = histoIEI(evOption)
  % Event option can be 'swr', 'sw', or 'r'

  evName = [evOption 'IEI'];
  
  [fileName1, pathName1] = uigetfile('.mat','Open slice 1 output matlab file'); 
  file1 = [pathName1 fileName1];
  s1 = load(file1, evName);
  
  [fileName2, pathName2] = uigetfile('.mat','Open slice 2 output matlab file'); 
  file2 = [pathName2 fileName2];
  s2 = load(file2, evName);
  
  figure('Name',evName,'units','normalized','outerposition',[0 0.05 0.45 0.6], 'Color', [1 1 1]);
  h_ax = subplot('Position',[0.10 0.12 0.88 0.86]);
  hold(h_ax, 'on')
  binVector = linspace(0, ceil(max(max(s1.(evName)), max(s2.(evName)))), 31);
  histogram(s1.(evName), binVector, 'Normalization', 'pdf');
  histogram(s2.(evName), binVector, 'Normalization', 'pdf');
  hold(h_ax, 'off')
  
  set(h_ax,'FontSize',14);
  xlabel(h_ax,'seconds','FontSize',20);
  ylabel(h_ax,'Norm. Freq.','FontSize',20);
  h_leg = legend(fileName1,fileName2);
  set(h_leg,'FontSize',20);
  
end