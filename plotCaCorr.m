function hand = plotCaCorr(data, hand, param)
%% hand = plotCaCorr(data, hand, param)
%  Function to plot event matrix, and SWR-SWR and Cell-Cell correlations from analyzeCaFile function

% Plot Event Matrix:
hand.CaEvFig = figure('Name',['SWR Event matrix for ' data.Ca.CaFile],'units','normalized');
hold on

if param.cellTypeOption
  indC  = zeros(param.nCellTypes, data.Ca.nChannels);
  for i = 1:param.nCellTypes
    indC(i,:) = (data.Ca.cellType == i);
    imagesc('XData', 1:size(data.SWR.Ca.evMatrixCorr,1), 'YData', find(indC(i,:)), 'CData', i*data.SWR.Ca.evMatrixCorr(:,logical(indC(i,:)))');
  end
else
  imagesc('XData', 1:size(data.SWR.Ca.evMatrixCorr,1), 'YData', 1:size(data.SWR.Ca.evMatrixCorr,2), 'CData', data.SWR.Ca.evMatrixCorr');
end

axis([0.5 size(data.SWR.Ca.evMatrixCorr,1) + 0.5 0.5 size(data.SWR.Ca.evMatrixCorr,2) + 0.5]);
caxis([0 256]);

evColMap = lines;
evColMap(1,:) = [1 1 1];
colormap(evColMap);

% Plot SWR-SWR Correlation Matrix:
hand.CaSWRCorrFig = figure('Name',['SWR-SWR correlation matrix for ' data.Ca.CaFile],'units','normalized');
imagesc('XData', 1:size(data.SWR.Ca.corrMatrix,1), 'YData', 1:size(data.SWR.Ca.corrMatrix,2), 'CData', data.SWR.Ca.corrMatrix');
axis([0.5 size(data.SWR.Ca.corrMatrix,1) + 0.5 0.5 size(data.SWR.Ca.corrMatrix,2) + 0.5]);
caxis([0 1]);
colormap(hot);
colorbar

% Plot Cell-Cell Correlation Matrix:
hand.CaCellCorrFig = figure('Name',['Cell-Cell correlation matrix for ' data.Ca.CaFile],'units','normalized');
imagesc('XData', 1:size(data.Ca.SWR.corrMatrix,1), 'YData', 1:size(data.Ca.SWR.corrMatrix,2), 'CData', data.Ca.SWR.corrMatrix');
axis([0.5 size(data.Ca.SWR.corrMatrix,1) + 0.5 0.5 size(data.Ca.SWR.corrMatrix,2) + 0.5]);
caxis([0 1]);
colormap(hot);
colorbar