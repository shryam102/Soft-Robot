
object_x = fpdockpredictedpositionssupershrunk1.x;
object_y = fpdockpredictedpositionssupershrunk1.y;

cobra_x  = rgbcobradocksynced.cobra_x;
cobra_y  = rgbcobradocksynced.cobra_y;

dock_x   = rgbcobradocksynced.dock_x;
dock_y   = rgbcobradocksynced.dock_y;


% figure('Color','w', 'Position',[300 300 600 600]);
% hold on;
% 

% plot(cobra_x, cobra_z, 'r--', 'LineWidth', 2);
% 

% plot(dock_x, dock_z, 'o', ...
%      'MarkerEdgeColor',[0 0.6 0], ...
%      'MarkerFaceColor','none', ...
%      'MarkerSize',12, ...
%      'LineWidth',2);
% 

% scatter(object_x, object_z, 10, 'b', 'filled', ...
%         'MarkerFaceAlpha',0.3, ...
%         'MarkerEdgeColor','b', ...
%         'MarkerEdgeAlpha',0.3, ...
%         'LineWidth',0.1);
% 
% hold off;
% 

% xlim([0.1   0.6]);     
% ylim([0.1   0.4]);  
% 

% axis square;
% grid on;
% 
% ax = gca;
% ax.LineWidth = 1.2;
% ax.FontSize  = 13;
% ax.GridAlpha = 0.3;
% ax.GridColor = [0 0 0];
% 

% yticks([0.1 0.25 0.4]);     
% xticks([0.1 0.35 0.6]);     
% 

% ax.XTickLabel = arrayfun(@num2str, ax.XTick, 'UniformOutput',false);
% ax.YTickLabel = arrayfun(@num2str, ax.YTick, 'UniformOutput',false);
% 

% xlabel('X','FontSize',14);
% ylabel('Z','FontSize',14);
% title('Docking Module Position Estimate in World Frame', ...
%       'FontSize',16, 'FontWeight','bold');
% legend({'Cobra Path','Dock Positions','Object Positions'}, ...
%        'Location','best', 'FontSize',12);






err_x = (object_x - dock_x).^2;
err_y = (object_y - dock_y).^2;
err_z = (object_z - dock_z).^2;
errors = [err_x(:), err_y(:), err_z(:)];


figure('Color','w','Position',[300 300 600 600]);


boxplot( ...
    errors, ...
    'Labels',{'X Variance','Y Variance','Z Variance'}, ...
    'Whisker',1.5, ...
    'Symbol','', ...            % no outlier markers
    'BoxStyle','outline', ...   % unfilled boxes
    'Colors',[0 0 0], ...       % black for all lines initially
    'Widths',0.6 ...
);
hold on;


boxes = findobj(gca,'Tag','Box');   % these are Line objects
set(boxes, 'Color', [0 0 1], 'LineWidth', 1.8);


medians = findobj(gca,'Tag','Median');
set(medians, 'Color', 'r', 'LineWidth', 2);



hold off;


axis square;
grid on;
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize  = 12;
ax.GridAlpha = 0.25;
ax.GridColor = [0 0 0];


xlabel('Coordinate','FontSize',14);
ylabel('Variance','FontSize',14);
title('Boxplot of Variance in XYZ Data Across Datasets','FontSize',15,'FontWeight','bold');

