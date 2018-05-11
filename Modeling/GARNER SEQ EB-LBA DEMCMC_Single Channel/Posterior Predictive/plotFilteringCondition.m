function plotFilteringCondition(mm, sm, pm, idx, plottitle, ylims, ylab)

% Filtering task
elines{1} = {' sk', ' dk', ' ok'};
elines{2} = {' sk', ' dk', ' ok'};
lines{1} = {':sk', '--dk', '-ok'};
lines{2} = {':sk', '--dk', '-ok'};
plines{1} = {' ok', ' ok', ' ok'};
plines{2} = {' ok', ' ok', ' ok'};

offsets = [-.15, .15];

items = [1 3 5];
figure('WindowStyle', 'docked')
subplotOrder = [1 4 2 5 3 6];
spcnt = 1;
for i = 1:numel(items)
    item = items(i);
    sets = [1:2:12; 2:2:12];
    
    for j = 1:2
        subplot(2,3,subplotOrder(spcnt)); spcnt = spcnt + 1;
        iset = mm{3}(mm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
        sset = sm{3}(mm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
        pset = pm{3}(pm{3}(:,2) == items(i) & ismember(mm{3}(:,1), sets(j,:)), :);
        
        
        if j == 1
            %             for pidx = 1:size(pset,2)-2
            %                 plot((1:6)+randn(1,6)*.1, pset(:,pidx+2), plines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10,...
            %                     'MarkerFaceColor', [.75 .75 .75], 'MarkerEdgeColor', [.75 .75 .75]);hold on
            %             end
            
            temp = pset(:,3:end);
            bandwidth = getbandwidth(temp(:));
            [h2,L2,MX2,MED2,bw2] = violin(pset(:,3:end)',...
                'facecolor', [1 1 1], 'edgecolor', 'k',...
                'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
            hold on
            
            errorbar((1:6), iset(:,3), sset(:,3)./sqrt(8), elines{j}{i}, 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
            
            h(i) = plot((1:6), iset(:,3), lines{j}{i}, 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
            
            
        else
            %             for pidx = 1:size(pset,2)-2
            %                 plot((1:6)+randn(1,6)*.1, pset(:,pidx+2), plines{j}{i}, 'LineWidth', 2, 'MarkerSize', 10,...
            %                     'MarkerFaceColor', [.75 .75 .75], 'MarkerEdgeColor', [.75 .75 .75]);hold on
            %             end
            temp = pset(:,3:end);
            bandwidth = getbandwidth(temp(:));
            [h2,L2,MX2,MED2,bw2] = violin(pset(:,3:end)',...
                'facecolor', [1 1 1], 'edgecolor', 'k',...
                'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
            hold on
            
            errorbar((1:6), iset(:,3), sset(:,3)./sqrt(8), elines{j}{i}, 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
            
            h(i) = plot((1:6), iset(:,3), lines{j}{i}, 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
            
        end
        
        set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
        xlabel('Previous Item', 'FontSize', 14)
        ylabel(ylab, 'FontSize', 14)
        if j == 1; title('Filtering\newlineSame Irrelevant Value', 'HorizontalAlignment', 'center'); else title('Filtering\newlineDifferent Irrelevant Value', 'HorizontalAlignment', 'center'); end
    end
end