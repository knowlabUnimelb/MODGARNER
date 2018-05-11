function plotCondition(mm, sm, pm, idx, plottitle, ylims, ylab)

figure('WindowStyle', 'docked')
elines = {' sk', ' dk', ' ok'};
lines = {':sk', '--dk', '-ok'};
plines = {' ok', ' ok', ' ok'};

items = [1 3 5];
for i = 1:numel(items)
    subplot(1,3,i)
    iset = mm{idx}(mm{idx}(:,2) == items(i), :); % Mean RT
    sset = sm{idx}(mm{idx}(:,2) == items(i), :); % Std RT
    pset = pm{idx}(pm{idx}(:,2) == items(i), :); % predicted RT

    temp = pset(:,3:end);
    bandwidth = getbandwidth(temp(:));
    [h2,L2,MX2,MED2,bw2] = violin(temp(:),...
        'facecolor', [1 1 1], 'edgecolor', 'k',...
        'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
    hold on
    
    errorbar(1:6, iset(:,3), sset(:,3)./sqrt(8), elines{i}, 'LineWidth', 1, 'MarkerSize', 10);
    h(i) = plot(1:6, iset(:,3), lines{i}, 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
    
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel(ylab, 'FontSize', 14)
    title(plottitle, 'HorizontalAlignment', 'center')
end
