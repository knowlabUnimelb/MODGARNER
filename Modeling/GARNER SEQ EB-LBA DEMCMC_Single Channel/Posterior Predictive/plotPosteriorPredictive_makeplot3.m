function plotPosteriorPredictive_makeplot3(mm, sm, pm, ma, sa, pa, varargin)
% Reorganize plot for posteriors
% Create separate figures for control, correlated, and filtering RT and
% accuracy.  Each figure is 1 x 3 subplots; each supblot contains different
% items (far, middle, near)
%
% pm should be cell matrix (cells are control, correlated, filtering)
%   each cell is 18 x NposteriorSamples (column 1 = current item, column 2
%   = previous item, columns 3-end are posterior samples)
%
% v3 - used N = 8 instead of N = 4
optargs = {[300 800], [0 1]};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[ylims, ylimsacc] = optargs{:}; % Place optional args in memorable variable names

%% Mean plots
plotCondition(mm, sm, pm, 1, 'Control', ylims, 'Mean RT');
plotCondition(mm, sm, pm, 2, 'Correlated', ylims, 'Mean RT')
plotFilteringCondition(mm, sm, pm, 3, 'Filtering', ylims, 'Mean RT');

%% Accuracy plots
plotCondition(ma, sa, pa, 1, 'Control', ylimsacc, 'Accuracy');
plotCondition(ma, sa, pa, 2, 'Correlated', ylimsacc, 'Accuracy')
plotFilteringCondition(ma, sa, pa, 3, 'Filtering', ylimsacc,'Accuracy');


%% Functions
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
    %
    %     for pidx = 1:size(pset,2)-2
    %         plot((1:6)+randn(1,6)*.1, pset(:,pidx+2), plines{i}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [.75 .75 .75], 'MarkerEdgeColor', [.75 .75 .75]); hold on
    %     end
    temp = pset(:,3:end);
    bandwidth = max(getbandwidth(temp(:)), .01);
    [h2,L2,MX2,MED2,bw2] = violin(pset(:,3:end)',...
        'facecolor', [1 1 1], 'edgecolor', 'k',...
        'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
    hold on
    
    errorbar(1:6, iset(:,3), sset(:,3)./sqrt(8), elines{i}, 'LineWidth', 1, 'MarkerSize', 10);
    hold on
    h(i) = plot(1:6, iset(:,3), lines{i}, 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
    
    
    set(gca,'XLim', [.5 6.5], 'XTick', 1:6, 'XTickLabel', {'far\newline', 'mid\newlineSame', 'near\newline', 'near\newline', 'mid\newlineOpp', 'far\newline'}, 'YLim', ylims)
    xlabel('Previous Item', 'FontSize', 14)
    ylabel(ylab, 'FontSize', 14)
    title(plottitle, 'HorizontalAlignment', 'center')
end

%%
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
            bandwidth = max(getbandwidth(temp(:)), .01);
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
            bandwidth = max(getbandwidth(temp(:)), .01);
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