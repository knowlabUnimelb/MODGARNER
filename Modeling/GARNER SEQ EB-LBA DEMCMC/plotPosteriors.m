plotBurnIn = false;
plotIndividualFigs = false;
% endi = i;
endi = size(weight, 2);
n.thin = 100;
names = fieldnames(theta);
hnames = fieldnames(phi);
%% Plot subject level parms
scnt = 1;
for j = 1:n.subjects
    eval(sprintf('fig%d = figure(''WindowStyle'', ''docked'');', scnt));
    eval(sprintf('fig%d = figure(''WindowStyle'', ''docked'');', scnt+1));
    
    [nr, nc] = nsubplots(n.subLevParms );
    cnt = 1;
    for i = 1:n.subLevParms
        eval(sprintf('figure(fig%d);', scnt))
        subplot(nr, nc, cnt); 
        
        plot(theta.(names{i}){j}')
        title(names{i})
        ylims = get(gca, 'YLim');
        line([n.burnin n.burnin], ylims, 'LineStyle', '--')
        set(gca,'XLim', [0 endi]);
        
        eval(sprintf('figure(fig%d);', scnt+1))
        subplot(nr, nc, cnt); cnt = cnt + 1;
        temp = theta.(names{i}){j}(:,n.burnin:n.thin:end);
        [c, e] = hist(temp(:), 20);
        h = bar(e, c./sum(c));
        set(h, 'FaceColor', 'w')
        title(names{i})
    end
    scnt = scnt + 2;
end

%% Plot one figure with all of the hyperparameters in it
    eval(sprintf('fig%d = figure(''WindowStyle'', ''docked'');', scnt));
    eval(sprintf('fig%d = figure(''WindowStyle'', ''docked'');', scnt+1));
[nr, nc] = nsubplots(n.hyperParms);
cnt = 1;
for i = 1:n.hyperParms
    eval(sprintf('figure(fig%d);', scnt))
    subplot(nr, nc, cnt); 
    %     plot(theta.(names{i})', ' .k')
    plot(phi.(hnames{i})')
    title(hnames{i})
    ylims = get(gca, 'YLim');
    line([n.burnin n.burnin], ylims, 'LineStyle', '--')
    set(gca,'XLim', [0 endi]);
    
    eval(sprintf('figure(fig%d);', scnt+1))
    subplot(nr, nc, cnt); cnt = cnt + 1;
    temp = phi.(hnames{i})(:,n.burnin:end);
    [c, e] = hist(temp(:), 20);
    h = bar(e, c./sum(c));
    set(h, 'FaceColor', 'w')
    title(hnames{i})
end

%% Plot one figure for each parm
if plotIndividualFigs
    cnt = 1;
    for i = 1:n.hyperParms
        figure('WindowStyle', 'docked')
        subplot(1, 2, 1); cnt = cnt + 1;
        %     subplot(nr, nc, cnt); cnt = cnt + 1;
        
        if plotBurnIn
            plot(1:n.burnin, phi.(hnames{i})(:,1:n.burnin)')
        else
            plot(1:n.mcsamples, phi.(hnames{i})(:,n.burnin+1:end)')
        end
        title(hnames{i})
        
        subplot(1,2,2);
        if plotBurnIn
            temp = phi.(hnames{i})(:,1:n.burnin);
        else
            temp = phi.(hnames{i})(:,n.burnin:end);
        end
        [c, e] = hist(temp(:), 20);
        h = bar(e, c./sum(c));
        set(h, 'FaceColor', 'w')
        
    end
end