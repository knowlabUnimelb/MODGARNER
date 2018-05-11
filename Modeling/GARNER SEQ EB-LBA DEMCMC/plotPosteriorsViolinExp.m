% Plot subject level violin plots
clear all
clc
close all

fitfns = {'GarnerWithin_Bri_Sat_meancentered.mat', 'GarnerSeparable_Bri_Sat_meancentered.mat', 'GarnerBoxcar_Line_Sat_meancentered.mat'};

fig1 = figure('WindowStyle', 'docked');      
fig2 = figure('WindowStyle', 'docked');
fig3 = figure('WindowStyle', 'docked');
fig4 = figure('WindowStyle', 'docked');
fig5 = figure('WindowStyle', 'docked');
      
for fidx = 1:numel(fitfns)      
    load(fullfile(pwd, 'Fits', fitfns{fidx}))

    n.mc = i-1; 
    n.burnin = n.mc - 750;

%% Plot violin plot for attention weights
figure(fig1)

% Baseline
w_a = phi.w_a(:,n.burnin:n.mc); 
w_b = phi.w_b(:,n.burnin:n.mc);
w = betastat(w_a(:), w_b(:));

wc_a = phi.wc_a(:,n.burnin:n.mc); 
wc_b = phi.w_b(:,n.burnin:n.mc);
wc = betastat(wc_a(:), wc_b(:));

a = [w'; wc']';

% Plot function
colors = [.5 0 0; 0 .5 0; 0 0 .5];
xidx = fidx; 
xval = [xidx - .15, xidx + .15];
bandwidth = max(getbandwidth(a(:)), .01);
[h1,L1,MX1,MED1,bw1] = violin(a, 'x', xval, ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
    hm1 = plot(xval(1), mean(a(:,1)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    hm1 = plot(xval(2), mean(a(:,2)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
figHand1(fidx) = h1(1);

if fidx == 3
xlabel('Parameter')
ylabel('Posterior Value')
title('Attention Weights')
set(gca, 'YLim', [0 .5], 'XLim', [0 4], 'XTick', sort([(1:3) - .15, (1:3) + .15]), 'XTickLabel', {'w', 'w_{cr}', 'w', 'w_{cr}', 'w', 'w_{cr}'})
legend(figHand1, 'Integral', 'Exp 1: Separable', 'Exp 2: Boxcars')
end


%% Plot violin plot for sensitivity
figure(fig2)

% Baseline
c = phi.c_mu(:,n.burnin:n.mc); c = c(:);
cf = phi.cf_mu(:,n.burnin:n.mc); cf = cf(:);

C = [c, cf];

% Plot function
xidx = fidx; 
xval = [xidx - .15, xidx + .15];
bandwidth = max(getbandwidth(C(:)), .01);
[h1,L1,MX1,MED1,bw1] = violin(C, 'x', xval, ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
    hm1 = plot(xval(1), mean(C(:,1)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    hm1 = plot(xval(2), mean(C(:,2)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
figHand2(fidx) = h1(1);

if fidx == 3
xlabel('Parameter')
ylabel('Posterior Value')
title('Sensitivity')
set(gca, 'YLim', [0 3], 'XLim', [0 4], 'XTick', sort([(1:3) - .15, (1:3) + .15]), 'XTickLabel', {'c', 'c_{f}', 'c', 'c_{f}', 'c', 'c_{f}'})
legend(figHand2,  'Integral', 'Exp 1: Separable', 'Exp 2: Boxcars')
end

%% Violin plot for biases
figure(fig3)
% % Baseline
b1_a = phi.b1_a(:,n.burnin:n.mc); 
b1_b = phi.b1_b(:,n.burnin:n.mc);
b1 = betastat(b1_a(:), b1_b(:));

b2_a = phi.b2_a(:,n.burnin:n.mc); 
b2_b = phi.b2_b(:,n.burnin:n.mc);
b2 = betastat(b2_a(:), b2_b(:));

b3_a = phi.b3_a(:,n.burnin:n.mc); 
b3_b = phi.b3_b(:,n.burnin:n.mc);
b3 = betastat(b3_a(:), b3_b(:));

b4_a = phi.b4_a(:,n.burnin:n.mc); 
b4_b = phi.b4_b(:,n.burnin:n.mc);
b4 = betastat(b4_a(:), b4_b(:));

b5_a = phi.b5_a(:,n.burnin:n.mc); 
b5_b = phi.b5_b(:,n.burnin:n.mc);
b5 = betastat(b5_a(:), b5_b(:));

b = [b1'; b2'; b3'; b4'; b5']';

% Plot function
xs = [2 5 8];
xidx = xs(fidx); 
xval = [xidx - .75, xidx - .375, xidx, xidx + .375, xidx + .75];
bandwidth = max(getbandwidth(b(:)), .01);
[h1,L1,MX1,MED1,bw1] = violin(b, 'x', xval, ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
figHand3(fidx) = h1(1);

if fidx == 3
xlabel('Parameter')
ylabel('Posterior Value')
title('Category Bias')
set(gca, 'YLim', [0 1], 'XLim', [0 10], 'XTick', sort([(xs) - .75, (xs) - .375, (xs), (xs) + .375, (xs) + .75]),...
    'XTickLabel', {'b_{1}', 'b_{2}', 'b_{3}', 'b_{4}', 'b_{5}', 'b_{1}', 'b_{2}', 'b_{3}', 'b_{4}', 'b_{5}', 'b_{1}', 'b_{2}', 'b_{3}', 'b_{4}', 'b_{5}', 'b_{1}', 'b_{2}', 'b_{3}', 'b_{4}', 'b_{5}'})
legend(figHand3,  'Integral', 'Exp 1: Separable', 'Exp 2: Boxcars')
end

%% Plot violin plot for sequential parameters
figure(fig4);

% % Baseline
wp_a = phi.wp_a(:,n.burnin:n.mc); 
wp_b = phi.wp_b(:,n.burnin:n.mc);
wp = betastat(wp_a(:), wp_b(:));
% wp(wp < (-3 * std(wp) + mean(wp)) | wp > (3 * std(wp) + mean(wp))) = [];
% 
bh = phi.bh_mu(:,n.burnin:n.mc); bh = bh(:);
% bh(bh < (-3 * std(bh) + mean(bh)) | bh > (3 * std(bh) + mean(bh))) = [];
% 
adj = [-.5 0 .5];
bandwidth = max(getbandwidth(wp(:)), .01);
[h1,L1,MX1,MED1,bw1] = violin(wp, 'x', adj(fidx)+1, ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(adj(fidx)+1, mean(wp), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

bandwidth = max(getbandwidth(bh(:)), .01);
[h2,L1,MX1,MED1,bw1] = violin(bh, 'x', adj(fidx)+3, ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(adj(fidx)+3, mean(bh), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
figHand4(fidx) = h1(1);

xlabel('Parameter')
ylabel('Posterior Value')
title('Sequential Parameters')
set(gca, 'YLim', [0 2], 'XLim', [0 4], 'XTick', [1 3], 'XTickLabel', {'Alpha', 'Beta'})
line([0 4], [1 1], 'LineStyle', '--', 'Color', 'k')
if fidx ==3
legend(figHand4,  'Integral', 'Exp 1: Separable', 'Exp 2: Boxcars')
end


%% LBA parms
figure(fig5)

A = phi.A_mu(:,n.burnin:n.mc); A = A(:); A(A < (-3 * std(A) + mean(A)) | A > (3 * std(A) + mean(A))) = [];
bMa1 = phi.bMa1_mu(:,n.burnin:n.mc); bMa1 = bMa1(:); bMa1(bMa1 < (-3 * std(bMa1) + mean(bMa1)) | bMa1 > (3 * std(bMa1) + mean(bMa1))) = [];
bMa2 = phi.bMa2_mu(:,n.burnin:n.mc);bMa2 = bMa2(:); bMa2(bMa2 < (-3 * std(bMa2) + mean(bMa2)) | bMa2 > (3 * std(bMa2) + mean(bMa2))) = [];
vs = phi.vs_mu(:,n.burnin:n.mc);vs = vs(:);  vs(vs < (-3 * std(vs) + mean(vs)) | vs > (3 * std(vs) + mean(vs))) = [];
t0 = phi.t0_mu(:,n.burnin:n.mc); t0 = t0(:); t0(t0 < (-3 * std(t0) + mean(t0)) | t0 > (3 * std(t0) + mean(t0))) = [];
 
adj = [-.5 0 .5];
bandwidth = max(getbandwidth(A(:)), .01);
[h1,L1,MX1,MED1,bw1] = violin(A(:), 'x', 1+adj(fidx), ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(1+adj(fidx), mean(A(:)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
figHand5(fidx) = h1(1);

bandwidth = max(getbandwidth(bMa1(:)), .01);
[h2,L1,MX1,MED1,bw1] = violin(bMa1(:), 'x', 3+adj(fidx), ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(3+adj(fidx), mean(bMa1(:)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

bandwidth = max(getbandwidth(bMa2(:)), .01);
[h3,L1,MX1,MED1,bw1] = violin(bMa2(:), 'x', 5+adj(fidx), ...
    'facecolor',  colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(5+adj(fidx), mean(bMa2(:)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

bandwidth = max(getbandwidth(vs(:)), .01);
[h4,L1,MX1,MED1,bw1] = violin(vs(:), 'x', 7+adj(fidx), ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(7+adj(fidx), mean(vs(:)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

bandwidth = max(getbandwidth(t0(:)), .01);
[h4,L1,MX1,MED1,bw1] = violin(t0(:), 'x', 9+adj(fidx), ...
    'facecolor', colors(fidx,:), 'edgecolor', 'k',...
    'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
hold on
hm1 = plot(9+adj(fidx), mean(t0(:)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

xlabel('Parameter')
ylabel('Posterior Value')
title('LBA Parameters')
set(gca, 'YLim', [0 1], 'XLim', [0 10], 'XTick', 1:2:9, 'XTickLabel', {'A', 'T_A-A', 'T_B-A', 's_{\nu}', 't0'})
if fidx == 3
legend(figHand5,  'Integral', 'Exp 1: Separable', 'Exp 2: Boxcars')
end
end
