clear all
clc

experiments = {'GarnerWithin_Bri_Sat_', 'GarnerSeparable_Bri_Sat_', 'GarnerBoxcar_Line_Sat_'}; % Fit experiments

fitfolder = fullfile(pwd, 'Fits');

dic1 = nan(1,numel(experiments));
dic2 = nan(1,numel(experiments));
for i = 1%:numel(experiments)
    if exist(fullfile(fitfolder, sprintf('%smeancentered.mat', experiments{i})), 'file') == 2
        load(fullfile(fitfolder, sprintf('%smeancentered.mat', experiments{i})), 'model', 'data', 'theta', 'weight', 'n')
        
        n.burnin = n.mc - 750;
        idx = n.burnin:n.mc;
        
        dic1(1,i) = computeDICgroup(data, theta, weight(:,n.burnin:end,:), 1, idx);
        dic2(1,i) = computeDICgroup(data, theta, weight(:,n.burnin:end,:), 2, idx);
        
        names = fieldnames(theta);
    end
end
