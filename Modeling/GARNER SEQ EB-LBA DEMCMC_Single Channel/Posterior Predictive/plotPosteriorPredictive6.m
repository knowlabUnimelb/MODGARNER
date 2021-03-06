% Plot seqEBLBA Posterior Predictive

% Ensure that R_METRIC and SHEPARD_EXPONENT in seqEBLBA.m are set
% appropriately for your experiment


clear all
clc
% close all
tic

experiment = 3; % 1 = Garner Within, 2 = Separable (Bri/Sat), 3 = Boxcar (Line/Sat), 4 = simulation
condition = 2; % 0 = all, 1 = first 4 subs, 2 = last 4 subs

%% Datafile information
% Folders
mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\GarnerDEMCMC_recode';
% mainfolder = 'C:\Users\littled.UNIMELB\Dropbox\Work\2017 Garner\GarnerDEMCMC_recode';
dataloc = fullfile(mainfolder, 'Data');
mdsloc = fullfile(mainfolder, 'MDS');
fitloc = fullfile(mainfolder, 'Fits'); 


%% Set up file information
switch experiment
    
    case 1
        subStr = {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'}; % GarnerWithin (1:4 - Brightness, 5:8 - Saturation)
        dataPrefix = 'Garner_Exp1_sub%s_ordered_tagged.dat'; % For sprintf
        if condition == 0
            error('Use plotPosteriorPredictive7.m')
        elseif condition == 1
            mdsfn = 'mds_Brightness - Constrained.dat'; 
            subIdxs = 1:4; % :4; 
%             fitfn  = 'GarnerWithin_Bri_EBGaussian.mat'; % Experiment 1 - Integral Brightness
            fitfn = 'GarnerWithin_Brightness_meancentered.mat';
        elseif condition == 2
            mdsfn = 'mds_Saturation - Constrained.dat'; 
            subIdxs = 5:8;
%             fitfn  = 'GarnerWithin_Sat_EBGaussian.mat'; % Experiment 1 - Integral Saturation
            fitfn = 'GarnerWithin_Saturation_meancentered.mat';
%             fitfn = 'GarnerWithin_Saturation.mat';
        elseif condition == 3
            mdsfn = 'mds_Saturation - Constrained.dat'; 
            subIdxs = 5;%:8;
            fitfn = 'GarnerWithin_Saturation_102.mat'; % Sep Experiment 1 - Brightness
        end
    case 2
        subStr = {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'}; % Separable (1:4 - Brightness, 5:8 - Saturation)
        dataPrefix = 'Garner_SEP_sub%s_ordered_tagged.dat'; % For sprintf
        mdsfn = 'mds_Separable.dat';
        if condition == 1
             subIdxs = 1:4;
            fitfn = 'GarnerSeparable_Brightness_meancentered.mat'; % Sep Experiment 1 - Brightness
        elseif condition == 2
             subIdxs = 5:8;
            fitfn = 'GarnerSeparable_Saturation_meancentered_mds.mat'; % Sep Experiment 1 - Brightness
        end
    case 3
        subStr = {'L1', 'L2', 'L3', 'L4', 'S1', 'S2', 'S3', 'S4'}; % Boxcar (1:4 - Line Position, 5:8 - Saturation)
        dataPrefix = 'Garner_BOXCAR_sub%s_ordered_tagged.dat'; % For sprintf
        mdsfn = 'mds_Separable.dat';
        if condition == 1
             subIdxs = 1:4;
            fitfn = 'GarnerBoxcar_Line_meancentered.mat';
        else
             subIdxs = 5:8;
            fitfn = 'GarnerBoxcar_Saturation_meancentered_mds.mat';
        end
    case 4
        subStr = {'S01', 'S02', 'S03', 'S04', 'S05'};
        dataPrefix = 'Garner_SIM2_sub%s_ordered_tagged.dat'; % For sprintf
        mdsfn = 'mds_Saturation - Constrained.dat'; 
        subIdxs = 1:5;
        fitfn = 'GarnerModelRecovery_new2.mat';
end

reversexy = false;
if (experiment == 1 && condition == 2 ) || (experiment == 4 && condition == 1 )
    reversexy = true;
end


nTaskConditions = 3; % Control, correlated, filtering
nPosteriorPredictionSamples = 200;

%% Get data averages for plotting
% See getSubjectDataAverages.m
% load('GarnerWithin_Sat.mat', 'mm', 'sm', 'cm', 'ma', 'sa', 'ca')
[means, accuracy, subcheck] = getSubjectDataAverages(dataloc, dataPrefix, subStr, subIdxs);

%% Load posterior samples
load(fullfile(fitloc, fitfn), 'theta', 'phi', 'n', 'names', 'hnames')

%% The current is that the posterior predictive distribution is too variable. This 
% seems to have arisen from two sources - the first is that there are
% actual individual differences in overall RT so mixing a set of posteriors
% with different modes seems more variable. This occurs regardless of
% whether I sample for each subject or for the overall posterior.
% To get around this, I could mean correct each individuals posterior and
% then average those. 

%% Sample from each subject level theta
% Theta is a structure with fields corresponding to each parameter
% Each field is a cell matrix with cells equal to the number of subjects
% Each cell contains a nChain x nPosteriorSamples matrix of sampled values
subjectNumberCellOrder = [1 3 5 7 2 4 6 8]; % Data from idx subjects are in the corresponding cells - this is condition order rather than subject number order
funIdx = @(x, start, finish, step)(reshape(x(:,(start+1):step:finish)', numel(x(:,(start+1):step:finish)), 1)); % Function to sample post burnin samples
funsamp = @(x)(datasample(x, nPosteriorPredictionSamples, 'Replace', true));                                    % Function to down sample the posterior to the number of predictive params

% Find the postburnin samples
postsamples = structfun(@(y)(cellfun(@(x)funIdx(x, n.burnin, n.mc, 10), y, 'UniformOutput', false)), theta, 'UniformOutput', false);

% Sample a small number from the postburnin samples
samples = structfun(@(y)(cellfun(@(x)funsamp(x), y, 'UniformOutput', false)), postsamples, 'UniformOutput', false);

% Set up matrix of parameters for each subject and repetition
subparms = nan(n.subLevParms, n.subjects, nPosteriorPredictionSamples);
for pidx = 1:nPosteriorPredictionSamples
    subparms(:,:,pidx) = cell2mat(struct2cell(getChain(samples, pidx, [], 'one'))); % parms x subject x repetitions matrix
end

%% Test model using subject trial information
model = @seqEBLBA;

% Loop through each subject
pred_means = cell(1,nTaskConditions);
pred_accuracy = cell(1,nTaskConditions);
for sidx = 1:numel(subIdxs);
    sub = subStr{subIdxs(sidx)};
    
    if iscell(mdsfn)
        [data, cols, mds] = getPosteriorPredictionTrialData(dataloc, sprintf(dataPrefix, sub),...
            mdsloc, mdsfn{sidx}, reversexy(sidx));
    else
        [data, cols, mds] = getPosteriorPredictionTrialData(dataloc, sprintf(dataPrefix, sub),...
            mdsloc, mdsfn, reversexy);
    end
     simrt   = nan(size(data,1), nPosteriorPredictionSamples);
     simacc = nan(size(data,1), nPosteriorPredictionSamples); 
     for pidx = 1:nPosteriorPredictionSamples
         [simrt(:,pidx), simacc(:,pidx)] =...
             model(subparms(:, sidx, pidx), mds, data, cols, true);
     end
     
     % Remove short and long simulated RTs
     simacc(simrt < .2 | simrt > 3) = nan;
     simrt(simrt < .2 | simrt > 3)  = nan;

%      simacc(simrt < .01 | simrt > 10) = nan;
%      simrt(simrt < .01 | simrt > 10)  = nan;
     
     % Get averages for each subjects
     [rtMeans, accuracyMeans] = getSubjectPredictionAverages(data, cols, simrt, simacc);
     
     % Append each cell as a 3-D matrix
     pred_means = cellfun(@(x, y)(cat(3, x, y)), pred_means, rtMeans, 'UniformOutput', false);
     pred_accuracy = cellfun(@(x, y)(cat(3, x, y)), pred_accuracy, accuracyMeans, 'UniformOutput', false);
end


%% Average data and predictions across subjects
favgSub = @(x)(mean(x, 3));
fstdSub = @(x)(std(x, [], 3));

% Average data
ma = cellfun(favgSub, accuracy, 'UniformOutput', false); % Accuracy 
sa = cellfun(fstdSub, accuracy, 'UniformOutput', false); % Accuracy std
ca = mat2cell(numel(subIdxs) * ones(1,nTaskConditions), 1, ones(1, nTaskConditions)); % Count
mm = cellfun(favgSub, means, 'UniformOutput', false); % RT 
sm = cellfun(fstdSub, means, 'UniformOutput', false); % RT std 
cm = ca; 

% Average predictions
pred_ma = cellfun(favgSub, pred_accuracy, 'UniformOutput', false); % Accuracy 
pred_sa = cellfun(fstdSub, pred_accuracy, 'UniformOutput', false); % Accuracy std
pred_mm = cellfun(favgSub, pred_means, 'UniformOutput', false); % RT 
pred_sm = cellfun(fstdSub, pred_means, 'UniformOutput', false); % RT std 

if mm{1}(1,3) < 10
    mm = cellfun(@(x)([x(:,1:2), x(:,3:end,:)*1000]), mm, 'UniformOutput', false); % Convert from sec to msec for plotting
    sm = cellfun(@(x)([x(:,1:2), x(:,3:end,:)*1000]), sm, 'UniformOutput', false); % Convert from sec to msec for plotting
end

if pred_mm{1}(1,3)< 10
    pred_mm = cellfun(@(x)([x(:,1:2), x(:,3:end,:)*1000]), pred_mm, 'UniformOutput', false); % Convert from sec to msec for plotting
    pred_sm = cellfun(@(x)([x(:,1:2), x(:,3:end,:)*1000]), pred_sm, 'UniformOutput', false); % Convert from sec to msec for plotting
end

plotPosteriorPredictive_makeplot2(mm, sm, pred_mm, ma, sa, pred_ma, [350 1200])