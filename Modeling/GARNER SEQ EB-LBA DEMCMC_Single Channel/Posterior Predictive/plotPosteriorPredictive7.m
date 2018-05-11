% Plot seqEBLBA Posterior Predictive

% v7 - set up to plot using different mds for each subjects, which can be
% used for plotting all subjects from two conditions

clear all
clc
% close all
tic

experiment = 3; % 1 = Garner Within, 2 = Separable (Sat/Bri), 3 = Boxcar (Line/Sat), 4 = simulation
condition = 0; % 0 = all, 1 = first 4 subs, 2 = last 4 subs

%% Datafile information
% Folders
mainfolder = 'C:\Users\littled\Dropbox\Work\2017 Garner\GARNER SEQ EB-LBA DEMCMC_Single Channel';
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
            mdsfn = {'mds_Brightness - Constrained.dat', 'mds_Brightness - Constrained.dat', 'mds_Brightness - Constrained.dat', 'mds_Brightness - Constrained.dat',...
                     'mds_Saturation - Constrained.dat', 'mds_Saturation - Constrained.dat', 'mds_Saturation - Constrained.dat', 'mds_Saturation - Constrained.dat'}; 
            subIdxs = 1:8; 
            fitfn  = 'GarnerWithin_Bri_Sat_singleChannel.mat'; % Experiment 1 - Integral Brightness
        end
    case 2
        subStr = {'B1', 'B2', 'B3', 'B4', 'S1', 'S2', 'S3', 'S4'}; % Separable (1:4 - Brightness, 5:8 - Saturation)
        dataPrefix = 'Garner_SEP_sub%s_ordered_tagged.dat'; % For sprintf
        mdsfn = 'mds_Separable.dat';
        if condition == 0
            mdsfn = {'mds_Separable.dat','mds_Separable.dat','mds_Separable.dat','mds_Separable.dat',...
                'mds_Separable.dat','mds_Separable.dat','mds_Separable.dat','mds_Separable.dat'};
            subIdxs = 1:8;
            fitfn  = 'GarnerSeparable_Bri_Sat_singleChannel.mat'; % Experiment 1 - Integral Brightness
        end
    case 3
        subStr = {'L1', 'L2', 'L3', 'L4', 'S1', 'S2', 'S3', 'S4'}; % Boxcar (1:4 - Line Position, 5:8 - Saturation)
        dataPrefix = 'Garner_BOXCAR_sub%s_ordered_tagged.dat'; % For sprintf
        mdsfn = 'mds_Separable.dat';
        if condition == 0
            mdsfn = {'mds_Separable.dat','mds_Separable.dat','mds_Separable.dat','mds_Separable.dat',...
                'mds_Separable.dat','mds_Separable.dat','mds_Separable.dat','mds_Separable.dat'};
            subIdxs = 1:8;
            fitfn  = 'GarnerBoxcar_Line_Sat_singleChannel.mat'; % Experiment 1 - Integral Brightness            
        end
end

reversexy = false(1,8);
if (experiment == 1 && condition == 2 ) || (experiment == 4 && condition == 1 )
    reversexy = true;
elseif experiment == 1 && condition == 0
    reversexy = [false(1,4), true(1,4)];
end

nTaskConditions = 3; % Control, correlated, filtering
nPosteriorPredictionSamples = 200;

%% Get data averages for plotting
% See getSubjectDataAverages.m
% load('GarnerWithin_Sat.mat', 'mm', 'sm', 'cm', 'ma', 'sa', 'ca')
[means, accuracy, subcheck] = getSubjectDataAverages(dataloc, dataPrefix, subStr, subIdxs);

%% Load posterior samples
load(fullfile(fitloc, fitfn), 'theta', 'phi', 'n', 'names', 'hnames')

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
model = @seqEBLBA_single;

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