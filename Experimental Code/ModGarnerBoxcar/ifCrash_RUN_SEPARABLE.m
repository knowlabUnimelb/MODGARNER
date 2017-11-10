%% GARNER WITHIN
% 160229 - DL: Modified to present separable dimensions
%
%
% JANUARY 2013 TONY WANG
%
% Categorisation of Integral Stimuli using multiple (>2) values on the
% saturation and brightneascas dimensions. This experiment uses the values of
% the Munsells colour space that was original proposed in the experimental
% design
%
% All five conditions of the Garner Design is implemented within each
% experimental session. The order in which participants experience each
% condition is randomly determine.
%
% Note, the experimenter must enter the schedule number (of 6 possible
% permutations) before each session. This determines the condition
% schedule/order. This could be randomly determined prior to the
% experiment.
%
% The experimental design of this experiment is similar to of previous
% Garner experiments. 12 stimuli are presented in each condition in order
% to avoid repetition.
%
% Created - 10/01/2013 Tony Wang
% Last Edited - 14/1/2013

%% Clear space
clear all
clc
close all
format short g
% global ptb3  % delete?
WaitSecs(1e-7); % Hack to load WaitSecs
warning('off', 'MATLAB:mex:deprecatedExtension')

% Get PsychtoolboxVersion and adjust display mode if ptb version = 3
x = PsychtoolboxVersion; if ischar(x); ptb = str2double(x(1)); else ptb = x; end; ptb3 = ptb >= 3;

debug = true; % If debug, do experiment with reduced number of trials
RTBox('fake', debug);

xup = -50; xleft = -18;  % center of screen offsets for showtext
bgcolor = [255 255 255]; % User specified background color
feedback = {'...Wrong...', '...Correct...'};
fbackLocation = [300, 0]; % Feedback Location

%% Experimental Variables
seed = 12699;                   % Seed for random number generator
fixationDuration = 1.5;         % Fixation cross presentation length
fbkDuration      = 2;           % Feedback presentation length
iti              = 1.5;         % Intertrial interval
timeout          = 5;

textsize  = 60;  % Text size
stimsize  = 100; % Stimulus size (height and width of circumscribing rectangle)
fixsize   = 50;  % Height and width of fixation cross (in pixels)
lineWidth = 5;   % Pixel width

% No of Experimental Trials
switch debug
    case true
        nPracTrials = 12;
        nExpTrials  = 12;
        nBlocks     = 6; % Refers to the total number of blocks in the exp (include practice trials)
    case false
        nPracTrials = 24; % Number of practice trials = 200
        nExpTrials = 120; % Number of experimental trials per block
        nBlocks    = 12;  % Number of blocks = 6
end


%% Open Experiment Window
% Present subject information screen until correct information is entered
subject = input('Enter Subject Number:');   % Subject no
condition = input('Enter Condition:');      % 1 = brightness relevant, 2 = saturation relevant
session = input('Enter Session:');          % Session no
startBlock = input('Start from block:');    % Which block to start from

outputfile = sprintf('cr2_GARNER_SEPARABLE_%03d_%02d_s%02d.dat', subject, condition, session);
tmpfile = sprintf('cr2_tmp_GARNER_SEPARABLE_%03d_%02d_s%02d.dat', subject, condition, session);
rng(seed + subject, 'twister'); % Seed random number generator for subject specific vars
brightleft = rand <= .5; % If true, then brightness varies on left with constant saturation, saturation varies on right with constant brightness; otherwise, reverse
sessionOrder = randperm(6);

rng(seed + subject * 2 + session, 'twister'); % Seed random number generator for subject/session vars

screenparms = prepexp(0, bgcolor); % Open onscreen window

%% Experimental Stimuli Values

% There are two conditions based on whether brightness (value) or 
% saturation (chroma) is relevant
% VALUE = BRIGHTNESS
% CHROMA = SATURATION

if condition == 1; % Brightness is relevant
    if brightleft
        left.brightness  = [3 4 5 6 7 8]; 
        right.brightness = 6;
        
        left.saturation  = 4;
        right.saturation = [6 8];
    else
        left.brightness  = 6;
        right.brightness = [3 4 5 6 7 8]; 
    
        left.saturation = [6 8];
        right.saturation = 4;
    end
else % Saturation is relevant
    if brightleft
        left.brightness = [5 6];
        right.brightness = 6;
        
        left.saturation = 4;
        right.saturation = [4 6 8 10 12 14];
    else
        left.brightness = 6;
        right.brightness = [5 6];
        
        left.saturation = [4 6 8 10 12 14];
        right.saturation = 4;
    end
end
left.colorCombinations = allcomb(left.brightness, left.saturation);
right.colorCombinations = allcomb(right.brightness, right.saturation);
% left.stimulusSet = allcomb(1:numel(left.brightness), 1:numel(left.saturation));
% right.stimulusSet = allcomb(1:numel(right.brightness), 1:numel(right.saturation));

% Convert numbers into RGB from MUNSELL colours for each combination
% This is important for ensuring that the stimuli are numbered
% appropriately commensurate with the integral stimulus experiment
if condition == 1 && brightleft == 1
    cnt = 1;
    for i = 1:size(left.colorCombinations, 1)
        for j = 1:size(right.colorCombinations, 1)
            left.colorSet(cnt,:) = munsell2rgb('5R', left.colorCombinations(i,1), left.colorCombinations(i,2));
            right.colorSet(cnt,:) = munsell2rgb('5R', right.colorCombinations(j,1), right.colorCombinations(j,2));
            
            left.colorComb(cnt,:)  = left.colorCombinations(i,:);
            right.colorComb(cnt,:) = right.colorCombinations(j,:);
            
            cnt = cnt + 1;
        end
    end
elseif condition == 1 && brightleft == 0
    cnt = 1;
    for j = 1:size(right.colorCombinations, 1)
        for i = 1:size(left.colorCombinations, 1)
            left.colorSet(cnt,:) = munsell2rgb('5R', left.colorCombinations(i,1), left.colorCombinations(i,2));
            right.colorSet(cnt,:) = munsell2rgb('5R', right.colorCombinations(j,1), right.colorCombinations(j,2));
            
            left.colorComb(cnt,:)  = left.colorCombinations(i,:);
            right.colorComb(cnt,:) = right.colorCombinations(j,:);
            
            cnt = cnt + 1;
        end
    end
elseif condition == 2 && brightleft == 1
    cnt = 1;
    for j = 1:size(right.colorCombinations, 1)
        for i = 1:size(left.colorCombinations, 1)
            left.colorSet(cnt,:) = munsell2rgb('5R', left.colorCombinations(i,1), left.colorCombinations(i,2));
            right.colorSet(cnt,:) = munsell2rgb('5R', right.colorCombinations(j,1), right.colorCombinations(j,2));
            
            left.colorComb(cnt,:)  = left.colorCombinations(i,:);
            right.colorComb(cnt,:) = right.colorCombinations(j,:);
            
            cnt = cnt + 1;
        end
    end
elseif condition == 2 && brightleft == 0
    cnt = 1;
    for i = 1:size(left.colorCombinations, 1)
        for j = 1:size(right.colorCombinations, 1)
            left.colorSet(cnt,:) = munsell2rgb('5R', left.colorCombinations(i,1), left.colorCombinations(i,2));
            right.colorSet(cnt,:) = munsell2rgb('5R', right.colorCombinations(j,1), right.colorCombinations(j,2));
            
            left.colorComb(cnt,:)  = left.colorCombinations(i,:);
            right.colorComb(cnt,:) = right.colorCombinations(j,:);
            
            cnt = cnt + 1;
        end
    end
end
feedbackprobs(:,1) = [1 1 1 1 1 1 0 0 0 0 0 0]';
feedbackprobs(:,2) = [0 0 0 0 0 0 1 1 1 1 1 1]';


%% Set up color space and load into textures
% Load the stimuli into textures
nstimuli = size(left.colorSet, 1);
stimTexture = cell(nstimuli,1);

horizontalSeparation = 15; % Pixels
[stimloc, dh, dv] = CenterRect([0,0, 2 * stimsize + horizontalSeparation, stimsize], screenparms.rect);
left.loc  = [stimloc([1 2]) stimloc(1) + stimsize stimloc(4)];
right.loc = [stimloc(1) + stimsize + horizontalSeparation stimloc([2 3 4])];

[w, h] = RectSize(stimloc); % Width and Height of stimulus location area
stimmat = ones(h, w, 3) * screenparms.color(1); % Preload matrix
for i = 1:nstimuli % Create and load textures for the entire stimulus set
    for j = 1:3
        stimmat(:, 1:stimsize, j) = repmat(left.colorSet(i,j), stimsize, stimsize);
        stimmat(:, stimsize+horizontalSeparation+1:w, j) = repmat(right.colorSet(i,j), stimsize, stimsize);
    end
    stimTexture{i} = Screen('MakeTexture', screenparms.window, stimmat);
end
fixcrossH = CenterRect([0 0 fixsize 1], screenparms.rect);
fixcrossV = CenterRect([0 0 1 fixsize], screenparms.rect);

%% Setting up exposure schedule
% Exposure order based on session
sessionOrder = perms(3:-1:1); % 1 = Control, 2 = Correlated, 3 = Filtration
currentSession = sessionOrder(sessionOrder(session),:);
sessionSchedule = [randperm(2) + (currentSession(1)-1)*2, randperm(2)+ (currentSession(2)-1)*2 , randperm(2)+ (currentSession(3)-1)*2];

while find(sessionSchedule == 5) > find(sessionSchedule == 6); % if 6 is before 5, resample
   sessionSchedule = [randperm(2) + (currentSession(1)-1)*2, randperm(2)+ (currentSession(2)-1)*2 , randperm(2)+ (currentSession(3)-1)*2];
end

%% Present Instructions
instructionFolder = 'Instructions';
breakImage = (fullfile(pwd, instructionFolder, 'Break.bmp'));
endImage = (fullfile(pwd, instructionFolder, 'Thanks.bmp'));
ExpInstImage = (fullfile(pwd, instructionFolder, 'ExperimentInstructions.bmp'));

instructionSet = {'PracticeInstructions.bmp'};
for i = 1:size(instructionSet,2)
    showInstructions(screenparms, fullfile(pwd, 'Instructions', instructionSet{i})) %This corresponds to a function in the current directory
end

%% Run the Experiment
cnt = 1;
ConditionCount = 1;
overallcorrect = []; % Preallocate
totaloutput = [];
for bblocks = startBlock:nBlocks
    % Block Variables
    e = 0;
    nTrialsPerBlock = nExpTrials; % Trials in each block
    currentCondition = sessionSchedule(ceil(bblocks/2)); % Selects Correct Condition
    
    % Set up schedule for each block
    switch currentCondition
        case 1 % Control 1
            subStimNum = [1 3 5 7 9 11];
        case 2 % Control 2
            subStimNum = [2 4 6 8 10 12];
        case 3 % Correlated 1 - Positive
            subStimNum = [1 3 5 8 10 12];
        case 4 % Correlated 2 - Negative
            subStimNum = [2 4 6 7 9 11];
        case 5 % Filtration Condition 1
            subStimNum = [1:12];
        case 6 % Filtration Condition 2 % Identical to the previous
            subStimNum = [1:12];
        otherwise
            closeexp(screenparms)
    end
    
    % Practice Trials for entire first session
    if session == 1
        PracSession = 1;
        timeout = 9999;
    else
        PracSession = 0;
        timeout = 5;
    end
    PracCount = 0;

    % Practice Trials before each experimental block; Participants don't know this
    if mod((bblocks - 1),2) == 0
        nTrialsPerBlock = nPracTrials;
        PracCount = 1;
        if currentCondition == 6; % skips practice trials for second filtration condition
            nTrialsPerBlock = 0;
%             PracCount = 0;
        end;
    end % for if
   
    stimNo = size(subStimNum,2);
    
    output = [];
    
    % Set up Experimental Trials within this block
    stimInOrder = repmat(subStimNum', nTrialsPerBlock/stimNo, 1);
    currentStimulus = stimInOrder(randperm(numel(stimInOrder)));
    
    response = nan(nTrialsPerBlock,1); 
    rt = nan(nTrialsPerBlock,1);
    categoryFlag = nan(nTrialsPerBlock,1);
    correctFlag = nan(nTrialsPerBlock,1);
    
%     priorityLevel = MaxPriority(screenparms.window,'WaitBlanking');
    
    % Start Experiment
    trialcnt = 1;
    if bblocks == startBlock
        WaitSecs(1)
        showtext(screenparms, 20, 'Press any BUTTON to start EXPERIMENTAL trials', 1, 0, xleft);
        Screen('Flip',screenparms.window);
        RTBox('clear'); % clear buffer and sync clocks before stimulus onset
        while ~any(RTBox('ButtonDown')); WaitSecs(0.01); end % Wait for any button press
        WaitSecs(1);
        Screen('Flip',screenparms.window);
    elseif bblocks ~= startBlock;
        WaitSecs(1);       
%         blkcorrect = totaloutput(:, 15);
%         blkcorrect(isnan(blkcorrect),:) = [];
%         blkcorrect(blkcorrect == 9,:) = [];
%         nblk = numel(blkcorrect);
%         blkpercent= 100 * sum(blkcorrect)/nblk;
% 
%         showtext(screenparms, 20, sprintf('Mean Accuracy = %4.2f percent', blkpercent), 1, -60, xleft);
        showtext(screenparms, 20,'Take a BREAK. Press any BUTTON to CONTINUE', 1, 0, xleft);

        Screen('Flip',screenparms.window);
        RTBox('clear'); % clear buffer and sync clocks before stimulus onset
        while ~any(RTBox('ButtonDown')); WaitSecs(0.01); end % Wait for any button press
        WaitSecs(1);
        Screen('Flip',screenparms.window);
    end
    
    for i = 1:nTrialsPerBlock

        showtext(screenparms, 20, sprintf('%d To Go', nTrialsPerBlock - i + 1), 1, 0, xleft);
         Screen('Flip', screenparms.window); % Flip the fixation cross to the front buffer
        WaitSecs(fixationDuration/2); % Wait    
        
        % Display fixation cross
        Screen('DrawLine', screenparms.window, screenparms.black, fixcrossH(1), fixcrossH(2), fixcrossH(3), fixcrossH(4), lineWidth)
        Screen('DrawLine', screenparms.window, screenparms.black, fixcrossV(1), fixcrossV(2), fixcrossV(3), fixcrossV(4), lineWidth)
        Screen('Flip', screenparms.window); % Flip the fixation cross to the front buffer
        WaitSecs(fixationDuration/2); % Wait

        
%         Priority(priorityLevel);
        
        % Present Stimulus
        Screen('DrawTexture', screenparms.window, stimTexture{currentStimulus(i,1)}, [], stimloc);
        Screen('FrameRect', screenparms.window, [0 0 0], GrowRect(stimloc, 20, 20), 1)
        % Present Instructions
        if debug == 1;
%             fprintf('Stimulus No: %02d,   Left RGB %3d %3d %3d, Right RGB %3d %3d %3d\n',...
%                 currentStimulus(i,1), left.colorSet(currentStimulus(i,1),1), left.colorSet(currentStimulus(i,1),2), left.colorSet(currentStimulus(i,1),3),...
%                 right.colorSet(currentStimulus(i,1),1), right.colorSet(currentStimulus(i,1),2), right.colorSet(currentStimulus(i,1),3));

            fprintf('Stim: %02d, LBri = %3d, LSat = %3d, RBri = %3d, RSat = %3d\n',...
                currentStimulus(i,1), left.colorComb(currentStimulus(i,1),1), left.colorComb(currentStimulus(i,1),2),...
                right.colorComb(currentStimulus(i,1),1), right.colorComb(currentStimulus(i,1),2));
        end
        % Present instructions (last two inputs manipulate location: first
        % is vertical offset (negative is above midline, positive is below
        % midline), second is horizontal offset (negative is left of
        % midline, positive is right of midline)
        showtext(screenparms, textsize/3, 'Which category does this object belong to?', 1, 250, 0);
        showtext(screenparms, textsize/3, '(CATEGORY A)',                                    1, 350, -screenparms.rect(4)/4);
        showtext(screenparms, textsize/3, '(CATEGORY B)',                                    1, 350,  screenparms.rect(4)/4);
        
        % Clear RT box buffer and sync clocks before stimulus onset
        RTBox('clear');
        vbl = Screen('Flip', screenparms.window);
        
        % Record response
        [cpuTime, buttonPress] = RTBox(timeout);  % computer time of button response
%         Priority(0);
        
        % Clear the screen
        FillScreen(screenparms);
        Screen('Flip', screenparms.window);
        
        if ~isempty(cpuTime)
            rt(i, 1)   = (cpuTime - vbl);
        end
        
        if ismember(buttonPress, {'1', '2'})
            response(i, 1) = 1;
        elseif ismember(buttonPress, {'3', '4'})
            response(i, 1) = 2;
        end
        
        % Display Feedback
        if feedbackprobs(currentStimulus(i,1),1);
            categoryFlag(i,1) = 1; % category is A
        else
            categoryFlag(i,1) = 2; % category is B
        end
        
        if any([response(i,1) == categoryFlag(i,1) && rt(i,1) < timeout && PracSession == 0,...
                response(i,1) == categoryFlag(i,1) && PracSession == 1])
            correctFlag(i,1) = 1;
            if PracSession == 1
                % Present stimulus
                Screen('DrawTexture', screenparms.window, stimTexture{currentStimulus(i,1)});
                Screen('FrameRect', screenparms.window, [0 0 0], GrowRect(stimloc, 20, 20), 1)

                %Show the feedback
                showtext(screenparms, textsize, feedback{2}, 1, fbackLocation(1), fbackLocation(2));
                Screen('Flip', screenparms.window);
                WaitSecs(fbkDuration);
            end
            FillScreen(screenparms);
            Screen('Flip', screenparms.window);

        elseif any([response(i,1) ~= categoryFlag(i,1) && rt(i,1) < timeout && PracSession == 0,...
                response(i,1) ~= categoryFlag(i,1) && PracSession == 1])
            correctFlag(i,1) = 0;
            
            % Present Stimulus
            Screen('DrawTexture', screenparms.window, stimTexture{currentStimulus(i,1)});
            Screen('FrameRect', screenparms.window, [0 0 0], GrowRect(stimloc, 20, 20), 1)

            showtext(screenparms, textsize, feedback{1}, 1, fbackLocation(1), fbackLocation(2));
            if PracSession == 1
                Screen('Flip', screenparms.window);
                WaitSecs(fbkDuration);
              
            else
                Screen('Flip', screenparms.window);
                WaitSecs(fbkDuration);
                
            end % for session..
            
            FillScreen(screenparms);
            Screen('Flip', screenparms.window);
            
        else
            correctFlag(i,1) = 9;
            showtext(screenparms, textsize, 'Too Slow!', 1, 0, 0);
            
            Screen('Flip', screenparms.window);
            WaitSecs(fbkDuration);
            
            FillScreen(screenparms);
            Screen('Flip', screenparms.window);
        end
        
        trialcnt = trialcnt + 1;
        cnt = cnt + 1;
    end
 
    output = [(1:nTrialsPerBlock)', currentStimulus,...
        left.colorComb(currentStimulus,1), left.colorComb(currentStimulus,2),...
        right.colorComb(currentStimulus,1), right.colorComb(currentStimulus,2),...
        categoryFlag, response, correctFlag, rt];
    
    %% Save data
    finaloutput = [...
        repmat(subject, nTrialsPerBlock,1),...         % Subject
        repmat(condition,nTrialsPerBlock,1),...        % Condition (1 = Brightness, 2 = Saturation)
        repmat(brightleft,nTrialsPerBlock,1),...       % Brightness on left, or Brightness on right
        repmat(session, nTrialsPerBlock,1),...         % Session 
        repmat(currentCondition,nTrialsPerBlock,1),... % Task
        repmat(PracCount,nTrialsPerBlock,1),...        % Practice?
        output];                                       % Output
    totaloutput = [totaloutput; finaloutput];
    
    dlmwrite(tmpfile,totaloutput,'-append');
    
    %% Prepare end of day summary
    %% CHANGED
        if bblocks == nBlocks
        
        %Output Warning Screen
        WaitSecs(fbkDuration);
        showtext(screenparms, 20, 'Preparing Output', 1, 0, 0);
        Screen('Flip', screenparms.window);
        
        dlmwrite(outputfile, totaloutput, '-append');
        
        finalData = totaloutput;
        
        finalData = totaloutput(nPracTrials+1:end,:);
        
        accuracyrate   = aggregate(finalData, 5, 10);               % Find the p(correct) in this block by each item
        nitemsperblock = aggregate(finalData, 5, 10, @count, 1);    % Find the total N per block for each item
        errorsperblock = nitemsperblock .* (1 - accuracyrate(:,2));  % Compute the number of errors in each block for each item
        
        fdata = finalData;
        fdata(fdata(:,10) == 0,:) = [];             % Now remove errors
        minrt = .2;
        means = aggregate(fdata, 5, 11);        % Compute the means for each item
        stds  = aggregate(fdata, 5, 11, @std);  % Compute the stds for each item
        
        % Iterate through the items and remove trials with RT's > 3 * std + mean or < .2
        for i = 1:size(means,1);
            outliers{i} = find(all([fdata(:,5) == i, any([fdata(:,11) < minrt, fdata(:,11) > means(i,2) + stds(i,2) * 3], 2)], 2));
            fdata(outliers{i}, :) = [];
        end
        mx = aggregate(fdata, 5, 11);           % Mean RT for each item
        nx = aggregate(fdata, 5, 11, @count);   % N for each item
        
        % TEMP - This is for testing only, at full trial numbers this should never happen
        if size(mx,1) ~= stimNo;  mx = [mx; ones(stimNo-size(mx,1),2)]; nx = [nx; ones(stimNo-size(nx,1),2)]; end
        if size(accuracyrate,1) ~= 12; accuracyrate =  [accuracyrate; ones(stimNo-size(accuracyrate,1),2)]; end
        
        finalsummary = [(1:stimNo); accuracyrate(:,2)'; mx(:,2)' * 1000; nx(:,2)'];
        finalsummary(2:3,:) = round(finalsummary(2:3,:) * 1000)/1000;
        fins{1} = ['End of day summary for subject number ' num2str(subject)];
        fins{2} = sprintf('%8d %8d %8d %8d %8d %8d %8d %8d %8d', finalsummary(1,:));
        fins{3} = sprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f', finalsummary(2,:));
        fins{4} = sprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f', finalsummary(3,:));
        fins{5} = sprintf('%8d %8d %8d %8d %8d %8d %8d %8d %8d', finalsummary(4,:));
        
        edata = finalData(finalData(:,5) ~= 1,:);
        bonusscore = mean(edata(:,12));
        bonusNtrials = size(edata(:,12), 1);
        
        totalcorrect=totaloutput(:,15);
        totalcorrect(isnan(totalcorrect),:) = [];
        totalcorrect(totalcorrect == 9,:) = [];
        ntotal = numel(totalcorrect);
        totalpercent= 100 * sum(totalcorrect)/ntotal;
           
        if totalpercent > 98
            fins{6} = sprintf('Bonus = Yes, Overall Accuracy= %4.2f percent', totalpercent);
        else
            fins{6} = sprintf('Bonus = No, Overall Accuracy = %4.2f percent', totalpercent);
        end
        disp(fins{6})
       
    end
end

%% Quitting the experiment
leftoffset = -100;
showtext(screenparms, 20, 'Thank you for participating. Please call experimenter.', 1, 200, leftoffset);
showtext(screenparms, 20, fins{6}, 1, -200, leftoffset);
Screen('Flip', screenparms.window);
pressEnter

% WaitSecs(2);
% closeexp(screenparms) %% Close the Experiment

%%if crashed, run these lines
%crashdata = importdata(tmpfile);
%crash_totalacc = 100 * sum(crashdata(:, 9))/ numel(crashdata(:, 9))