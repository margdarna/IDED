%----------------------------------------------------------------------
%                       Description
%----------------------------------------------------------------------
%
% This experiment encompasses an adapted version of the Intra- and Extra-dimensional
% Shifting Task from the CANTAB Battery
% This version of the ASST was created based on the following paper:
% Oh, A., Vidal, J., Taylor, M. J., & Pang, E. W. 
% (2014). Neuromagnetic
% correlates of intra- and extra-dimensional set-shifting. Brain and Cognition,
% 86, 90-97. https:\\doi.org/10.1016/j.bandc.2014.02.006
%
% Version 2:
% 1. Added is_break in all variables that get saved
% 2. variables get saved in all_variables
% 3. Changed \\ to /
% 4. changed end of practice text to two separate screens. This way the
% participant cannot accidentally start with the actual experiment without
% getting the chance of asking questions
% 5. Added option to choose subject when performance mat files are missing
% but a subject folder already exists
% 6. Change event locking after response: 30 + feedback where 1: correct
% response, 2: incorrect/miss
%
% Version 3
% 1. fixed minor bugs
% 2. Added event to signal last trial before switch
%
% Version 4
% 1. changed / to \\ for compatibility in MATLAB in EEG LAB
% 2. changed input to mouse (from keyboard)
% 3. change instructions to fit the new input type
% 4. made fixation cross bigger
% 5. added line to stop synchronization test in case that exp does not run
%
% Version 5 and 6
% 1. adaptations made for new lab
% 2. added mouse as input device
%
% Version 7
% 1. adapted triggers to work in new lab
%
% Version 8
% 1. Changed minor details so that triggers are sent correctly
%
% COMMENT:
% triggers are right now commented oud. They only work on a PC with
% signalling port set up. To enable triggers with io64 remove the comment
% from all trig lines.
%
% Written by: Margarita Darna
% margarita.darna@lin-magdeburg.de
% Leibniz Institute for Neurobiology, Magdeburg
%
% Created:     24.11.2021
% Last edited: 01.09.2023

%%
%%----------------------------------------------------------------------
%                       Preamble
%----------------------------------------------------------------------
% Clear the workspace
close all;
clear all;
sca;

% Define directories
% from EEG Lab
project_dir = 'C:\your_project_directory\';
% remaining directories
exp_dir = sprintf('%sExpCode\\', project_dir);
img_dir = sprintf('%sStimuli\\', exp_dir);
raw_data_dir = sprintf('%sData\\Raw_Data\\', project_dir);
subj_stim_dir = sprintf('%sSubject_Stimuli\\', exp_dir);

% io64 directory, directory where mexfile is located
io64_dir = sprintf('%sExpCode\\\trigger\\', project_dir);

%Screen('Preference', 'SkipSyncTests', 1);%
%----------------------------------------------------------------------
%                       Subject Information
%----------------------------------------------------------------------

% Ask for subject number
prompt = {'Enter subject number(Must be a number between 001 and 500):'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'500'}; % default subject number for testing
% all numbers higher than 400 can be used for testing
subject = inputdlg(prompt,dlgtitle,dims,definput);
subject = subject{1,1};

% check if subject number already exists
if strcmp(subject, '500') == 0
    exist = dir(sprintf('%sS%s', raw_data_dir, subject));    
    if isempty(exist) == 0
        mat_exist = dir(sprintf('%sS%s\\Behavior\\*_task-IDED_all_var.mat', raw_data_dir, subject));
        if isempty(mat_exist) == 0
            prompt = {sprintf('Performance file of subject number %s already exists. Enter new subject number:', subject)};
            dlgtitle = 'Input';
            dims = [1 35];
            subject = inputdlg(prompt,dlgtitle,dims);
            subject = subject{1,1};
            mat_exist = dir(sprintf('%sS%s\\Behavior\\*all_var.mat', raw_data_dir, subject));
            exist = dir(sprintf('%sS%s', raw_data_dir, subject));
        else
            dlgTitle    = 'User Question';
            dlgQuestion = 'No mat files for this task where found. However, a folder for this subject already exists. Continue with experiment?';
            choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
            dims = [1 35];
            if strcmp(choice, 'No')
                fprintf('Experiment stopped by user!\n')
                return
            end
        end
    end
end

% Create folder to save subject information into
if not(isfolder(sprintf('%sS%s', raw_data_dir, subject)))
    mkdir(sprintf('%sS%s', raw_data_dir, subject))
    mkdir(sprintf('%sS%s\\Behavior', raw_data_dir, subject))
    mkdir(sprintf('%sS%s\\EEG', raw_data_dir, subject))
    mkdir(sprintf('%sS%s\\Exp_Files', raw_data_dir, subject))
end

% Load mat file with stimulus order for that specific participant
try
    load(sprintf('%sS%s_stimuli.mat', subj_stim_dir, subject));
    % Paste Experiment File in Raw Data folder for safe keeping in case the
    % originals accidentally get changed
    copyfile(sprintf('%sS%s_stimuli.mat', subj_stim_dir, subject), sprintf('%sS%s\\Exp_Files', raw_data_dir, subject))
catch
    fprintf('Could not find stimuli mat file for subject %s\n', subject);
    return;
end

%%
%----------------------------------------------------------------------
%                       IO64
%----------------------------------------------------------------------
% initialize port (port adress C0B0)

% comment these lines out when the code is ran on a PC without port
% signalling. Otherwise code does not run
% addpath(io64_dir);
% trig = triggerhandler();
% trig.send_trigger(0);


% Port signalling:
% Fixation Cross:            10
% Stimulus Presentation:
%             practice:      210 + Target position
%                  pre:      220 + Target position
%               repeat:      230 + Target position
%                   ID:      240 + Target position
%                   ED:      250 + Target position
%                 last:      290 + Target position   
% Button press:              30  + correctness
% Feedback(practice trials): 40  + Feedback
% Break start:               50  + number of break
% Break end:                 60  + number of break

%%
%----------------------------------------------------------------------
%                       Psychtoolbox
%----------------------------------------------------------------------

% Setup PTB with some default values
PsychDefaultSetup(2);

%rand('seed', sum(100 * clock)); % for older syntax
rng('shuffle')

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Flip to clearewe
Screen('Flip', window);

% Set the text size
Screen('TextSize', window, 60);

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set the blend function for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set language settings for Psychtoolbox
Screen('Preference', 'TextEncodingLocale', 'UTF-8');

%%
%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Query the frame duration
ifi = Screen('GetFlipInterval', window);
slack = ifi/2;

%%
%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

% Define the keyboard keys that are listened for. We will be using the escape key as
% a exit/reset key

escapeKey = KbName('ESCAPE');

% in case we use keyboard as input device instead of mouse
%leftKey = KbName('LeftArrow');
%rightKey = KbName('RightArrow');

%[nclick,xmouse,ymouse, button] = GetClicks(window);

%%
%----------------------------------------------------------------------
%                        Fixation Cross
%----------------------------------------------------------------------

% Here we set the size of the arms
fixCrossDimPix = 20;

% Set the line width for our fixation cross
lineWidthPix = 3;

% Now we set the coordinates (relative to zero)
xCoords = [-fixCrossDimPix fixCrossDimPix fixCrossDimPix -fixCrossDimPix];
yCoords = [-fixCrossDimPix fixCrossDimPix -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

%%
%----------------------------------------------------------------------
%                        Feedback
%----------------------------------------------------------------------
% set feedback size
dotSizePix = 30;

%%
%----------------------------------------------------------------------
%                           Text
%----------------------------------------------------------------------
welcome_text1 = 'Wilkommen!\n';
welcome_text2 = sprintf('Gleich beginnen wir mit den Übungsdurchgängen.\n ');

finish_practice_text1 = 'Ende der Übungsrunde.\n';
finish_practice_text2 = 'Gibt es noch Fragen?\n';

break_text   = 'PAUSE\n\n';
press_any_button_text = 'Drücke eine beliebige Maustaste zum Weitermachen.\n';

end_text = 'ENDE\nVielen Dank für die Teilnahme.\nBitte wende dich an die Versuchsleiterin';

%%
%----------------------------------------------------------------------
%                           Stimuli
%----------------------------------------------------------------------
% We will set the height of each drawn image to a fraction of the screens
% height

scale_factor = 0.1;
imageHeight = screenYpixels .* scale_factor;
imageWidth = imageHeight; % Stimuli have the same width and height

% Creating Scaling Rectangle
y_pos    = 100;
x_pos    = 200;
theRect  = [0 0 imageWidth imageHeight];
dstRectC = CenterRectOnPointd(theRect, xCenter, yCenter + 2 * y_pos);
dstRectL = CenterRectOnPointd(theRect, xCenter - x_pos, yCenter - y_pos);
dstRectR = CenterRectOnPointd(theRect, xCenter + x_pos, yCenter - y_pos);

% Stimuli position, make sure that the correct stimulus appear on the left
% and right side an equal amount of times
num_pos = round(numel(trial_num)/2);
pos_Tar = [ones(num_pos,1); repmat(2,num_pos,1)]; % left = 1, right = 1

%%
%----------------------------------------------------------------------
%                     Make a response matrix
%----------------------------------------------------------------------

% This a matrix that records the responses of participants
% In columns:
% 1: trial number
% 2: Target position (1: left, 2: right)
% 3: respond key (1: left, 2: right)
% 4: reaction time
respMat = nan(numel(trial_num), 4);

%%
%----------------------------------------------------------------------
%                     Make a timing log matrix
%----------------------------------------------------------------------

% This a matrix that logs the timing information of the experiment
% In columns:
% 1: onset of fixation cross
% 2: onset of stimulus
% 3: Press of Button
% 4: onset of black screen after response
t_Mat = nan(numel(trial_num), 4);

%%
%----------------------------------------------------------------------
%                         Welcome Screen
%----------------------------------------------------------------------

% Draw all the text in one go
Screen('TextSize', window, 70);
DrawFormattedText(window, [welcome_text1 welcome_text2], ...
    'center', yCenter, white);

% Flip to the screen
HideCursor;
Screen('Flip', window);

% Now we wait for any keyboard button press to move to the practice trials
KbStrokeWait;

%%
%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------
try    
    % shuffle rows pos_Tar
    sh_pos_Tar = Shuffle(pos_Tar, 2);
    % set break number as zero;
    break_num = 0;
    for trial = 1:numel(trial_num) 
        % Single Trial
        % Prepare Images
        try
            [Target, ~, alpha] = imread(sprintf('%s%s%s.png', img_dir, target_col(trial), target_sh(trial)));
            Target(:, :, 4) = alpha;
            TargetTexture = Screen('MakeTexture', window, Target);
            [Match, ~, alpha] = imread(sprintf('%s%s%s.png', img_dir, match_col(trial), match_sh(trial)));
            Match(:, :, 4) = alpha;
            MatchTexture = Screen('MakeTexture', window, Match);
            [NoMatch, ~, alpha] = imread(sprintf('%s%s%s.png', img_dir, no_match_col(trial), no_match_sh(trial)));
            NoMatch(:, :, 4) = alpha;
            NoMatchTexture = Screen('MakeTexture', window, NoMatch);
        catch
            fprintf('Could not load images\n');
            sca;
            return;
        end
        % Present Fixation Cross
        Screen('DrawLines', window, allCoords, lineWidthPix, grey, ...
            [xCenter yCenter], 2);

        if trial == 1 || trial == practice_trials + 1 || is_break(trial - 1)
            t_fixation = Screen('Flip', window);
        elseif trial_type(trial) == 'practice'
            t_fixation = Screen('Flip', window, t_feedback + 1 - slack );
        else
            t_fixation = Screen('Flip', window, t_black + 1 - slack );
        end

        % port signal for fixation cross
        %trig.send_trigger(10);

        % Draw Images
        Screen('DrawLines', window, allCoords, lineWidthPix, grey, ...
            [xCenter yCenter], 2);
        Screen('DrawTextures', window, TargetTexture, [], dstRectC);
        if sh_pos_Tar(trial) == 1
            Screen('DrawTextures', window, MatchTexture, [], dstRectL);
            Screen('DrawTextures', window, NoMatchTexture, [], dstRectR);
        elseif sh_pos_Tar(trial) == 2
            Screen('DrawTextures', window, MatchTexture, [], dstRectR);
            Screen('DrawTextures', window, NoMatchTexture, [], dstRectL);
        else
            fprintf("Target position not found. Check sh_pos_Tar\n")
            sca;
            return;
        end
        
        % preparing port signal for stimulus presentation
         if trial_type(trial) == 'practice'
             signal = 10;
         elseif trial_type(trial) == 'pre' 
             signal = 20;
         elseif is_last(trial) 
             signal = 90;
         elseif trial_type(trial) == 'repeat'
             signal = 30;
         elseif trial_type(trial) == 'ID'
             signal = 40;
         elseif trial_type(trial) == 'ED'
             signal = 50;
         else
             fprintf('Trial type not found. Please check trial_type\n')
             sca;
             return
         end

        % Flip to the screen
        t_stimuli = Screen('Flip', window, t_fixation + 0.800 + rand(1)/2 - slack);
        
        % port signal
        % fprintf('%d\n', signal);
        %trig.send_trigger(200 + signal + sh_pos_Tar(trial));

        % wait for press of button
        % Cue to determine whether a response has been made
        respToBeMade = true;
        tStart = GetSecs;
        while respToBeMade == true
            tNow = GetSecs;
            % Check the mouse. The person should press the
            [~,~, button] = GetMouse(window);
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                fprintf('Experiment abrupted through press of Escape key\n');
                abrupted = true;
                last_stage = stage_names{k};
                sca;
                return;
            elseif button(1)
                response = 1;
                respToBeMade = false;
            elseif button(3)
                response = 2;
                respToBeMade = false;
            elseif tNow - tStart > 4
                response = 3; % miss
                respToBeMade = false;
            end
        end

        if response == sh_pos_Tar(trial)
                feedback = 1;       % correct
        else
                feedback = 0;       % incorrect
        end

        % port signal for response
        %trig.send_trigger(30 + feedback);

        tEnd = GetSecs;
        rt = tEnd - tStart;

        % switching to black screen
        Screen('FillRect', window, [0 0 0]);
        t_black = Screen('Flip', window);

        % Give feedback if trials are practice trials
        if trial_type(trial) == 'practice'
            % Give feedback
            if response == sh_pos_Tar(trial)
                dotColor = [0 1 0]; % green
            else
                dotColor = [1 0 0]; % red
            end
            Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
            t_feedback = Screen('Flip', window, t_black + 0.800 + rand(1)/5 - slack);
            % port signal for feedback
            %trig.send_trigger(40 + feedback);
        end
        
        % Logging trial information
        respMat(trial, 1) = trial;
        respMat(trial, 2) = sh_pos_Tar(trial);
        respMat(trial, 3) = response;
        respMat(trial, 4) = rt;

        % Logging timing information
        t_Mat(trial, 1) = t_fixation;
        t_Mat(trial, 2) = t_stimuli;
        t_Mat(trial, 3) = tEnd;
        t_Mat(trial, 4) = t_black;

        % save data for backup
        % saving subject preformance
        save(sprintf('%sS%s\\Behavior\\S%s_task-IDED_all_var.mat', raw_data_dir, subject, subject));

        % Show end of practice session if trial is the last practice trial
        if trial == practice_trials
            % Presenting Text to finish practice
            Screen('TextSize', window, 50);
            DrawFormattedText(window, [finish_practice_text1 finish_practice_text2], ...
                'center', yCenter, white);
            Screen('Flip', window, t_feedback + 1 - slack);
            KbStrokeWait;
        end

        % Show break screen is the trial is set as a practice trial
        if is_break(trial)
            % change break num
            break_num = break_num + 1;
            % Present break text
            Screen('TextSize', window, 50);
            DrawFormattedText(window, break_text, 'center', yCenter, white);
            t_break = Screen('Flip', window, t_black + 1 - slack);
            % port signal for start of break
            %trig.send_trigger(50 + break_num);

            % Change break text after 2 seconds
            Screen('TextSize', window, 50);
            DrawFormattedText(window, [break_text press_any_button_text], 'center', yCenter, white);
            Screen('Flip', window, t_break + 2 - slack);

            respToBeMade = true;
            while respToBeMade
                [~,~, button] = GetMouse(window);
                if any(button)
                    response = 1;
                    respToBeMade = false;
                end
            end
            % port signal for end of break
            %trig.send_trigger(60 + break_num);
        end

    end
     % Present end screen
    Screen('TextSize', window, 70);
    DrawFormattedText(window, end_text,'center', yCenter, white);
    Screen('Flip', window);
    
    try
        % save subject performance and all variables
        save(sprintf('%sS%s\\Behavior\\S%s_task-IDED_all_var.mat', raw_data_dir, subject, subject));
        % saving as .csv in a MATLAB 2018 compatible way
        csvwrite(sprintf('%sS%s\\Behavior\\S%s_task-IDED_respMat.csv', raw_data_dir, subject, subject),respMat);
        csvwrite(sprintf('%sS%s\\Behavior\\S%s_task-IDED_tMat.csv', raw_data_dir, subject, subject),t_Mat);
        dlmwrite(sprintf('%sS%s\\Behavior\\S%s_task-IDED_respMat.tsv', raw_data_dir, subject, subject), respMat, 'delimiter', '\t');
        dlmwrite(sprintf('%sS%s\\Behavior\\S%s_task-IDED_tMat.tsv', raw_data_dir, subject, subject), t_Mat, 'delimiter', '\t');
    catch
        fprintf('Participant Performance Data could not be saved correctly')
        sca;
        return;
    end
    % Wait for press of button to close the screen
    KbStrokeWait;

    % Clear the screen
    sca;
catch
    sca
    fprintf('An error occured. Please check your code in the experimental loop!\n');
    return;
end
