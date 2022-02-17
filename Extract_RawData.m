%% Clean Workspace
clear;
close all;
clc;
%% Loop all participants
who = "amp";
rest = 0;                           % 1: include resting / 0: discard resting
numCommGes = 4;                     % Number of Common Gestures
GestList = ["Flexion" "Extension" "Pronation" "Supination"];    % Mapping Gestures Table
if who == "able"
    from = 1;    fin = 10;    fl = 3; 	titl = "Able-bodied";    numGestures = 4 + rest;
else
    from = 11;   fin = 13;   fl = 3;  	titl = "Amputees";       numGestures = 4 + rest;
end
Raw_Data = cell(1,fin-from+1);
%% DATA PROCESSING - Loop through participants
for part=from:fin
    %% General variables
    partNum = part;                    	% Participant ID number
    partID = int2str(partNum);          % Participant ID as string
    numberTrials = 5;                   % Number of Trials
    musclepairs = 6;
    %% Data Processing variables
    fs = 1000;                          % sample frequency
    hp = 1;                             % high pass frequency
    lp = 150;                           % low pass frequency
    Gestures = cell(3,fl*numGestures*numberTrials);
    %% Load marked labels file
    load(strcat('../Participants/Par',partID,'/Marks',partID,'.mat'));
    g_off = 0;
    count = 0;
    %% Process Participant Data
    for file = 1:2
        try
            M = csvread([strcat('..\..\Data Library\ParticipantData\F',partID) Marks{file,1}(1,5:end)]);
        catch
            M = csvread([strcat('..\..\Data Library\ParticipantData\P',partID) Marks{file,1}(1,5:end)]);
        end
        M(1,36) = 0;
        M = M(:,20:36);
        [x,y] = size(M);
        N = zeros(x,(y-1)/2+1);
        for p = 1:x
            for q = 1:8
                N(p,q) = (M(p,((q*2)-1))*255)+(M(p,((q*2))));
            end
            N(p,9) = p/fs;
        end
        %% Data Preprocessing
        progest = N(:,1:8);
        %% Data Segmentation Variables
        offset = 1000;                              % Offset from marks
        offset1 = 750;
        windowLength = 3000-1;                      % Lenght of data to extract
        %% Segment Location Variables
        TotTrials = 5;                              % Total number of trials
        numchannels = 4;                            % Number of active channels
        fileID = 5*count+(file-1)*numberTrials*2;   % Column position of next trial   
        g_off = 0;
        hol = zeros(windowLength+1,numchannels);
        if (mod(file,((numGestures-1)/2))==1 && rest==1)
            ge = 2;
            %% Gesture 0 - Rest
            for g = 1:numberTrials
                hol(:,:) = progest(Marks{file,2}(1,ge)+offset1:Marks{file,2}(1,ge)+offset1+windowLength,1:numchannels);
                Gestures(:,g+fileID) = {hol() Marks{file,4} 0};
                ge = ge + 2;
            end
            g_off = 5;
            count = count+1;
        end
        ge = 1;                                     % Starting position of marker
        %% 1-5 Trials
        for g = 1:numberTrials
            hol(:,:) = progest(Marks{file,2}(1,ge)+offset:Marks{file,2}(1,ge)+offset+windowLength,1:numchannels);
            Gestures(:,g+fileID+g_off) = {hol() Marks{file,4} Marks{file,5}};
            ge = ge + 2;
        end
        ge = ge + (TotTrials-numberTrials)*2;
        %% 6-10 Trials
        for g = (numberTrials+1):numberTrials*2
            hol(:,:) = progest(Marks{file,2}(1,ge)+offset:Marks{file,2}(1,ge)+offset+windowLength,1:numchannels);
            Gestures(:,g+fileID+g_off) = {hol() Marks{file,4} Marks{file,6}};
            ge = ge + 2;
        end
    end
    keep = any(~cellfun('isempty',Gestures), 1);
    Gestures = Gestures(:,keep);
    Raw_Data{part} = Gestures;
end
keep = any(~cellfun('isempty',Raw_Data), 1);
Raw_Data = Raw_Data(:,keep);