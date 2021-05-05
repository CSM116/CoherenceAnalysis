%% Clean Workspace
clear;
close all;
clc;
%% Loop all participants
% whos = ["able","amp"];
who = "able";
if who == "able"
    from = 1;    fin = 7;    fl = 3;    titl = "Able-bodied";    numGestures = 7;
else
    from = 11;   fin = 13;   fl = 3;    titl = "Amputees";       numGestures = 5;
end
coh_part = cell(1,fin-from+1);       	% Coherence cell array of Participants
%% Mapping Gestures Table
if from >= 11 && from <= 20
    GestList = ["Rest" "Flexion" "Extension" "Pronation" "Supination"];
else
    GestList = ["Rest" "Flexion" "Extension" "Pronation" "Supination" "Adduction" "Abduction"];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Coherence Results / Process Data
try
    load(strcat(titl,'-coh.mat'));
    load(strcat(titl,'-nf.mat'));
catch
    %% DATA PROCESSING - Loop through participants
    for part=from:fin
        %% General variables
        partNum = part;                    	% Participant ID number
        partID = int2str(partNum);          % Participant ID as string
        numberTrials = 5;                   % Number of Trials
        musclepairs = 6;
        %% Data Processing variables
        fs = 1000;                          % sample frequency
        hp = 0.5;                           % high pass frequency
        lp = 150;                           % low pass frequency
%         num_hz = 200;
        order = 4;                          % order
        Gestures = cell(3,fl*numGestures*numberTrials);
        %% Load marked labels file
        load(strcat('../Participants/Par',partID,'/Marks',partID,'.mat'));
        g_off = 0;
        count = 0;
        %% Process Participant Data
        for file = 1:size(Marks,1)
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
            %{
    %             [b2, a2] = butter(order, [(2*hp)/fs (2*lp)/fs], 'bandpass');
    %             filtered_data = filtfilt(b2, a2, R);
    %             hilbert_data = hilbert(filtered_data);
    %             progest = dataProcessor(R(:,1:8),fs,hp,lp,order);
    %             progest = abs(progest);
    %             progest = filtered_data;
            %}
            %% Data Segmentation Variables
            offset = 1000;                          % Offset from marks
            offset1 = 750;
            windowLength = 3000-1;                  % Lenght of data to extract
            %% Segment Location Variables
            TotTrials = 5;                          % Total number of trials
            numchannels = 4;                        % Number of active channels
            fileID = 5*count+(file-1)*numberTrials*2; % Column position of next trial   
            g_off = 0;
            hol=zeros(windowLength+1,numchannels);
            if (mod(file,((numGestures-1)/2))==1)
                ge = 2;
                %% Gesture 0 - Rest
                for g=1:numberTrials
                    hol(:,:) = progest(Marks{file,2}(1,ge)+offset1:Marks{file,2}(1,ge)+offset1+windowLength,1:numchannels);
                    Gestures(:,g+fileID) = {hol() Marks{file,4} 0};
                    ge = ge + 2;
                end
                g_off = 5;
                count = count+1;
            end
            ge = 1;                                 % Starting position of marker
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

        %% Separate Trials by Gestures and Force Level
        gestF = cell(fl,numGestures);       % Rows: Force Level - Columns: Gestures
        [b2, a2] = butter(order, [(2*hp)/fs (2*lp)/fs], 'bandpass');
        for m=1:fl
            for w=1:numGestures
                for i=1:numberTrials
                    a = Gestures{1,i+numberTrials*(w-1)+numGestures*numberTrials*(m-1)};
                    filt_data = filtfilt(b2, a2, a);
                    hilb_data = hilbert(filt_data);
    %                 hilb_data = abs(hilb_data);
                    gestF{m,w}(:,:,i) = hilb_data;
                end
            end
        end
    %{
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %% Wavelet convolution
        %{

        sig = gestF{1,1}(:,:,:);
        % frequency parameters
        frex = linspace(hp,lp,num_hz);

        % other wavelet parameters
        range_cycles = [ 4 10 ];
    %     s = logspace(log10(.6),log10(.3),num_hz);
        s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_hz) ./ (2*pi*frex);
        wavtime = -2:1/fs:2;
        half_wave = (length(wavtime)-1)/2;

        % FFT parameters
        nWave = length(wavtime);
        nData = size(sig,1)*size(sig,3);
        nConv = nWave + nData - 1;

        % now compute the FFT of all trials concatenated
        alldata = reshape(gestF{1,1}(:,1,:),[],1,1);
        dataX   = fft(alldata, nConv);

        % initialize output time-frequency data
        tf = zeros(size(sig,1),num_hz);

        %% Perform convolution - loop over frequencies
        for fi=1:length(frex)
            % create wavelet and get its FFT
            wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
            waveletX = fft(wavelet,nConv);
            waveletX = waveletX ./ max(waveletX);

            % now run convolution in one step
            as = ifft(waveletX.' .* dataX);
            as = as(half_wave+1:end-half_wave);

            % and reshape back to time X trials
            as = reshape( as, [], size(sig,3) );

            % compute power and average over trials
            tf(:,fi) = mean( abs(as).^2 ,2);
        end

        %% Plot Heatmap
        figure;
        contourf(1:size(sig,1),frex,tf.',40,'linecolor','none');
        % set(gca,'clim',[1 4],'xlim',[-500 1300])
        colormap hot;
        xlabel('Time (ms)'), ylabel('Frequency (Hz)'); 
        %}

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% FFT - amplitude and power of signal
         %{
        %% Plot amplitude and power of signal
        signal = reshape(gestF{1,1}(:,:,:),[],4,1);
        signal1 = signal(:,1);
        signal2 = signal(:,2);

        % variables
        N = length(real(signal1));             	% length of sequence
        time  = (0:N-1)/fs;
        nyquist     = fs/2;                     % Nyquist frequency
        hz = linspace(0,nyquist,floor(N/2)+1);  % Frequencies in Hz

        % Plot the time-domain signal
        figure;
        plot(time,real(signal1),'b');hold on;
        % Detrended signal
        d_signal1 = detrend(signal1,1,3001);
        plot(time,real(d_signal1),'--r');
        % Downsampling
        down_sig = downsample(d_signal1,4);
        plot(downsample(time,4),real(down_sig),'k','linew',1.1);hold off;
        xlabel('Time (s)'), ylabel('Voltage (\muV)');
        legend({'signal1','detrended','downsampled'});

        % Power and amplitude
        figure;
        for i = 1:numchannels
            subplot(2,numchannels/2,i);
            ampl = abs(fft(signal(:,i))/N);
    %         ampl = abs(signal(:,i));
            powr = ampl.^2;

            plot(hz,ampl(1:length(hz)),'--k','linew',1.1);hold on;
            plot(hz,powr(1:length(hz)),'r','linew',1);hold off;
            set(gca,'xlim',[0 30])

            title(strcat('Channel-',int2str(i)));
            xlabel('Frequency (Hz)');
            ylabel('Amplitude or power');
            legend({'Amplitude';'Power'});
        end
         %}

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %% Coherence & Welch's method animation
         %{
        %% Signal declaration
        signal = reshape(gestF{1,1}(:,:,:),[],4,1);
        signal1 = signal(:,3);
        signal2 = signal(:,4);  
        %% Parameters
        N = length(real(signal1));             	% length of sequence
        nyquist	= fs/2;                         % Nyquist frequency
        hz = linspace(0,nyquist,floor(N/2)+1);  % Frequencies in Hz
        %% Plot power spectrum of each signal
        figure; 
        pow1 = abs(fft(signal1)/N).^2;
        pow2 = abs(fft(signal2)/N).^2;
        subplot(211); plot(hz,pow1(1:length(hz)),'linew',1.1); set(gca,'xlim',[0 20]);
        xlabel('Frequency (Hz)'); ylabel('Power');
        title('Frequency vs Power'); legend({'signal1'});
        subplot(212); plot(hz,pow2(1:length(hz)),'linew',1.1); set(gca,'xlim',[0 20]);
        xlabel('Frequency (Hz)'); ylabel('Power');
        title('Frequency vs Power');legend({'signal2'});
        %% Manual Welch's method
        % parameters   
        winlen = floor(fs/2);                   % window length in seconds*fs
        ovelp = 0;%floor(winlen/2);                % number of points of overlap    
        winon = 1:winlen-ovelp:N-winlen;        % window onset times
        hzW = linspace(0,nyquist,floor(winlen/2)+1);    % diff-length signal needs a diff Hz vector
        hannw = .5-cos(2*pi*linspace(0,1,winlen))./2;   % Hann taper window
        welchspec = zeros(length(hzW),1);	% initialize the power matrix (windows x frequencies)
        figure;
        plot(real(signal1));
        [b2, a2] = butter(order, (2*hp)/fs, 'high');            % filter coefficients
        for wi=1:length(winon)                                  % loop over frequencies
            datachunk = signal1(winon(wi):winon(wi)+winlen-1);  % chunk of data from time window
            datachunk = filtfilt(b2, a2, datachunk);            % high-pass filter data
            datachunk = datachunk .* hannw';                    % apply Hann taper to data
            hold on;
            plot(winon(wi):winon(wi)+winlen-1,real(datachunk));
            plot(winon(wi):winon(wi)+winlen-1,imag(datachunk));
            pause(.1)     
            tmppow = abs(fft(datachunk)/winlen).^2;             % compute its power
            welchspec = welchspec  + tmppow(1:length(hzW));     % enter into matrix
        end
        welchspec = welchspec/length(winon);                    % divide by number of windows
        fftspec = abs( fft(signal1)/N ).^2;                 	% generate FFT of entire signal
        %% Plotting Methods
        figure; subplot(211); hold on;
        plot(hz,fftspec(1:length(hz)),'--b','linew',1.1);
        plot(hzW,welchspec,'k','linew',1.15);
        set(gca,'xlim',[0 30]);
        xlabel('Frequency (Hz)');
        legend({'"Static FFT';'Welch''s method'});
        title('Using FFT and Welch''s');
        %% Pwelch Matlab function
        subplot(212);
        winsize = floor(fs/2.5);                            % 400 samples
        wind = .5 - cos(2*pi*linspace(0,1,winsize))./2;     % create Hann window
        spectres = 0.5;
        nfft = round(fs/spectres);                          % number of DFT points
        noverlap = floor(winsize/2);                    	% 75% overap
        pwelch(signal1,wind,noverlap,nfft,fs);                % Apply Welch method
         %}

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Coherence Analysis
        %{
        %% Signal declaration
        signal = reshape(gestF{1,1}(:,:,:),[],4,1);
        signal1 = signal(:,1);
        signal2 = signal(:,3);  
        %% Parameters
        N = length(real(signal1));             	% length of sequence
        nyquist	= fs/2;                         % Nyquist frequency
        hz = linspace(0,nyquist,floor(N/2)+1);  % Frequencies in Hz
        %% Plot power spectrum of each signal
        figure(1); 
        pow1 = abs(fft(signal1)/N).^2;
        pow2 = abs(fft(signal2)/N).^2;
        subplot(211); plot(hz,pow1(1:length(hz)),'linew',1.1); set(gca,'xlim',[0 20]);
        xlabel('Frequency (Hz)'); ylabel('Power');
        title('Frequency vs Power'); legend({'signal1'});
        subplot(212); plot(hz,pow2(1:length(hz)),'linew',1.1); set(gca,'xlim',[0 20]);
        xlabel('Frequency (Hz)'); ylabel('Power');
        title('Frequency vs Power');legend({'signal2'});
        %% Manual Welch's method
        % Parameters   
        winlen = floor(fs/(1/5));                       % window length in seconds*fs
        ovelp = floor(winlen*0.5);                      % number of points of overlap    
        winon = 1:winlen-ovelp:N-winlen;                % window onset times
        hzW = linspace(0,nyquist,floor(winlen/2)+1);    % diff-length signal needs a diff Hz vector
        hannw = .5-cos(2*pi*linspace(0,1,winlen))./2;   % Hann taper window
        [b2, a2] = butter(order, (2*hp)/fs, 'high');	% filter coefficients
        welch_cross = zeros(length(hzW),1);
        welch_powX = zeros(length(hzW),1);
        welch_powY = zeros(length(hzW),1);
    %     figure(2);
        for wi=1:length(winon)                          % loop over frequencies
            datachunkX = signal1(winon(wi):winon(wi)+winlen-1); % data-chunk - time window - Signal 1
            datachunkY = signal2(winon(wi):winon(wi)+winlen-1); % data-chunk - time window - Signal 1
            %% Filter DC Offset
            datachunkX = filtfilt(b2, a2, datachunkX);      	% high-pass filter data
            datachunkX = datachunkX .* hannw';                  % apply Hann taper to data
            datachunkY = filtfilt(b2, a2, datachunkY);          % high-pass filter data
            datachunkY = datachunkY .* hannw';                  % apply Hann taper to data
            %% Compute power Spectrums
            fft_x = fft(datachunkX);
            fft_y = fft(datachunkY);
            cross_spc = (fft_x.*conj(fft_y))./winlen^2;         % Cross-spectrum
            pow_x = (fft_x.*conj(fft_x))./winlen^2;             % Power Spectrum X
            pow_y = (fft_y.*conj(fft_y))./winlen^2;             % Power Spectrum Y
    %         pow_x = abs(fft_x./winlen).^2;                    % Power Spectrum X
    %         pow_y = abs(fft_y./winlen).^2;                    % Power Spectrum Y
            welch_cross = welch_cross + cross_spc(1:length(hzW));
            welch_powX = welch_powX + pow_x(1:length(hzW));
            welch_powY = welch_powY + pow_y(1:length(hzW)); 
            %% Plot Coherence Steps
            %{
            figure;hold on;
            plot(welch_powX);plot(welch_powY);
            plot(welch_powX.*welch_powY);
            plot(abs(welch_cross).^2);
            set(gca,'xlim',[0 20]);
            legend({'X','Y','Spec-Prod','Cross'});
            figure;
            plot((abs(welch_cross).^2)./(welch_powX.*welch_powY));
            title('MSCoherence');
            %}
            %% Plot Cross-Spectrum Addition
            %{
            clf;
            plot(abs(welch_cross),'linew',1.1);
            title('Crss - Spectrum X-Y');
            set(gca,'xlim',[0 20]);
            pause(0.1);
            %}
        end
        %% Average Spectrums
        welch_cross = abs(welch_cross/length(winon)).^2; 
        welch_powX = welch_powX/length(winon); 
        welch_powY = welch_powY/length(winon);
        figure(3); hold on;
        plot(welch_powX(1:length(hzW)));plot(welch_powY(1:length(hzW)));
        plot(welch_cross(1:length(hzW)));
        plot((welch_powX(1:length(hzW)).*welch_powY(1:length(hzW))))
        set(gca,'xlim',[0 20]);
        legend({'Pow-X','Pow-Y','Cross-Spectrum X-Y','Spec-Prod'});
        title('Power&Cross-Spectrums X-Y');
        %% Calculate and Plot Magnitude Square Coherence
        power_specs = (welch_powX(1:length(hzW)).*welch_powY(1:length(hzW)));
        power_specs( power_specs < 1 ) = 1;                     % Set values less than 1 to one to avoid ripple
        msc = welch_cross(1:length(hzW))./power_specs;
        figure(4); 
        subplot(211); plot(hzW,msc,'linew',1.1); set(gca,'xlim',[0 30]);
        title('Magnitude Square Coherence - Manual Method');
        [Cxy,nf] = mscohere(signal1,signal2,winlen,ovelp,winlen,fs);
        subplot(212); plot(nf,Cxy,'linew',1.1); 
        title('ms-coherence - Matlab');
         %}

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %% Save Trial Space to File
%         save(strcat('Par',partID,'/',partID,'gestF1-',int2str(numberTrials),'.mat'),'gestF');
    %}
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Coherence Calculation - mean/median of Trials
        %{
        coh = cell(fl,numGestures);
        winsize = floor(fs/(1/2));                      	% window length in seconds*fs
        ovelp = floor(winsize*0.9);                     	% number of points of overlap 
        % Loop over Force Levels and Gestures
        for m=1:fl
            for w=1:numGestures
                coh{m,w} = zeros(((winsize/2)+1),6,numberTrials);
                for in=1:numberTrials
                    conc_trials = gestF{m,w}(:,:,in);
                    conc_trials = [complex(zeros(winsize/2,4),zeros(winsize/2,4)) ;...
                    conc_trials ; complex(zeros(winsize/2,4),zeros(winsize/2,4))];
                    [Cxy,nf] = manual_coherence(conc_trials(:,1),conc_trials(:,2),winsize,ovelp,fs);
                    coh{m,w}(:,1,in) = Cxy;
                    clear Cxy; clear nf;
                    [Cxy,nf] = manual_coherence(conc_trials(:,1),conc_trials(:,3),winsize,ovelp,fs);
                    coh{m,w}(:,2,in) = Cxy;
                    clear Cxy; clear nf; 
                    [Cxy,nf] = manual_coherence(conc_trials(:,1),conc_trials(:,4),winsize,ovelp,fs);
                    coh{m,w}(:,3,in) = Cxy;
                    clear Cxy; clear nf;
                    [Cxy,nf] = manual_coherence(conc_trials(:,2),conc_trials(:,3),winsize,ovelp,fs);
                    coh{m,w}(:,4,in) = Cxy;
                    clear Cxy; clear nf;
                    [Cxy,nf] = manual_coherence(conc_trials(:,2),conc_trials(:,4),winsize,ovelp,fs);
                    coh{m,w}(:,5,in) = Cxy;
                    clear Cxy; clear nf;
                    [Cxy,nf] = manual_coherence(conc_trials(:,3),conc_trials(:,4),winsize,ovelp,fs);
                    coh{m,w}(:,6,in) = Cxy; 
                    clear Cxy;
                end
                coh{m,w} = mean(coh{m,w},3);
            end
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Coherence Calculation - Concatenation of Trials
%         %{
        coh = cell(fl,numGestures);
        winsize = floor(fs/(2/5));                      	% window length in seconds*fs
        ovelp = floor(winsize*0.9);                     	% number of points of overlap
        % Loop over Force Levels and Gestures
        for m=1:fl
            for w=1:numGestures
%                 conc_trials = reshape(gestF{m,w},[],4);
                conc_trials = num2cell(gestF{m,w}, [1 2]); %split A keeping dimension 1 and 2 intact
                conc_trials = vertcat(conc_trials{:});
                conc_trials = [complex(zeros(winsize/2,4),zeros(winsize/2,4)) ;...
                    conc_trials ; complex(zeros(winsize/2,4),zeros(winsize/2,4))];
                [Cxy,nf] = manual_coherence(conc_trials(:,1),conc_trials(:,2),winsize,ovelp,fs);
                coh{m,w}(:,1) = Cxy;
                clear Cxy; clear nf;
                [Cxy,nf] = manual_coherence(conc_trials(:,1),conc_trials(:,3),winsize,ovelp,fs);
                coh{m,w}(:,2) = Cxy;
                clear Cxy; clear nf;
                [Cxy,nf] = manual_coherence(conc_trials(:,1),conc_trials(:,4),winsize,ovelp,fs);
                coh{m,w}(:,3) = Cxy;
                clear Cxy; clear nf;
                [Cxy,nf] = manual_coherence(conc_trials(:,2),conc_trials(:,3),winsize,ovelp,fs);
                coh{m,w}(:,4) = Cxy;
                clear Cxy; clear nf;
                [Cxy,nf] = manual_coherence(conc_trials(:,2),conc_trials(:,4),winsize,ovelp,fs);
                coh{m,w}(:,5) = Cxy;
                clear Cxy; clear nf;
                [Cxy,nf] = manual_coherence(conc_trials(:,3),conc_trials(:,4),winsize,ovelp,fs);
                coh{m,w}(:,6) = Cxy;
                clear Cxy;
            end
        end
        %}
    %{
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plot Coherence Results
        %{
        figure;
        for j=1:numGestures
            subplot(2,numGestures/2,j);
            for i=1:6   % Muscle pairs
                plot(nf,coh{1,j}(:,i),'LineWidth',1.1);hold on;
            end
            hold off;
            legend({'FCR-BRA','FCR-EDC','FCR-FCU','BRA-EDC','BRA-FCU','EDC-FCU'});
            title(strcat(" MS-Coherence - Gesture: ",GestList{j}));
            ylabel('Magnitude-Squared Coherence');
            xlabel('Frequency [Hz]');
        end
        %}

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save Cohere Matrix
%         save(strcat('Par',partID,'/',partID,'coherence.mat'),'coh');
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Average of Coherence
        %{  
        copges = cell(fl,numGestures);
        suma = [0];
        musclePairs = 6;
        for m=1:fl
            for w=1:numGestures
                for i=1:musclePairs
                    for t=1:numberTrials
                        suma = suma + coh{m,w}{1,t}(:,i);
                    end
                    copges{m,w}(:,i) = suma/numberTrials;
                    suma = [0];
                end
            end
        end
        %}
    %}
        coh_part{1,part} = coh;
    end
    %% Save Feature Space to File
    save(strcat(titl,'-coh.mat'),'coh_part');
    save(strcat(titl,'-nf.mat'),'nf');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coherence - Gesture
coh_ges = cell(fl,numGestures);
for m=1:fl
    for w=1:numGestures
        coh_temp = zeros(length(nf),6,fin-from+1);
        ParList = [];
        for i=from:fin
            coh_temp(:,:,i-from+1) = coh_part{1,i}{m,w};
            ParList = [ParList "Participant "+int2str(i)];
        end
        coh_ges{m,w} = mean(coh_temp,3);
%         plt_coh_temp(GestList(w),cellstr(ParList),nf,reshape(coh_temp(:,2,:),[],fin-from+1));    % Plot coh_temp
    end
end
%% Coherence - Group
coh_grp = cell(fl,1);
for m=1:fl                              % Force Levels
    for i=1:6                           % Muscle pairs
        coh_temp = zeros(length(nf),6,numGestures);
        for j=1:numGestures             % Number of gestures
            coh_temp(:,:,j) = coh_ges{m,j};
        end
        coh_grp{m} = mean(coh_temp,3);
%         plt_coh_temp(strcat("MusclePair: "+int2str(i)),cellstr(GestList),nf,reshape(coh_temp(:,1,:),[],6));    % Plot coh_temp
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Coherence - Gesture Results
%{
figure;%('DefaultAxesFontSize',17);
for m=1:fl
    for j=1:numGestures
        subplot(fl,numGestures,j+numGestures*(m-1));
        plot(nf,coh_ges{m,j},'LineWidth',1.1);
        set(gca,'xlim',[0 40]);
%         legend({'FCR-BRA','FCR-EDC','FCR-FCU','BRA-EDC','BRA-FCU','EDC-FCU'});
        title({titl,strcat(" Force:",int2str(m)," - ",GestList(j))});
        ylabel('MSC');
        xlabel('Frequency [Hz]');
    end
end
%}
%% Plot Coherence of Gestures at each Force Level
%{
for w=1:numGestures
    figure;
    for m=1:fl
        subplot(fl,1,m);
        plot(nf,coh_ges{m,w},'LineWidth',1.1);
        set(gca,'xlim',[0 50]);
        legend({'FCR-BRA','FCR-EDC','FCR-FCU','BRA-EDC','BRA-FCU','EDC-FCU'});
        title({strcat("MSCoherence - Average of ",titl),strcat("Gesture: ",GestList(w)," - ForceLevel: ",int2str(m))});
        ylabel('MSC');
        xlabel('Frequency [Hz]');
    end
end
%}
%% Plot Coherence - per Participant
%{
for par=from:fin
    figure;
    for j=1:numGestures
        subplot(1,numGestures,j);
        plot(nf,coh_part{1,par}{1,j},'LineWidth',1.1);hold on;
        set(gca,'xlim',[0 50]);
        legend({'FCR-BRA','FCR-EDC','FCR-FCU','BRA-EDC','BRA-FCU','EDC-FCU'});
        if par>=11, lab = "Amputee: ";
        else, lab = "Participant: "; end
        title(strcat(lab,int2str(par-from+1)," - MSC - ",GestList{j}));
        ylabel('Magnitude-Squared Coherence');
        xlabel('Frequency [Hz]');
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Append Coherence - All gestures - All participants
%{
conc_coh_grp = zeros(length(coh_ges{1,1}),numGestures*6,fl);
for m=1:fl
    for w=1:numGestures
        conc_coh_grp(:,1+6*(w-1):6*w,m) = coh_ges{m,w};
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-negative Matrix Factorization - Group Results
% %{
%%%%%% NNMF Group %%%%%%
k = 3;                          % Number of components
freqs = [4,15,150];
nnmf_grp = cell(fl,1);
nnmf_ges = cell(fl,numGestures);
nnmf_part = cell(1,fin-from+1);           
for m=1:fl
    [nnmf_grp{m,1}.H,nnmf_grp{m,1}.W] = nnmf_mat(m,1,0,k,coh_grp,freqs,nf,fl,numGestures,titl,GestList,0);
end
%%%%%% NNMF gesture %%%%%%
for m=1:fl                              
    for w=1:numGestures
        [nnmf_ges{m,w}.H, nnmf_ges{m,w}.W] = nnmf_mat(m,w,0,k,coh_ges,freqs,nf,fl,numGestures,titl,GestList,0);
    end
end
%%%%% NNMF Participant %%%%%%
for i=from:fin
%     figure;
    for m=1:fl
        for w=1:numGestures
            [nnmf_part{1,i}(m,w).H ,nnmf_part{1,i}(m,w).W] = nnmf_mat(m,w,i,k,coh_part,freqs,nf,fl,numGestures,titl,GestList,0);
        end
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Connectivity Matrix and Metrics - Group Results
% %{
try
    ad = addpath("C:\Users\sebas\Downloads\BCT\2019_03_03_BCT");
    delete(ad);
catch
    addpath("C:\Users\csm116\Downloads\BCT\2019_03_03_BCT");
end
clc;
s = {'FCR' 'FCR' 'FCR' 'BRA' 'BRA' 'EDC'};
t = {'BRA' 'EDC' 'FCU' 'EDC' 'FCU' 'FCU'};
varnames = {'G', 'A', 'BC', 'CC', 'GE', 'ME', 'NS', 'ED','SD'};
varDescr = {'Connection Matrix', 'Adjancency Matrix', 'Betweenness Centrality',...
    'Clustering Coefficient', 'Global Efficiency', 'Median of Edges', 'Node Strength',...
    'Edges','Standard Deviation'};
th_op = 0;      % 0:no_thresholding / 1:std-method / 2:top%-method / 3:absoluteVal-method
th_val = 0.1;   % for 1:Number of stds / 2:% of top edges / 3: min thres value
adjmat_grp = cell(fl,1);
for m=1:fl
    adjmat_grp{m,1} = connect_mat(m,1,0,k,s,t,nnmf_grp,th_op,th_val,varnames);
end
%% Connectivity Matrix and Metrics - per Gesture
adjmat_ges = cell(fl,numGestures);
for m=1:fl
    for w=1:numGestures
        adjmat_ges{m,w} = connect_mat(m,w,0,k,s,t,nnmf_ges,th_op,th_val,varnames);
    end
end
%% Connectivity Matrix and Metrics - per Participant
adjmat_part = cell(1,fin-from+1);
for i=from:fin
    for m=1:fl
        for w=1:numGestures
            adjmat_part{1,i}{m,w} = connect_mat(m,w,i,k,s,t,nnmf_part,th_op,th_val,varnames);
        end
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boxplots
numCommGes = 5;
%% Differential Boxplots - Upload amp adjmat_ges
%{
try
    adjmatAmp_ges = load(strcat('Amputees','-adjmat-ges.mat'));
    adjmatAmp_ges = adjmatAmp_ges.adjmat_ges;
catch
    if (who=="amp")
        save(strcat(titl,'-adjmat-ges.mat'),'adjmat_ges');
    end
end
% Generate Differential Boxplot
[bx_CC,bx_NS,bx_ED] = bx_plt_prep2(numCommGes,k,adjmat_ges,adjmatAmp_ges);
bx_metrics = {bx_CC,bx_NS,bx_ED};
type = 'Differential Plot';
bx_plt(bx_metrics,k,varDescr,GestList,titl,1,'t');
%}
%% Boxplots Average of Coherence across Participants
% %{
[bx_CC,bx_NS,bx_ED] = bx_plt_prep2(numCommGes,k,adjmat_ges);
bx_metrics = {bx_CC,bx_NS,bx_ED};
bx_plt(bx_metrics,k,varDescr,GestList,titl,0,'t');
%}
%% Boxplots Average of Network across Participants 
% %{
[bx_CC,bx_NS,bx_ED] = bx_plt_prep(numCommGes,from,fin,adjmat_part);
bx_CC = mean(bx_CC(:,:,:,from:fin),4);bx_CC(bx_CC==0) = NaN;        
bx_NS = mean(bx_NS(:,:,:,from:fin),4);bx_NS(bx_NS==0) = NaN;
bx_ED = mean(bx_ED(:,:,:,from:fin),4);bx_ED(bx_ED==0) = NaN;
bx_metrics = {bx_CC,bx_NS,bx_ED};
bx_plt(bx_metrics,k,varDescr,GestList,titl,0,'t');
%}
%% Boxplots All Participants using Participants
%{
% GE
[bx_CC,bx_NS,bx_ED,bx_GE] = bx_plt_prep(numCommGes,from,fin,adjmat_part);
tmp_bx = num2cell(bx_GE(:,:,:,from:fin), [1 2 3]);
bx_GE = vertcat(tmp_bx{:});bx_GE(bx_GE==0) = NaN;
% CC
tmp_bx = num2cell(bx_CC(:,:,:,from:fin), [1 2 3]);
bx_CC = vertcat(tmp_bx{:});bx_CC(bx_CC==0) = NaN;
% NS
tmp_bx = num2cell(bx_NS(:,:,:,from:fin), [1 2 3]);
bx_NS = vertcat(tmp_bx{:});bx_NS(bx_NS==0) = NaN;
% ED
tmp_bx = num2cell(bx_ED(:,:,:,from:fin), [1 2 3]);
bx_ED = vertcat(tmp_bx{:});bx_ED(bx_ED==0) = NaN;
% Plot
bx_metrics = {bx_CC,bx_NS,bx_ED,bx_GE};
bx_plt(bx_metrics,k,varnames,GestList,titl,0);
%}
%% Boxplots Average of Nodes using Participants
%{
[bx_CC,bx_NS,bx_ED] = bx_plt_prep(numCommGes,from,fin,adjmat_part);
% CC
bx_CC = mean(bx_CC(:,:,:,from:fin),1); 
bx_CC(bx_CC==0) = NaN; bx_CC = squeeze(bx_CC);
bx_CC = permute(bx_CC, [3 1 2]);
% NS
bx_NS = mean(bx_NS(:,:,:,from:fin),1);
bx_NS(bx_NS==0) = NaN; bx_NS = squeeze(bx_NS);
bx_NS = permute(bx_NS, [3 1 2]);
% ED
bx_ED = mean(bx_ED(:,:,:,from:fin),1);
bx_ED(bx_ED==0) = NaN; bx_ED = squeeze(bx_ED);
bx_ED = permute(bx_ED, [3 1 2]);
% Plot Boxplot
bx_metrics = {bx_CC,bx_NS,bx_ED};
bx_plt(bx_metrics,k,varnames,GestList,titl,0);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group Connectivity Matrix 
% plt_connmat_grp(titl,adjmat_grp);
%% per Gesture
% %{
for gest=1:numGestures
    plt_connmat_grp(GestList{gest},adjmat_ges(:,gest)); 
end
%}
%% per Participant
%{
par = 11;           % Participant
gest = 3;          % Gesture
plt_connmat_par(GestList{gest},par,adjmat_part{1,par}(:,gest));   
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Heatmap Matrix of Average of Coherence
for w=1:numCommGes
    figure;
    for i=1:k
        subplot(1,k,i);
        h = heatmap(adjmat_ges{1,w}.A{i,1},'Colormap', flipud(hot));
        caxis([0 round(max(cell2mat(adjmat_ges{1,w}.A),[],'all'),2)]);
        h.Title = {strcat(titl,' - ',GestList(w)),strcat('Adjancency Matrix: k',int2str(i))};
        h.YData = {'FCR','BRA','EDC','FCU'};
        h.XData = h.YDisplayData;
    end
end
%% Heatmap Matrix of Average of Network
tmp = 0;
A_av = cell(k,numGestures);
for w=1:numGestures
    for i=1:k
        for m=from:fin
            tmp = tmp + adjmat_part{1,m}{1,w}.A{i,1};
        end
        tmp = tmp/(fin-from+1);
        A_av{i,w} = tmp;
    end
end
% Generate Heatmap Matrix
for w=1:numCommGes
    figure;
    for i=1:k
        subplot(1,k,i);
        h = heatmap(A_av{i,w},'Colormap', flipud(hot)); %flipud(jet(25)
        caxis([0 round(max(cell2mat(A_av(:,w)),[],'all'),2)]);
        h.Title = {strcat(titl,' - ',GestList(w)),strcat('Adjancency Matrix: k',int2str(i))};
        h.YData = {'FCR','BRA','EDC','FCU'};
        h.XData = h.YDisplayData;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BarPlot Network Metrics - Group Results 
%{
figure;
X = categorical(varnames(3:end));
X = reordercats(X,varnames(3:end));
for m=1:fl
    subplot(fl,1,m);
    tbl = adjmat_grp{m,1};
    Y = [cellfun(@mean, tbl.BC)';cellfun(@mean, tbl.CC)';tbl.GE';tbl.ME';....
        cellfun(@mean, tbl.NS)';cellfun(@mean, tbl.ED)';tbl.SD'];
    b = bar(X,Y);
    legnd = cell.empty(0,k);
    for i=1:k
        legnd = [legnd,{strcat('Component-',int2str(i))}];
        xtips1 = b(i).XEndPoints;
        ytips1 = b(i).YEndPoints;
        labels1 = string(b(i).YData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
    legend(legnd,'Location','Northwest');
    title({"ForceLevel: " + int2str(m),"Network Metrics - " + titl});
end
%}

%% BarPlot Network Metrics - per Gesture
%{
legnd = cell.empty(0,k);
for i=1:k
    legnd = [legnd,{strcat('Component-',int2str(i))}];
end
a = cell(1,k);
X = categorical(legnd);
X = reordercats(X,legnd);
figure;
for m=1:fl
    for i=1:k
        tbl_ges = cell2table(cell(0,length(varnames(3:end))), 'VariableNames', varnames(3:end));
        for j=1:numGestures
            tbl = adjmat_ges{m,j};
            tbl_ges = [tbl_ges;{cellfun(@mean, tbl.BC(i)),cellfun(@mean, tbl.CC(i)),...
                tbl.GE(i),tbl.ME(i),cellfun(@mean, tbl.NS(i)),cellfun(@mean, tbl.ED(i)),tbl.SD(i)}];
        end
        a{1,i} = tbl_ges;
    end
    BC = double.empty(0,numGestures);
    CC = double.empty(0,numGestures);
    GE = double.empty(0,numGestures);
    ME = double.empty(0,numGestures);
    NS = double.empty(0,numGestures);
    ED = double.empty(0,numGestures);
    SD = double.empty(0,numGestures);
    for i=1:k
        BC = [BC;a{1,i}.BC'];
        CC = [CC;a{1,i}.CC'];
        GE = [GE;a{1,i}.GE'];
        ME = [ME;a{1,i}.ME'];
        NS = [NS;a{1,i}.NS'];
        ED = [ED;a{1,i}.ED'];
        SD = [SD;a{1,i}.SD'];
    end
    % Generate Plot
    met = length(varnames(3:end));
    for i=1:met
        subplot(fl,met,i+met*(m-1));
        if(i==1);bar(X,BC(:,1:4)); ylim([0 2.5]); tl='BC';end
        if(i==2);bar(X,CC(:,1:4)); ylim([0 0.5]); tl='CC';end
        if(i==3);bar(X,GE(:,1:4)); ylim([0 0.5]); tl='GE';end
        if(i==4);bar(X,ME(:,1:4)); ylim([0 0.5]); tl='ME';end
        if(i==5);bar(X,NS(:,1:4)); ylim([0 1.5]); tl='NS';end
        if(i==6);bar(X,ED(:,1:4)); ylim([0 0.5]); tl='ED';end
        if(i==7);bar(X,SD(:,1:4)); ylim([0 0.4]); tl='SD';end
%         legend(GestList)
        title({titl + " -  ForceLevel: " + int2str(m),tl});
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


