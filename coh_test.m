close all;
clear all;
clc;
%% Parameters
fs = 1000;                      % Smapling rate
t = 0:1/fs:15;          
a1 = 50;                        % Amplitude of signals
winlen = 1024;                  % window length in seconds*fs
noverlap = floor(winlen/2);     % number of points of overlap
%% Plot sinewaves with noise
f1=4;f2=12;f3=25;f4=9;          % Frequencies
% Sine waves
y1 = a1*sin(2*pi*t*f1); 
y2 = a1*sin(2*pi*t*f2);
y3 = a1*sin(2*pi*t*f3);
y4 = a1*sin(2*pi*t*f4);
% Signals
[b2, a2] = butter(4, 30/(fs/2), 'low');
sig1 = y1+y2+y3;sig1 = awgn(sig1,10,'measured'); sig1 = sig1';
sig2 = y1+y4+y3;sig2 = awgn(sig2,20,'measured'); sig2 = sig2';
filt_test_sig1 = filtfilt(b2, a2, sig1);
filt_test_sig2 = filtfilt(b2, a2, sig2);
% Plot signals
figure; tiledlayout(3,2,'Padding','compact','TileSpacing','tight');
nexttile();plot(sig1,'color','red');   hold on;   plot(sig2,'color','blue'); set(gca, 'xlim', [0 500]);xlabel('samples');
title(strcat("Signals with common freqs. at ",strcat(int2str(f1)," and ",int2str(f3)," Hz")));
nexttile(); plot(filt_test_sig1,'color','red');   hold on;   plot(filt_test_sig2,'color','blue'); set(gca, 'xlim', [0 500]);
title("Signals low-pass filtered at 30 Hz");xlabel('samples');
% Plot Welch's Power Spectral Density      
nexttile();pwelch(sig1,winlen,noverlap,[],fs);   h = get(gca, 'Children');set(h(1), 'Color', 'red','linew',1.1);hold on;
pwelch(sig2,winlen,noverlap,[],fs); h = get(gca, 'Children');set(h(1), 'Color', 'blue', 'lineStyle','--','linew',1.1); set(gca, 'xlim', [0 100])
title('Welch Power Spectral Density Estimate');
legend({strcat("signal_1 - freqs: ",int2str(f1),",",int2str(f2),",",int2str(f3)," Hz"),...
    strcat("signal_2 - freqs: ",int2str(f1),",  ",int2str(f4),",",int2str(f3)," Hz")},'FontSize',10,'location','northEast');
nexttile();pwelch(filt_test_sig1,winlen,noverlap,[],fs);   h = get(gca, 'Children');set(h(1), 'Color', 'red','linew',1.1);hold on;
pwelch(filt_test_sig2,winlen,noverlap,[],fs); h = get(gca, 'Children');set(h(1), 'Color', 'blue', 'lineStyle','--','linew',1.1);set(gca, 'xlim', [0 100])
title('Welch Power Spectral Density Estimate');
%% Plot Coherence
% Plot Magnitude Square Coherence built in function
nexttile()
[Cxy,nf] = mscohere(sig1,sig2,nuttallwin(winlen),noverlap,[],fs);
plot(nf,Cxy,'linew',1.25,'Color', '#EDB120'); hold on;
[mmsc, freq] = manual_coherence(sig1,sig2,winlen,noverlap,fs);
plot(freq,mmsc,'linew',1.75,'lineStyle','-.','Color', '#7E2F8E'); set(gca, 'xlim', [0 200]);
title('Coherene Estimate');ylabel('Magnitude-Square Cohrence');xlabel('Frequency (Hz)');
nexttile()
[Cxy,nf] = mscohere(filt_test_sig1,filt_test_sig2,nuttallwin(winlen),noverlap,[],fs);
plot(nf,Cxy,'linew',1.25,'Color', '#EDB120'); hold on;
[mmsc, freq] = manual_coherence(filt_test_sig1,filt_test_sig2,winlen,noverlap,fs);
plot(freq,mmsc,'linew',1.75,'lineStyle','-.','Color', '#7E2F8E'); set(gca, 'xlim', [0 200]);
title('Coherene Estimate');ylabel('Magnitude-Square Coherence');xlabel('Frequency (Hz)');
legend({'mscohere - Matlab function','Modified-MSC'},'FontSize',11);


%% White noise
L = 15000;
mu = 0;
sigma = 50;
figure; tiledlayout(3,2,'Padding','compact','TileSpacing','tight');
test_sig1 = sigma*randn(L,1)+mu;
test_sig2 = sigma*randn(L,1)+mu;
[b2, a2] = butter(4, 30/(fs/2), 'low');             % filter coefficients
filt_test_sig1 = filtfilt(b2, a2, test_sig1);
filt_test_sig2 = filtfilt(b2, a2, test_sig2);
% Plot Signals
nexttile();plot(test_sig1,'color','#0072BD');hold on; plot(test_sig2,'color','#D95319');set(gca, 'xlim', [0 15000]);
title('White noise data');xlabel('samples');
nexttile();plot(filt_test_sig1,'color','#0072BD');hold on; plot(filt_test_sig2,'color','#D95319');set(gca, 'xlim', [0 15000]);
title('White noise data low-pass filtered at 30 Hz');xlabel('samples');
% Plot Welch's Power Spectral Density
nexttile();pwelch(test_sig1,winlen,noverlap,[],fs);hold on; pwelch(test_sig2,winlen,noverlap,[],fs);
h = get(gca, 'Children');set(h(1), 'Color', '#0072BD','linew',1.1);set(h(2), 'Color', '#D95319','linew',1.1);
nexttile();pwelch(filt_test_sig1,winlen,noverlap,[],fs);hold on; pwelch(filt_test_sig2,winlen,noverlap,[],fs);
h = get(gca, 'Children');set(h(1), 'Color', '#0072BD','linew',1.1);set(h(2), 'Color', '#D95319','linew',1.1);
% Plot MSC of signals
nexttile();[Cxy,nf]=mscohere(test_sig1,test_sig2,nuttallwin(winlen),noverlap,[],fs);plot(nf,Cxy,'linew',1.1,'Color', '#EDB120');
hold on;[mmsc, freq] = manual_coherence(test_sig1,test_sig2,winlen,noverlap,fs);
plot(freq,mmsc,'linew',1.1, 'Color', '#7E2F8E');
title('Coherence Estimate via Welch');ylabel('Magnitude-Square Cohrence');xlabel('Frequency');
nexttile();[Cxy,nf]=mscohere(filt_test_sig1,filt_test_sig2,nuttallwin(winlen),noverlap,[],fs);plot(nf,Cxy,'Color', '#EDB120');
hold on;[mmsc, freq] = manual_coherence(filt_test_sig1,filt_test_sig2,winlen,noverlap,fs); 
plot(freq,mmsc,'linew',1.1, 'Color', '#7E2F8E');
title('Coherence Estimate via Welch');ylabel('Magnitude-Square Cohrence');xlabel('Frequency');
legend({'mscohere - Matlab function','Modified-MSC'},'FontSize',11);


