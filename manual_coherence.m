function [msc, freq] = manual_coherence(signal1,signal2,winlen,ovelp,fs)
    % Manual_coherence estimates the Magnitude Square Coherence of two functions using Welch's method.
    N = length(signal1);                            % length of sequence
    nyquist	= fs/2;                                 % Nyquist frequency
    %% Manual Welch's method
    % Parameters
    winon = 0:winlen-ovelp:N-winlen;                % window onset times
    hzW = linspace(0,nyquist,floor(winlen/2)+1);    % diff-length signal needs a diff Hz vector
%     hannw = .5-cos(2*pi*linspace(0,1,winlen))./2;   % Hann taper window
%     hannw = blackmanharris(winlen)';
    hannw = nuttallwin(winlen)';
%     hannw = tukeywin(winlen,0.05)';
%     hannw = hamming(winlen)';
    [b2, a2] = butter(4, 1/(fs/2), 'high');             % filter coefficients
    welch_cross = zeros(length(hzW),1);
    welch_powX = zeros(length(hzW),1);
    welch_powY = zeros(length(hzW),1);
    % Loop over frequencies
    for wi=1:length(winon)
        datachunkX = signal1(winon(wi)+1:winon(wi)+winlen); % data-chunk - time window - Signal 1
        datachunkY = signal2(winon(wi)+1:winon(wi)+winlen); % data-chunk - time window - Signal 2
        %% Filter DC Offset
        datachunkX = filtfilt(b2, a2, datachunkX);      	% high-pass filter data 1
        datachunkX = datachunkX .* hannw';                  % apply Hann taper to data 1
        datachunkY = filtfilt(b2, a2, datachunkY);          % high-pass filter data 2
        datachunkY = datachunkY .* hannw';                  % apply Hann taper to data 2
        %% Compute Power Spectrums
        fft_x = fft(datachunkX);
        fft_y = fft(datachunkY);
        cross_spc = (fft_x.*conj(fft_y))./winlen^2;         % Cross-spectrum
        pow_x = (fft_x.*conj(fft_x))./winlen^2;             % Power Spectrum X
        pow_y = (fft_y.*conj(fft_y))./winlen^2;             % Power Spectrum Y
        % Add the spectrums
        welch_cross = welch_cross + cross_spc(1:length(hzW));
        welch_powX = welch_powX + pow_x(1:length(hzW));
        welch_powY = welch_powY + pow_y(1:length(hzW));
    end
    %{
    welch_cross = welch_cross/length(winon);       	% Average cross-spectrum
    welch_powX = welch_powX/length(winon);      	% Average power-spectrum X
    welch_powY = welch_powY/length(winon);          % Average power-spectrum Y
    % Compute MSC
    power_specs = sqrt(welch_powX(1:length(hzW))).*sqrt(welch_powY(1:length(hzW)));
    power_specs( power_specs < 0.1 ) = 1;             % Set small values <1 to 1 to avoid ripple
    msc = welch_cross(1:length(hzW))./power_specs;
    freq = hzW;
end
%}
% %{
    welch_cross = abs(welch_cross/length(winon)).^2;       	% Average cross-spectrum
    welch_powX = welch_powX/length(winon);                  % Average power-spectrum X
    welch_powY = welch_powY/length(winon);                  % Average power-spectrum Y
    % Compute MSC
    power_specs = (welch_powX(1:length(hzW)).*welch_powY(1:length(hzW)));
    power_specs( power_specs < 1 ) = 1;                     % Set small values <1 to 1 to avoid ripple
    msc = welch_cross(1:length(hzW))./power_specs;
    freq = hzW;
end
%}