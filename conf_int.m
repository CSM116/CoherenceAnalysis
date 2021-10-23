ci_method = 1;
ovlp = 0.825;
alpha = 0.05;
nfft = 512;
L = ceil(15000/nfft);

if (ci_method==0 || ovlp==0)
    Zci = 1 - alpha^(1/(L-1));
else
    % win = hamming(nfft);
    win = nuttallwin(nfft);
    % win = hann(nfft);
    den = sum(win.^2);
    num = 0;
    ww = 0;
    if (ovlp>0.75)
        for j=1:4
            for i=1:floor(nfft*(1-j*(1-ovlp)))
                num = num + win(i)*win(i+floor((1-ovlp)*nfft));
            end
            ww = ww + (num/den)^2;
        end
        w = 1/(1+(2*ww));
    else
        for i=1:floor(nfft*ovlp)
            num = num + win(i)*win(i+(1-ovlp)*nfft);
        end
        ww = (num/den)^2;
        w = 1/(1+(2*ww));
    end
    Lstar = floor((L-1)/(1-ovlp))+1;
    Zci = 1 - 0.05^(1/(w*Lstar-1));
end
