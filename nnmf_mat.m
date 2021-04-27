function [H_temp,W_temp] = nnmf_mat(m,w,i,k,coh_x,freqs,nf,fl,numGestures,titl,GestList,a)
    if (i==0)
        coh_tmp = coh_x{m,w};
        if(size(coh_x,2)==1)
            numGestures = 1;
            GestList = {'All'};
        end
    else
        coh_tmp = coh_x{1,i}{m,w};
    end
    % NNMF Initialization
    unit = length(coh_tmp)/500;
    w0 = unifrnd(0,0.0001,length(coh_tmp),k);
    h0 = unifrnd(0,0.05,k,size(coh_tmp,2));
    indf = 0;
    for j=1:k
        ind = indf;
        indf = ceil(unit*freqs(j));
        mn = median(coh_tmp(ind+1:indf,:),2);
        r = unifrnd(max(mn),max(mn)+0.1,1,indf-ind);
%                 r = unifrnd(0.9,1,1,indf-ind);
        w0(ind+1:ind+length(r),j) = r;
    end
    % Run NNMF
    opt = statset('MaxIter',100);
    [W,H] = nnmf(coh_tmp,k,'W0',w0,'H0',h0,'options',opt,'algorithm','mult');
    % Rearrange
    [W_temp,H_temp] = rearrange(k,unit,freqs,W,H);
    % Plot NNMF
	if (a==1)
        subplot(fl,numGestures,w+(m-1)*numGestures);
        plot(nf,W_temp,'linew',1.1);
        set(gca,'xlim',[0 35]);
        title({titl+" - ForceLev. "+int2str(m),"Gesture: "+GestList(w)});
        xlabel('Frequency [Hz]');
    end
end