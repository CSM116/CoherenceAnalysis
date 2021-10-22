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
%         r = unifrnd(max(mn),max(mn)+0.1,1,indf-ind);
        r = unifrnd(0.9,1,1,indf-ind);
        w0(ind+1:ind+length(r),j) = r;
    end
    % Run NNMF
    dprev = 1;
    for i=1:1
        opt = statset('MaxIter',100);
        [Wtmp,Htmp,dnow] = nnmf(coh_tmp,k,'W0',w0,'H0',h0,'options',opt,'algorithm','mult');
        if (dnow<dprev)
            W = Wtmp;
            H = Htmp;
            dprev = dnow;
%             disp(dnow);
        end
    end
    % Rearrange
    [W_temp,H_temp] = rearrange(k,unit,freqs,W,H);
    % Plot NNMF
	if (a==1)
        newcolors = {[0.75 0.17 0.17],...
                     [0.80 0.42 0.00],...
                     [0.37 0.15 0.70],...
                     [0.15 0.55 0.10]};
        color = newcolors{w};                
        nexttile(w+(m-1)*numGestures);
%         subplot(fl,numGestures,w+(m-1)*numGestures);
        h = plot(nf,W_temp.');
        set(h,{'LineWidth'},{2.5;2.0;1.75});
        [col1, col2, col3] = deal(color,color+0.25,color+0.45);
        col1(col1>=1)=1; 
        col2(col2>=1)=1;
        col3(col3>=1)=1;
        set(h,{'color'},{col1;col2;col3});
        if (w~=1)
            set(gca, 'box','off','YTickLabel',[],'YTick',[]);
        else
            set(gca, 'box','off');
            xlabel('Frequency [Hz]');
        end
        if (titl=="Amputees");limy = 0.925;else; limy = 0.8;end
        set(gca,'xlim',[0 35],'ylim',[0 limy]);
%         title({titl+" - ForceLev. "+int2str(m),GestList(w)});
        title({titl,GestList(w)});
        legend(["1-5 Hz","5-12 Hz","12-40 Hz"]);
        set(gcf,'color','w'); % Set background colour to white
    end
end