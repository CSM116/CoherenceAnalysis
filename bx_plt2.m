function [] = bx_plt2(bx_XX,i,a,type)
    tmp = squeeze(bx_XX(:,i,:));
    boxplot(tmp,'whisker', inf);
    set(gca,'TickLength',[0 0])
    bajo = min(bx_XX,[],'all');
    alto = max(bx_XX,[],'all');
    diff = abs(alto-bajo);
    if(alto<=0.6);add = 0.05;else; add=0;end
    ylim([bajo-diff*0.17 alto+diff*0.38+add]);   
    hold on;
    yspace1 = diff*0.05; yspace2 = yspace1*0.6;
    % Condition for differential Plot
    if (a==1)
        % t-test
        [h,p] = ttest(tmp);
        txt_bajo = bajo-diff*0.07;
        text(1.7,txt_bajo-diff*0.035,'One-sample t-test p-values','FontSize',9,'fontweight', 'bold');
        for ii=1:size(bx_XX,3)
            if (p(ii)<0.05); col ='red'; else; col ='black'; end
            txt = [num2str(p(ii), 2)];
            text(ii-0.1,txt_bajo, txt, 'FontSize',9, 'Color',col);
        end
    end
    lim = size(tmp,2);
    off = 0;
    for j=1:lim-1
        ylimpos = max(tmp(:,lim-j:lim),[],'all');
        for ii=lim-j:lim-1
            if (type=='t')
                [h,p] = ttest2(tmp(:,lim-j),tmp(:,ii+1));
            elseif (type=='w')
                [p,h] = ranksum(tmp(:,lim-j),tmp(:,ii+1));
            end
            if (h); col ='red'; else; col ='black'; end
            ylimpos = ylimpos+yspace1 + off;
            bx_connector([lim-j ii+1],[ylimpos ylimpos],col);
            ypos = ylimpos + yspace2;
            txt = ['p: ' num2str(p, 2)];
            text((((lim-j)+(ii+1))/2)-0.25,ypos,txt,'FontSize',9,'Color',col);
            off = off + 0.0015;
        end
    end
    % Text two sample t-test
    if (type=='t')
        texto = 'Two-sample t-test p-values';
    elseif (type=='w')
        texto = 'Wilcoxon ranksum p-values';
    end
    text((lim/2)-0.3,ypos+yspace1,texto,'FontSize',9,'fontweight','bold');                
    hold off;
end