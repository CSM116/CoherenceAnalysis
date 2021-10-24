function [delta, probs, stats] = bx_plt2(bx_XX,i,a,type)
    delta = zeros(1,6);
    probs = zeros(1,6);
    tmp = squeeze(bx_XX(:,i,:));
    newcolors = [0.75 0.17 0.17;0.80 0.42 0.00;0.37 0.15 0.70;0.15 0.55 0.10];
    boxplot(tmp,'whisker', inf, 'Colors', newcolors);
    %Plot individual points for amputee participants
    %{
    hold on
    for i=1:size(tmp,2)
        scatter(i*ones(size(tmp,1),1),tmp(:,i),40,1.2*newcolors(i,:),'filled','MarkerEdgeColor',0.8*newcolors(i,:));
        hold on;
    end
    %}
    set(findobj(gca,'type','line'),'linew',1.1)
    med_lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(med_lines, 'linew', 1.75);
    set(gca,'TickLength',[0 0])
    bajo = min(bx_XX,[],'all');
    alto = max(bx_XX,[],'all');
    diff = abs(alto-bajo);
    if(alto<=0.6);add = 0.05;else; add=0;end
%     ylim([bajo-diff*0.17 alto+diff*0.38+add]);   
    ylim([bajo-diff*0.17 0.49]);
    hold on;
    yspace1 = diff*0.05; yspace2 = yspace1*0.6;
    % Condition for differential Plot
    if (a==1)
        % t-test
        [h,p] = ttest(tmp);
        txt_bajo = bajo-diff*0.07;
        text(1.7,txt_bajo-diff*0.035,'One-sample t-test p-values','FontSize',9,'fontweight', 'bold');
        for ii=1:size(bx_XX,3)
            if (p(ii)<0.0027);col='red'; elseif(p(ii)<0.05);col=[0.9500 0.5250 0.0980]; else; col ='w'; end
            txt = [num2str(p(ii), 2)];
            text(ii-0.1,txt_bajo, txt, 'FontSize',9, 'Color',col);
        end
    end
    lim = size(tmp,2);
    off = 0;
    jj = 1;
    for j=1:lim-1
        ylimpos = max(tmp(:,lim-j:lim),[],'all');
        for ii=lim-j:lim-1
            if (type=='t')
%                 'Alpha',0.0027
                [h,p,ci,st] = ttest2(tmp(:,lim-j),tmp(:,ii+1),'Alpha',0.05);
                delta(jj) = -100*((mean(tmp(:,lim-j))-mean(tmp(:,ii+1)))/mean(tmp(:,lim-j)));
                probs(jj) = p;
                stats(jj) = st;
            elseif (type=='w')
                [p,h,st] = ranksum(tmp(:,lim-j),tmp(:,ii+1));
                delta(jj) = -100*((mean(tmp(:,lim-j))-mean(tmp(:,ii+1)))/mean(tmp(:,lim-j)));
                probs(jj) = p;
                stats(jj) = st;
            end
            jj = jj+1;
            if (p<0.0027); col = 'red'; elseif (p<0.05); col=[0.9500 0.5250 0.0980]; else; col ='w'; end
            ylimpos = ylimpos+yspace1 + off;
            bx_connector([lim-j ii+1],[ylimpos ylimpos],col);
            ypos = ylimpos + yspace2;
%             txt = ['p: ' num2str(p, 2)];
%             text((((lim-j)+(ii+1))/2)-0.25,ypos,txt,'FontSize',9,'Color',col);
            txt = "*";
            if (h); text((((lim-j)+(ii+1))/2),ypos,txt,'FontSize',11,'Color',col);end
            off = off + 0.0015;
        end
    end
    % Text two sample t-test
    if (type=='t')
        texto = 'Two-sample t-test';
    elseif (type=='w')
        texto = 'Wilcoxon ranksum';
    end
%     text((lim/2)-0.3,ypos+yspace1,texto,'FontSize',9,'fontweight','bold');                
    hold off;
end