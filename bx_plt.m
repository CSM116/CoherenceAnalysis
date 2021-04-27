function [] = bx_plt(bx_metrics,k,VarNames,GestList,titl,a,type)    % a: 0-normal 1-differential / type: t-ttest w-wilcoxon
    bx_CC = bx_metrics{1};
    bx_NS = bx_metrics{2};
    bx_ED = bx_metrics{3};
    if(length(bx_metrics)>3);bx_GE = bx_metrics{4};end
    for j=1:length(bx_metrics)
        figure;
        for i=1:k   % Loop through 3 components
            subplot(1,k,i);
            if (j==1)
                bx_plt2(bx_CC,i,a,type);
                ti = VarNames(4);
            elseif (j==2)
                bx_plt2(bx_NS,i,a,type);
                ti = VarNames(7);
            elseif (j==3)
                bx_plt2(bx_ED,i,a,type);
                ti = VarNames(8);
            else
                bx_plt2(bx_GE,i,a,type);
                ti = VarNames(5);
            end
            xticklabels(GestList(1,1:size(bx_CC,3)));
            if (a==1); titl = 'Differential Plot'; end
    %         xtickangle(45);
            title({titl, ti + " - k" + i});
        end
    end
end