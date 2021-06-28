function [delta, probs] = bx_plt(bx_metrics,k,VarNames,GestList,titl,a,type)    % a: 0-normal 1-differential / type: t-ttest w-wilcoxon
    bx_CC = bx_metrics{1};
    bx_NS = bx_metrics{2};
    bx_ED = bx_metrics{3};
    delta = cell(k,length(bx_metrics));
    probs = cell(k,length(bx_metrics));
    if(length(bx_metrics)>3);bx_GE = bx_metrics{4};end
    for j=1:length(bx_metrics)
%         figure;
        figure;tiledlayout(k,1,'Padding','compact','TileSpacing','tight');
        for i=1:k   % Loop through 3 components
%             subplot(1,k,i);
            nexttile(i)
            if (j==1)
                [delta{i,j},probs{i,j}] = bx_plt2(bx_CC,i,a,type);
                ti = VarNames(4);
            elseif (j==2)
                [delta{i,j},probs{i,j}] = bx_plt2(bx_NS,i,a,type);
                ti = VarNames(7);
            elseif (j==3)
                [delta{i,j},probs{i,j}] = bx_plt2(bx_ED,i,a,type);
                ti = VarNames(8);
            else
                [delta{i,j},probs{i,j}] = bx_plt2(bx_GE,i,a,type);
                ti = VarNames(5);
            end
            if(i==3)
                set(gca, 'box','off');
            else
                set(gca, 'box','off','XTickLabel',[],'XTick',[]);
            end
            % Generate x tick labels
            gestos = upper(convertStringsToChars(GestList));
            gestos1 = cell(1,size(bx_CC,3));
            for ge=1:size(bx_CC,3)
                gestos1{ge} = gestos{ge}(1:3);
            end
            xticklabels(gestos1);
%             xticklabels(GestList(1,1:size(bx_CC,3)));  

            if (a==1); titl = 'Differential Plot'; end
%             xtickangle(45);
%             title({titl, ti + " - k" + i});
            component = ["1-5 Hz","5-15 Hz","15-35 Hz"];
            ylabel(strcat("Component: ", component(i)),'FontSize',13.5,'Fontweight', 'bold');
        end
        set(gcf,'color','w'); % Set background colour to white
    end
end