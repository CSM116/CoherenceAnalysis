%% Clean Workspace
% clear;
% close all;
clc;
%% Global variables
GestList = ["Flexion" "Extension" "Pronation" "Supination"];
alpha = 0.05;
adj_alpha = 0.0028;
%% Titles
figure;
tiledlayout(1,1,'Padding','tight');
% Titles Horizontal
title_pos_x = 0.225;
title_pos_y = 0.90;
for i=1:4
    text(title_pos_x,title_pos_y,GestList(i),'fontweight','bold','fontsize',15);
    title_pos_x = title_pos_x + 0.2;
end
% Titles Vertical
title_pos_x = 0.025;
title_pos_y = 0.70;
for i=1:4
    text(title_pos_x,title_pos_y,GestList(i),'fontweight', 'bold','fontsize',15);
    title_pos_y = title_pos_y - 0.2;
end
xlim([0 1]);
ylim([0 1]);
hold on;
%% Box
plot([0 1],[1 1],'k');
plot([1 1],[0 1],'k');
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
%% Edges
% Edges x
for i=0.2:0.2:0.8
    plot([i i],[0 1],'k');
end
% Edges y
for i=0.2:0.2:0.8
    plot([0 1],[i i],'k');
end
%% Edges components - right
c = 0.6;
a = 0.2/3;
pos_x = 0.4;
pos_y = 0.8 - a;
for j=1:3
    for i=1:2
        plot([pos_x 1],[pos_y pos_y],'color',[c c c]);
        pos_y = pos_y - a;
    end
    pos_x = pos_x + 0.2;
    pos_y = pos_y - a;
end
%% Edges components - left
a = 0.2/3;
pos_x = 0.4;
pos_y = 0.6 - a;
for j=1:3
    for i=1:2
        plot([0.2 pos_x],[pos_y pos_y],'color',[c c c]);
        pos_y = pos_y - a;
    end
    pos_x = pos_x + 0.2;
    pos_y = pos_y - a;
end

%% Add Probabilities and Delta - Right
for kk=1:k  % index component
    % Last position
    x_pos = 0.81;
    y_pos = 0.365-(kk-1)*0.065;
    prb = probs{kk,1}(1);
    delt = delta{kk,1}(1);
    sta = stats{kk,1}(1).tstat;
    if (prb<alpha)
        if(prb<adj_alpha); col='red'; else; col='orange'; end
        if (delt<0); sym = '\color{red}\downarrow'; else ; sym = '\color{blue}\uparrow'; end
        str = strcat('\color{black}','\Delta:'," ",num2str(delt,3),"% ",sym);
    else
        col = 'black';
        str = '';
    end
    text(x_pos,y_pos,strcat('k',int2str(kk)," - ",'{\itt-val}:'," ",num2str(sta,3)," ",strcat('\color{',col,'}'),'{\itp}:'," ",num2str(prb,2)," ",str),'FontSize',13);
    % Middle position
    x_pos = x_pos - 0.2;
    y_pos = y_pos + 0.2;
    for i = 1:2
        prb = probs{kk,1}(i+1);
        delt = delta{kk,1}(i+1);
        sta = stats{kk,1}(i+1).tstat;
        if (prb<alpha)
            if(prb<adj_alpha); col='red'; else; col='orange'; end
            if (delt<0); sym = '\color{red}\downarrow'; else ; sym = '\color{blue}\uparrow'; end
            str = strcat('\color{black}','\Delta:'," ",num2str(delt,3),"% ",sym);
        else
            col = 'black';
            str = '';
        end
        text(x_pos,y_pos,strcat('k',int2str(kk)," - ",'{\itt-val}:'," ",num2str(sta,3)," ",strcat('\color{',col,'}'),'{\itp}:'," ",num2str(prb,2)," ",str),'FontSize',13);
        x_pos = x_pos + 0.2;
    end
    % First positions
    x_pos = x_pos - 0.6;
    y_pos = y_pos + 0.2;
    for i = 1:3
        prb = probs{kk,1}(i+3);
        delt = delta{kk,1}(i+3);
        sta = stats{kk,1}(i+3).tstat;
        if (prb<alpha)
            if(prb<adj_alpha); col='red'; else; col='orange'; end
            if (delt<0); sym = '\color{red}\downarrow'; else ; sym = '\color{blue}\uparrow'; end
            str = strcat('\color{black}','\Delta:'," ",num2str(delt,3),"% ",sym);
        else
            col = 'black';
            str = '';
        end
        text(x_pos,y_pos,strcat('k',int2str(kk)," - ",'{\itt-val}:'," ",num2str(sta,3)," ",strcat('\color{',col,'}'),'{\itp}:'," ",num2str(prb,2)," ",str),'FontSize',13);
        x_pos = x_pos + 0.2;
    end
end
%% Add Probabilities and Delta - Left
for kk=1:k  % index component
    % Last position
    x_pos = 0.61;
    y_pos = 0.165-(kk-1)*0.065;
    prb = probs{kk,1}(1);
    delt = -delta{kk,1}(1);
    sta = stats{kk,1}(1).tstat;
    if (prb<alpha)
        if(prb<adj_alpha); col='red'; else; col='orange'; end
        if (delt<0); sym = '\color{red}\downarrow'; else ; sym = '\color{blue}\uparrow'; end
        str = strcat('\color{black}','\Delta:'," ",num2str(delt,3),"% ",sym);
    else
        col = 'black';
        str = '';
    end
    text(x_pos,y_pos,strcat('k',int2str(kk)," - ",'{\itt-val}:'," ",num2str(sta,3)," ",strcat('\color{',col,'}'),'{\itp}:'," ",num2str(prb,2)," ",str),'FontSize',13);
    % Middle position
    x_pos = x_pos - 0.2;
    y_pos = y_pos + 0.2;
    for i = 1:2
        prb = probs{kk,1}(i+1);
        delt = -delta{kk,1}(i+1);
        sta = stats{kk,1}(i+1).tstat;
        if (prb<alpha)
            if(prb<adj_alpha); col='red'; else; col='orange'; end
            if (delt<0); sym = '\color{red}\downarrow'; else ; sym = '\color{blue}\uparrow'; end
            str = strcat('\color{black}','\Delta:'," ",num2str(delt,3),"% ",sym);
        else
            col = 'black';
            str = '';
        end
        text(x_pos,y_pos,strcat('k',int2str(kk)," - ",'{\itt-val}:'," ",num2str(sta,3)," ",strcat('\color{',col,'}'),'{\itp}:'," ",num2str(prb,2)," ",str),'FontSize',13);
        y_pos = y_pos - 0.2;
    end
    % First positions
    y_pos = y_pos + 0.6;
    x_pos = x_pos - 0.2;
    for i = 1:3
        prb = probs{kk,1}(i+3);
        delt = -delta{kk,1}(i+3);
        sta = stats{kk,1}(i+3).tstat;
        if (prb<alpha)
            if(prb<adj_alpha); col='red'; else; col='orange'; end
            if (delt<0); sym = '\color{red}\downarrow'; else ; sym = '\color{blue}\uparrow'; end
            str = strcat('\color{black}','\Delta:'," ",num2str(delt,3),"% ",sym);
        else
            col = 'black';
            str = '';
        end
        text(x_pos,y_pos,strcat('k',int2str(kk)," - ",'{\itt-val}:'," ",num2str(sta,3)," ",strcat('\color{',col,'}'),'{\itp}:'," ",num2str(prb,2)," ",str),'FontSize',13);
        y_pos = y_pos - 0.2;
    end
end
set(gcf,'color','w'); % Set background colour to white
