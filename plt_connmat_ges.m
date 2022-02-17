function plt_connmat_ges(gest_list,connmat_grp)   % args: connect_matrix, component
    figure;
    max_lw = 12;     min_lw = 1.5;
    max_ns = 22;    min_ns = 5;
    k = size(connmat_grp{1,1},1);
    numges = size(connmat_grp,2);
    tiledlayout(k,size(connmat_grp,2),'Padding','compact','TileSpacing','none');
    
    for j=1:k                                       % loop through number of components
        for m=1:numges                              % loop through gestures          
            
            H = zeros(numges,6);                         % Values of edges
            N = zeros(numges,4);                         % Values of nodes
            a = connmat_grp{1,m}.G{j,1}.Edges.Weight';
            a = [a, NaN(1, length(H(j,:)) - length(a))];
            H(m,:) = a;
            for i=1:size(connmat_grp{1,m}.G{j,1}.Nodes,1)
                N(m,i) = sum(connmat_grp{1,m}.G{j,1}.Edges.Weight(outedges(connmat_grp{1,m}.G{j,1},connmat_grp{1,m}.G{j,1}.Nodes.Name(i))));
            end
            maxv = 1; %max(H,[],'all');
            minv = 0; %min(H,[],'all');
            maxn = 2; %max(N,[],'all');
            minn = 0; %min(N,[],'all');
            
            nexttile(m+numges*(j-1));
            % Load forearm cross-section image
            axis([-2 2 -2 2]);
            I = imread('Cross-Forearm.png');
            image('XData',[-1.75 1.75],'YData',[1.75 -1.75],'CData',I);
            hold on;
            % Plot Graph
            G = connmat_grp{1,m}.G{j,1};
            h = plot(G,'Layout','subspace','EdgeCData', G.Edges.Weight,...
                'EdgeLabel',G.Edges.Weight,'NodeLabel',{},'EdgeLabel',{},'Edgealpha',0.85);
%                 'Edgecolor','#0072BD','NodeColor','#0072BD');

            if (j==k&&m==numges)
                cb = colorbar;
                set(cb,'visible','on','Fontsize',11);
                cb.Layout.Tile = 'east';
%                 J = customcolormap([0 0.15 0.35 0.65 1], [1,1,0.5; 0.9,0.9,0; 0.9,0.69,0.15; 0.35,0.75,0.45; 0,0,1]);
                J = customcolormap([0 1], [0.98,0.98,0.98; 0,0.3,1]);
                colormap(flip(J));
%                 colormap(flip(parula));
                cb.Ticks = linspace(0.05,0.65,7);
                cb.TickLabels = num2cell(0.05:0.1:0.65);
                caxis([0.1 0.65]);
            else
                colorbar('off');
            end
            set(gca,'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
            set(gca,'XColor','none','YColor','none');
            set(get(gca,'YLabel'),'visible','on','color','k');
            set(gcf,'color','w');
            % Position of nodes
            h.XData = [0.7,  1.8,  0.5,  -1.65];
            h.YData = [1.71,  0.35, -1.79,   0.85];
            hold off;
            % Change width of Edges
            for i=1:size(G.Edges,1)
                path = G.Edges.EndNodes(i,:);
                linw = G.Edges.Weight(i);
                lw = (linw-minv)/(maxv-minv) * (max_lw-min_lw) + min_lw;
                highlight(h,path,'LineWidth',lw);
            end
            % Change size of Nodes
            for i=1:size(G.Nodes,1)
                wd = N(m,i);
                ns = (wd-minn)/(maxn-minn) * (max_ns-min_ns) + min_ns;
                highlight(h,G.Nodes.Name(i),'MarkerSize',ns);
                if ns>=0;highlight(h,G.Nodes.Name(i),'NodeColor',[0,0.3,1]);end
            end
            if((m+numges*(j-1))==1||(m+numges*(j-1))==5||(m+numges*(j-1)==9))
%             title({titl+" - ForceLevel: "+int2str(m),"Connectivity Matrix - Component: "+int2str(j)});
                component = ["1-5 Hz","5-12 Hz","12-40 Hz"];
                ylabel(strcat("Component: ", component(j)),'FontSize',12,'Fontweight', 'bold');
%             title({titl,component(j)});
            end
            if((m+numges*(j-1))>8)
                title(gest_list((m+numges*(j-1))-8));
                set(get(gca,'title'),'FontSize',12,'Fontweight','bold','Position',...
                     [0.15, -2.75], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            end
        end
    end
end