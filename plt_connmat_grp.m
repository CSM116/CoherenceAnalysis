function plt_connmat_grp(titl_ges,connmat_grp)   % args: connect_matrix, component
%     figure;
    for m=1:size(connmat_grp,1)                 % Levels of Force
        max_lw = 12;     min_lw = 1.5;
        max_ns = 22;    min_ns = 5;
        k = size(connmat_grp{m,1},1);
        H = zeros(k,6);
        N = zeros(k,4);
        for j=1:k
            a = connmat_grp{m,1}.G{j,1}.Edges.Weight';
            a = [a, NaN(1, length(H(j,:)) - length(a))];
            H(j,:) = a;
            for i=1:size(connmat_grp{m,1}.G{j,1}.Nodes,1)
                N(j,i) = sum(connmat_grp{m,1}.G{j,1}.Edges.Weight(outedges(connmat_grp{m,1}.G{j,1},connmat_grp{m,1}.G{j,1}.Nodes.Name(i))));
            end
        end
        maxv = 1; %max(H,[],'all');
        minv = 0; %min(H,[],'all');
        maxn = 2; %max(N,[],'all');
        minn = 0; %min(N,[],'all');
        tiledlayout(k,size(connmat_grp,1),'TileSpacing','none');
        for j=1:k
            nexttile(j+k*(m-1));
%             subplot(size(connmat_grp,1),k,j+k*(m-1));
            % Load forearm cross-section image
%             %{
            axis([-2 2 -2 2]);
            I = imread('Cross-Forearm.png');
            image('XData',[-1.75 1.75],'YData',[1.75 -1.75],'CData',I) ;
%             alpha(0.45);
%             colormap(gcf, gray(256));
            hold on;
            %}
            % Plot Graph
            G = connmat_grp{m,1}.G{j,1};
            h = plot(G,'Layout','subspace','EdgeCData', G.Edges.Weight,...
                'EdgeLabel',G.Edges.Weight,'NodeLabel',{},'EdgeLabel',{},'Edgealpha',0.85);
%                 'Edgecolor','#0072BD','NodeColor','#0072BD');
            cb = colorbar;
            set(cb,'visible','off');
            cb.Layout.Tile = 'east';
            J = customcolormap([0 1], [1 1 1; 0 0 1]);
            colormap(flip(J));
            caxis([0.1 0.55]);
%             cb.Ticks = linspace(0.05,0.65,7);
%             cb.TickLabels = num2cell(0.05:0.1:0.65);
            if (j==k)
                cb = colorbar;
                set(cb,'visible','on','Fontsize',11);
                cb.Layout.Tile = 'east';
%                 cb.Ticks = linspace(0.05,0.65,7);
%                 cb.TickLabels = num2cell(0.05:0.1:0.65);
                caxis([0.1 0.55]);
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
%                 if linw>0;highlight(h,path,'EdgeColor','b');end
%                 if linw>0.5;highlight(h,path,'EdgeColor','b');end
%                 if linw<=0.2;highlight(h,path,'EdgeColor','c');end
            end
            % Change size of Nodes
            for i=1:size(G.Nodes,1)
                wd = N(j,i);
                ns = (wd-minn)/(maxn-minn) * (max_ns-min_ns) + min_ns;
                highlight(h,G.Nodes.Name(i),'MarkerSize',ns);
                if ns>=0;highlight(h,G.Nodes.Name(i),'NodeColor','b');end
%                 if ns>=12;highlight(h,G.Nodes.Name(i),'NodeColor','b');end
%                 if ns<=9; highlight(h,G.Nodes.Name(i),'NodeColor','c');end
            end
%             title({titl+" - ForceLevel: "+int2str(m),"Connectivity Matrix - Component: "+int2str(j)});
            component = ["1-5 Hz","5-12 Hz","12-35 Hz"];
            ylabel(strcat("Component: ", component(j)),'FontSize',13.5,'Fontweight', 'bold');
%             title({titl,component(j)});
        end
    end
end