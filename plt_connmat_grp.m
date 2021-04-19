function plt_connmat_grp(titl,connmat_grp)   % args: connect_matrix, component
    figure;
    for m=1:size(connmat_grp,1)                 % Levels of Force
        max_lw = 8;     min_lw = 0.01;
        max_ns = 18;    min_ns = 0.5;
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
        for j=1:k
            subplot(size(connmat_grp,1),k,j+k*(m-1));
            % Load forearm cross-section image
            %{
            axis([-2 2 -2 2]);
            I = imread('cross_forearm.png');
            I = I(:,:,1);
            image('XData',[-2 2],'YData',[2 -2],'CData',I) ;
            colormap(gcf, gray(256));
            hold on;
            %}
            % Plot Graph
            G = connmat_grp{m,1}.G{j,1};
            h = plot(G,'Layout','subspace','EdgeLabel',G.Edges.Weight);
            set(gca,'XColor', 'none','YColor','none');
            % Position of nodes
            h.XData = [-1,0,1,0.5];
            h.YData = [0,0.8,0.6,-1];
            hold off;
            % Change width of Edges
            for i=1:size(G.Edges,1)
                path = G.Edges.EndNodes(i,:);
                linw = G.Edges.Weight(i);
                lw = (linw-minv)/(maxv-minv) * (max_lw-min_lw) + min_lw;
                highlight(h,path,'LineWidth',lw);
                if linw>0.5;highlight(h,path,'EdgeColor','k');end
                if linw<=0.2;highlight(h,path,'EdgeColor','r');end
            end 
            % Change size of Nodes
            for i=1:size(G.Nodes,1)
                wd = N(j,i);
                ns = (wd-minn)/(maxn-minn) * (max_ns-min_ns) + min_ns;
                highlight(h,G.Nodes.Name(i),'MarkerSize',ns);
            end
            title({titl+" - ForceLevel: "+int2str(m),"Connectivity Matrix - Component: "+int2str(j)});
        end
    end
end