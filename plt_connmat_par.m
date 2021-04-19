function plt_connmat_par(titl,par,connmat_part)   % args: connect_matrix, component
    figure;
    set(gcf,'color','w');
    for m=1:size(connmat_part,1)
        max_lw = 6;    min_lw = 0.5;
        max_ns = 13;    min_ns = 4;
        k = size(connmat_part{m,1},1);
        H = zeros(k,length(connmat_part{m,1}.G{1,1}.Edges.Weight));
        N = zeros(k,size(connmat_part{m,1}.G{1,1}.Nodes,1));
        for j=1:k
            a = connmat_part{m,1}.G{j,1}.Edges.Weight';
            a = [a, NaN(1, length(H(j,:)) - length(a))];
            H(j,:) = a;
            for i=1:size(connmat_part{m,1}.G{j,1}.Nodes,1)
                N(j,i) = sum(connmat_part{m,1}.G{j,1}.Edges.Weight(outedges(connmat_part{m,1}.G{j,1},connmat_part{m,1}.G{j,1}.Nodes.Name(i))));
            end
        end
        maxv = max(H,[],'all');
        minv = min(H,[],'all');
        maxn = max(N,[],'all');
        minn = min(N,[],'all');
        for j=1:k
            subplot(size(connmat_part,1),k,j+3*(m-1));
            G = connmat_part{m,1}.G{j,1};
            h = plot(G,'Layout','subspace','EdgeLabel',G.Edges.Weight);
            set(gca,'XColor', 'none','YColor','none');
            % Position of nodes
            h.XData = [-1,0,1,0.5];
            h.YData = [0,0.8,0.6,-1];
            % Change width of Edges
            for i=1:size(G.Edges,1)
                path = G.Edges.EndNodes(i,:);
                linw = G.Edges.Weight(i);
                lw = (linw-minv)/(maxv-minv) * (max_lw-min_lw) + min_lw;
                highlight(h,path,'LineWidth',lw);
            end
            % Change size of Nodes
            for i=1:size(G.Nodes,1)
                wd = N(j,i);
                ns = (wd-minn)/(maxn-minn) * (max_ns-min_ns) + min_ns;
                highlight(h,G.Nodes.Name(i),'MarkerSize',ns);
            end
            if (par>10); part=par-10;str="Amputee: ";else ;part=par; str="Participant: ";end
            title({titl+" - ForceLevel: "+int2str(m),str + int2str(part)},{"Connectivity Matrix - Component: "+int2str(j)});
        end
    end
end