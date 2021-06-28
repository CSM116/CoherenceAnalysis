function [adjmat] = connect_mat(m,w,i,k,s,t,nnmf_x,th_op,th_val,varnames)
    tbl0 = cell2table(cell(0,length(varnames)), 'VariableNames', varnames);
    for j=1:k
        try
            if (i==0)
                G = graph(s,t,nnmf_x{m,w}.H(j,:));
            else
                G = graph(s,t,nnmf_x{1,i}(m,w).H(j,:));
            end
        catch
            if (i==0)
                G = graph(s,t,nnmf_x{m,w});
            else
                G = graph(s,t,nnmf_x{1,i}{m,w});
            end
        end
        ED = G.Edges.Weight';
        % Threshold
%         %{
        if th_op==0         % for NO thresholding
            ind = 1;
        elseif th_op==1     % for std thresholding
%             ind = ED < (median(ED)+std(ED)*th_std) & ED > (median(ED)-std(ED)*th_std);	% index of values inside 1 std from median
            ind = ED > (max(ED)-std(ED)*th_val);	% index of values higher than 1 std from median

        elseif th_op==2     % for percentage value thresholding
            ind = zeros(1,6);
            [B,I] = maxk(ED,round(6*th_val/100));   
            ind(I) = 1;
        else                % for absolute value thresholding
            ind = ED > th_val;  
        end
        G = rmedge(G,find(~ind));                                               % remove edges outside range
        ED = G.Edges.Weight';
%         %}
        A = adjacency(G,'weighted');
        BC = betweenness_wei(A);
        CC = clustering_coef_wu(A);
        GE = efficiency_wei(A,0);
        ME = median(G.Edges.Weight);
        SD = std(G.Edges.Weight);
        NS = full(sum(A));
        tbl0 = [tbl0;{G,A,{BC},{full(CC')},GE,ME,{NS},{ED},SD}];
    end
    adjmat = tbl0;
end