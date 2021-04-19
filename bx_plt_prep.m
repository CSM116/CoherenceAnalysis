function [bx_CC,bx_NS,bx_ED,bx_GE] = bx_plt_prep(numCommGes,from,fin,adjmat_part)
    bx_CC = zeros(0,0,0,0);
    bx_NS = zeros(0,0,0,0);
    bx_ED = zeros(0,0,0,0);
    bx_GE = zeros(0,0,0,0);
    for m=from:fin      
        for w=1:numCommGes
            bx_CC(:,:,w,m) = cell2mat(adjmat_part{1,m}{1,w}.CC)';
            bx_NS(:,:,w,m) = cell2mat(adjmat_part{1,m}{1,w}.NS)';
            bx_GE(:,:,w,m) = adjmat_part{1,m}{1,w}.GE';
            dimensions = cell2mat(cellfun(@size,adjmat_part{1,m}{1,w}.ED,'UniformOutput',false));
            tmp = 0;
            for jj=1:size(dimensions,1)
                if (size(dimensions(jj,:),2))
                    tmp = 1;
                end
            end
            if (length(unique(dimensions))>2 || tmp == 1)
                maxdim = 6;
                for ii=1:size((adjmat_part{1,m}{1,w}.ED),1)
                    if(max(size(adjmat_part{1,m}{1,w}.ED{ii,1}))<maxdim)
                        adjmat_part{1,m}{1,w}.ED{ii,1}(max(size(adjmat_part{1,m}{1,w}.ED{ii,1}))+1:maxdim) = NaN;
                    end
                end
            end
            bx_ED(:,:,w,m) = cell2mat(adjmat_part{1,m}{1,w}.ED)';
        end
    end
end