function [bx_CC,bx_NS,bx_ED,bx_GE] = bx_plt_prep2(numCommGes,k,adjmat_ges,adjmatAmp_ges)
    tres=0;
    if (nargin<4)
        adjmatAmp_ges=cell(size(adjmat_ges));
        tres=1;
    end
    bx_CC = zeros(0,0,0);
    bx_NS = zeros(0,0,0);
    bx_ED = zeros(0,0,0);
    bx_GE = zeros(0,0,0);
    for w=1:numCommGes
        dimensions = cell2mat(cellfun(@size,adjmat_ges{1,w}.ED,'UniformOutput',false));
        tmp = 0;
        for jj=1:size(dimensions,1)
            if (size(dimensions(jj,:),2))
                tmp = 1;
            end
        end
        if (length(unique(dimensions))>2 || tmp == 1)
            maxdim = 6;
            for ii=1:size((adjmat_ges{1,w}.ED),1)
                if(max(size(adjmat_ges{1,w}.ED{ii,1}))<maxdim)
                    adjmat_ges{1,w}.ED{ii,1}(max(size(adjmat_ges{1,w}.ED{ii,1}))+1:maxdim) = NaN;
                end
            end
        end
        if (tres==1)
            CC = zeros(size(cell2mat(adjmat_ges{1,w}.CC))); CC = mat2cell(CC,ones(1,k));
            NS = zeros(size(cell2mat(adjmat_ges{1,w}.NS))); NS = mat2cell(NS,ones(1,k));
            ED = zeros(size(cell2mat(adjmat_ges{1,w}.ED))); ED = mat2cell(ED,ones(1,k));
            GE = zeros(size(adjmat_ges{1,w}.GE)); 
            adjmatAmp_ges{1,w} = table(CC,NS,ED,GE); 
        end
        bx_CC(:,:,w) = cell2mat(adjmat_ges{1,w}.CC)'-cell2mat(adjmatAmp_ges{1,w}.CC)';
        bx_NS(:,:,w) = cell2mat(adjmat_ges{1,w}.NS)'-cell2mat(adjmatAmp_ges{1,w}.NS)';
        bx_GE(:,:,w) = adjmat_ges{1,w}.GE'-adjmatAmp_ges{1,w}.GE';
        bx_ED(:,:,w) = cell2mat(adjmat_ges{1,w}.ED)'-cell2mat(adjmatAmp_ges{1,w}.ED)';
    end
end