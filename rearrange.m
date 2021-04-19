function [W_temp,H_temp] = rearrange(k,unit,freqs,W,H)
            a = [0 freqs];
            indic = zeros(1,k);
            W_temp = zeros(size(W));
            H_temp = zeros(size(H));
            for ii=1:length(a)-2
                re_ind1 = ceil(unit*a(ii))+1;
                re_ind2 = ceil(unit*a(ii+1));
                fin_maxin = find(W(re_ind1:re_ind2,:)==max(W(re_ind1:re_ind2,:),[],'all'));
                if (fin_maxin<=(re_ind2-re_ind1)+1)
                    W_temp(:,ii) = W(:,1);
                    H_temp(ii,:) = H(1,:);
                    indic(1) = 1;
                elseif (fin_maxin>((re_ind2-re_ind1)+1)*2)
                    W_temp(:,ii) = W(:,3);
                    H_temp(ii,:) = H(3,:);
                    indic(3) = 1;
                else
                    W_temp(:,ii) = W(:,2);
                    H_temp(ii,:) = H(2,:);
                    indic(2) = 1;
                end  
            end
            W_temp(:,3) = W(:,find(indic==0));
            H_temp(3,:) = H(find(indic==0),:);
end