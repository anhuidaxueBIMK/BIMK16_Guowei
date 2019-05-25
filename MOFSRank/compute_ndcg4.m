function meanndcg = compute_ndcg4(Y,Yt,p)
conv = [0 1 3 7 15];
ind = 0;
ndcg = zeros(length(Yt),p);
for i=1:length(Yt)
    ind = ind(end)+[1:length(Yt{i})];
    if( p > length(ind))
        continue;
    end
    q = min(p,length(ind));
    %     disc = [1 log2(2:q)];
    [~,ind2] = sort(-Yt{i});
    % best_dcg = sum(conv(Yt{i}(ind2(1:q))+1) ./ disc) + eps;
    %     best_dcg = sum(conv(Yt{i}(ind2(1:q))+1) ./ disc);
    best_level = conv(Yt{i}(ind2(1:q))+1);
    best_dcg =zeros(1,q);
    %     best_dcg(1) = (2^best_level(1) - 1);
    best_dcg(1) = best_level(1);
    for j=2:q
        if j<3
            best_dcg(j) = best_dcg(j-1)+best_level(j);
        else
            best_dcg(j) = best_dcg(j-1)+ best_level(j)*(log(2)/log(j+1-1));
        end
    end
    
    if( best_dcg(1) == 0)
        continue;
    else
        [~,ind2] = sort(-Y(ind));
        
        current_level = conv(Yt{i}(ind2(1:q))+1);
        current_dcg =zeros(1,q);
        current_dcg(1) = current_level(1);
      
        for j=2:q
            if j<3
                current_dcg(j) = current_dcg(j-1)+current_level(j);
            else
                current_dcg(j) = current_dcg(j-1)+ current_level(j)*(log(2.0)/log(j+1-1));
            end
        end
        ndcg(i,:) =current_dcg./best_dcg;
    end
end;

ndcg10 = ndcg(:,10);

meanndcg = mean(ndcg10);
