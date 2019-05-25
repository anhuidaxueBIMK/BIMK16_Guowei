function [Weight] = Compute_Weight( param )
%COMPUTE_WEIGHT 此处显示有关此函数的摘要
%   此处显示详细说明
X = param.X;
Y = param.Y;
m = size(X,2);
Weight = zeros(1,m);
feature = zeros(1,m);
%pairwiseFeature = param.pairwiseFeature;
% X = param.X;

p=10;
for i=1:m
    %   Cdim = train_feature(:,1);
    C_feature = feature;
    C_feature(i) = 1;
    %   A = sortrows([Cdim,train_label]);
    ndcg1 = Compute_ndcg(X*C_feature',Y,p);
    C_feature(i) = -1;
    ndcg2 = Compute_ndcg(X*C_feature',Y,p);
    %     B =  sortrows([Cdim,train_label],-1);
    %    ndcg2 = Compute_ndcg(B(:,2),Y,p);
    Weight(i) = max(ndcg1,ndcg2);
    %  Weight(i)  = ndcg1;
end
 QR= sum(abs(X));
 Weight(QR==0)=0;
Weight  =Weight./max(Weight);

%  Graph = [Weight',(1:m)'];
%  S = [];
%  c =0.01;%平衡系数
%  for i=1:m
%      Graph = sortrows(Graph,-1);
%         Vk = Graph(1,2);
%         for j=1:size(Graph,1)
%             if(j~=Vk)
%             Weight(j) =   Weight(j) - c*CorrP(Vk,j);
%             end
%         end
%         S = [S,Vk];
%         Graph = Graph(2:end,:);
%  end


end


function ndcg = Compute_ndcg(Y,Yt,p)
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
    best_dcg(1) = (2^best_level(1) - 1);
    for j=2:q
        best_dcg(j) = (2^best_level(j)-1)/log2(j+1) +best_dcg(j-1);
    end
    
    if( best_dcg(1) == 0)
        continue;
    else
        [~,ind2] = sort(-Y(ind));
        
        current_level = conv(Yt{i}(ind2(1:q))+1);
        current_dcg =zeros(1,q);
        current_dcg(1) = (2^current_level(1) - 1);
        for j=2:q
            current_dcg(j) = (2^current_level(j)-1)/log2(j+1) +current_dcg(j-1);
        end
        ndcg(i,:) =current_dcg./best_dcg;
    end
end;
ndcg=mean(mean(ndcg));
end
