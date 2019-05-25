function [index] = selectedExplanation( param)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
Xt = param.Xt;
Yt = param.Yt;
pairwiseLabel = param.pairwiseLabel;
pairwiseFeature = param.pairwiseFeature;
poplength = size(Xt,2);
intanceFront1_PopSet = param.intanceFront1_Pop;
featureFront1_PopSet = param.featureFront1_Pop;
intanceFront1_Pop = intanceFront1_PopSet(:,3:end);
featureFront1_Pop = featureFront1_PopSet(:,3:end);
intanceFront1_PopLength = size(intanceFront1_PopSet,1);
featureFront1_PopLength = size(featureFront1_PopSet,1);
allNDCG = zeros(intanceFront1_PopLength*featureFront1_PopLength,1);

for j=1:featureFront1_PopLength
    currentFeatureFront1_Pop = logical(featureFront1_Pop(j,:));
    for i=1 : intanceFront1_PopLength
        currentIntanceFront1_Pop = logical(intanceFront1_Pop(i,:));
        
        feature = sparse(pairwiseFeature(currentIntanceFront1_Pop,currentFeatureFront1_Pop));
        label   = pairwiseLabel(currentIntanceFront1_Pop);
        model = train(label,feature,'-s 2');
        W = model.w;
        complete_W = zeros(1,poplength);
        complete_W(currentFeatureFront1_Pop==1) = W;
        NDCG =  compute_ndcg3(complete_W*Xt',Yt,10);
        allNDCG((i-1)*featureFront1_PopLength+j) = NDCG;
    end
end
% [maxNDCG] = max(max(allNDCG));
% [sitei,sitej] = find(allNDCG==maxNDCG);
% bestFeatureNum = featureFront1_PopSet(sitej,1);
% bestIstanceNum = intanceFront1_PopSet(sitei,1);
[~,index] = sort(allNDCG);
bestSite = zeros(length(index),2);
for i = 1:length(index)
    a = index(i);
    bestSite(i,1) = ceil(a/intanceFront1_PopLength);
    bestSite(i,2) = rem(a,intanceFront1_PopLength);
end
end

