function [ istanceAim2 ] = istanceAim2Compute( param,Population  )
%ISTANCEAIM2COMPUTE 此处显示有关此函数的摘要
%   此处显示详细说明
X = param.X;
Y = param.Y;
% sparseX = sparse(X);
% pairwiseLabel = param.pairwiseLabel;
pairwiseFeature = param.pairwiseFeature;
PopulationLength = size(Population,1);
istanceAim2 = zeros(PopulationLength,1);

for i=1:PopulationLength
    
    currentPop = logical(Population(i,:));
    if sum(currentPop)==0
        istanceAim2(i)=1;
        continue; 
    end
    currentpairwiseFeature = pairwiseFeature(currentPop,:);
    currentpairwiseLabel = param.pairwiseLabel(currentPop);
    currentpairwiseFeature = sparse(currentpairwiseFeature);
  %  model = train(currentpairwiseLabel,currentpairwiseFeature,['-s 2 -c ' num2str(param.BestParameterC) ' -q']);
     model = train(currentpairwiseLabel,currentpairwiseFeature,'-s 2 -c 1 -q');
   %    model = fitcsvm(currentpairwiseFeature,currentpairwiseLabel,'-s 3 -t 2 -c 2.2 -g 2.8 -p 0.01');
    %model2 = train(currentpairwiseLabel,currentpairwiseFeature,['-s 2 -c 1 -q']);
%           [regressWX]= predict(param.Y_,sparseX,model,'-q');
  W = model.w;
   NDCG =  compute_ndcg4(W*X',Y,10);
%     NDCG =  compute_ndcg4(regressWX,Y,10);
    istanceAim2(i) =1- NDCG; 
end

end

