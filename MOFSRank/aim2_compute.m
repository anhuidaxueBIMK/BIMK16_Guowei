function [ aim2 ] = aim2_compute(  param,newpopulation )
%AIM2_COMPUTE_NDCG2 此处显示有关此函数的摘要
%   此处显示详细说明
intanceFront1_Pop = param.intanceFront1_Pop(:,3:end);
X = param.X;
Y = param.Y;
 pairwiseLabel = param.pairwiseLabel;
pairwiseFeature = param.pairwiseFeature;
poplength = param.FSnumber;
N = size(newpopulation,1);
aim2 = zeros(N,1);
intanceFront1_PopLength = size(intanceFront1_Pop,1);
% classifier.W= zeros(intanceFront1_PopLength,poplength);
% classifier = repmat(classifier,N,1);

for i=1:N
    C_NDCG = zeros(intanceFront1_PopLength,1);
    Feature_selected = logical(newpopulation(i,1:end));
    %     symbolFlag = newpopulation(i,end-4);
    %     if symbolFlag==0
    %     parameterC = 2^(- bin2dec(num2str(newpopulation(i,end-3:end))));
    %     else
    %     parameterC = 2^( bin2dec(num2str(newpopulation(i,end-3:end))));
    %     end
    for j=1:intanceFront1_PopLength
        currentIstanceSelected = logical(intanceFront1_Pop(j,:));
        
        Indices = crossvalind('Kfold', size(param.countData,1),5);
        Indices(Indices==5)=0;
        crossoverFoldNDCG = zeros(5,1);
        for I=1:5
            
            CIndices = true(1,size(Indices,1));
            I_1 = rem(I,5);
            I_2 = rem((I+1),5);
            I_3 = rem((I+2),5);
            CIndices(Indices==I_1|Indices==I_2|Indices==I_3)=0;
            
            for K = 1:size(param.countData,1)
                if K==1 &&  CIndices(K)==0
                    currentIstanceSelected(1:param.countData(K,3))=0;
                elseif CIndices(K)==0
                    currentIstanceSelected(param.countData(K-1,3)+1:param.countData(K,3))=0;
                end
            end
            if sum(currentIstanceSelected)==0
                  continue;
            end 
            feature = sparse(pairwiseFeature(currentIstanceSelected,Feature_selected));
            label   =pairwiseLabel(currentIstanceSelected);
          %  model = train(label,feature,['-s 2 -c 1 -q']);
             model = train(label,feature,'-s 2 -c 1 -q');
            W = model.w;
            complete_W = zeros(1,poplength);
          %  fprintf(' %d_%d_%d\n',length(complete_W),sum(Feature_selected),length(W));
            complete_W(Feature_selected) = W;
            
            % [regressWX]= predict(param.Y_,sparseX,model,'-q');
            
            NDCGCount = 0;
            NDCG = 0;
            for K=1:size(param.countData,1)
                if K==1 && CIndices(K)==1
                    dataX = X(1:param.countData(1,2),:);
                    dataY{1} = Y{1};
                    NDCGCount = NDCGCount+1;
                    NDCG = NDCG+compute_ndcg4(complete_W*dataX',dataY,10);
                elseif CIndices(K)==1
                    dataX = X(param.countData(K-1,2)+1:param.countData(K,2),:);
                    dataY{1} = Y{K};
                    NDCGCount = NDCGCount+1;
                    NDCG = NDCG+compute_ndcg4(complete_W*dataX',dataY,10);
                end
                
            end
            crossoverFoldNDCG(I) = NDCG/NDCGCount;
        end
        %     feature = sparse(pairwiseFeature(currentIstanceSelected,Feature_selected));
        %     label   = pairwiseLabel(currentIstanceSelected);
        %     model = train(label,feature,['-s 2 -c ' num2str(param.BestParameterC) ' -q']);
        %     W = model.w;
        %     complete_W = zeros(1,poplength);
        %     complete_W(Feature_selected) = W;
        % %     classifier(i).W(j,:) = complete_W;
        %     NDCG =  compute_ndcg4(complete_W*X',Y,10);
        if sum(crossoverFoldNDCG~=0)==0
         C_NDCG(j) = 0;
         continue;
        end
        C_NDCG(j) = sum(crossoverFoldNDCG)/sum(crossoverFoldNDCG~=0);
    end
    aim2(i) =1- max(C_NDCG);
end
end




