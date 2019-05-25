function [ pairwiseFeature,pairwiseLabel,param] = getPairwise( X,Y,param )
%GETPAIRWISE 此处显示有关此函数的摘要
%   此处显示详细说明
Ylength = length(Y);
XDim = size(X,2);
pairwiseFeature = zeros(1000000,XDim);
pairwiseLabel = zeros(1000000,1);
pairwiseNum = 0;
pairwiseFeatureIndix = zeros(1000000,1);
FeatureIndix = zeros(size(X,1),1);
cunrrentXSite = 0;
countData = zeros(Ylength,3);
for i=1:Ylength
    currentQidY = Y{i};
    currentQidYLength = length(currentQidY);
    currentQidX = X(cunrrentXSite+1:cunrrentXSite+currentQidYLength,:);
    FeatureIndix(cunrrentXSite+1:cunrrentXSite+currentQidYLength) =i;
    for j=1:currentQidYLength
%         if (currentQidY(j)==0),continue;end 
        for k=j+1:currentQidYLength
            if currentQidY(j)~=currentQidY(k)
                pairwiseNum = pairwiseNum+1;
                flag = currentQidY(j) - currentQidY(k);
                if flag>0
                    pairwiseLabel(pairwiseNum) =1;
                    if flag==2
                    pairwiseFeature(pairwiseNum,:) = 3.*(currentQidX(j,:)-currentQidX(k,:));
                    else
                    pairwiseFeature(pairwiseNum,:) = currentQidX(j,:)-currentQidX(k,:);
                    end
                else
                     pairwiseLabel(pairwiseNum) =-1;
                    if flag==-2
                    pairwiseFeature(pairwiseNum,:) = 3.*(currentQidX(j,:)-currentQidX(k,:));
                    else
                    pairwiseFeature(pairwiseNum,:) = currentQidX(j,:)-currentQidX(k,:);
                    end
                end
                pairwiseFeatureIndix(pairwiseNum) = i;
            end
        end
    end
    cunrrentXSite = cunrrentXSite+currentQidYLength;
    
    if i==1
    countData(i,1) = i;
    countData(i,2) =   currentQidYLength;
    countData(i,3) =  pairwiseNum;
    else
    countData(i,1) = i;
    countData(i,2) =   countData(i-1,2) + currentQidYLength;
    countData(i,3) =  pairwiseNum;
    end
    
end
pairwiseFeature = pairwiseFeature(1:pairwiseNum,:);
pairwiseLabel = pairwiseLabel(1:pairwiseNum,:);
param.countData = countData;




end

