function [ surrogateTrainSet ] = slecteDB_Point( Offspring,DB,surrogateTrainSet )
%SLECTEDB_POINT 此处显示有关此函数的摘要
%   此处显示详细说明
addSize = 50;
epilon = 0.3;
effectiveData.X = DB.X(1:DB.vernier,:);
effectiveData.Y = DB.Y(1:DB.vernier,:);
centerPoint = mean(Offspring);
maxDistance=max(pdist2(centerPoint,Offspring));
DBdistance = pdist2(centerPoint,effectiveData.X);
selectPointIndex = find(DBdistance<=(maxDistance+epilon));
if length(selectPointIndex)>addSize

    selectePartPointIndex = randperm(length(selectPointIndex));
%     selectePartPointIndex = randi([1 sum(selectPointIndex)],1,addSize);
    surrogateTrainSet.TrainX(end+1:end+addSize,:) = DB.X(selectPointIndex(selectePartPointIndex(1:addSize)),:);
    surrogateTrainSet.TrainY(end+1:end+addSize) = DB.Y(selectPointIndex(selectePartPointIndex(1:addSize)),2);
else
    surrogateTrainSet.TrainX(end+1:end+length(selectPointIndex),:) = DB.X(selectPointIndex,:);
    surrogateTrainSet.TrainY(end+1:end+length(selectPointIndex)) = DB.Y(selectPointIndex,2);
end

end

