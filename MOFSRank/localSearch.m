function [ FunctionValue, Population] = localSearch( param,FunctionValue,Population )
%LOCALSEARCH 此处显示有关此函数的摘要
%   此处显示详细说
PopNum = size(Population,1);
PopSize = size(Population,2);
AI = zeros(5000,PopSize+2);
AINum = 0;
for i=1: PopNum
    currentPop = Population(i,:);
    currentAim1 = FunctionValue(i,1);
    currentAim2 = FunctionValue(i,2);
     flag = true;
    for k=1:1
        j = randi([1 PopSize]);
        changePop = currentPop;
        changePop(j) = 1-changePop(j);
        changeAim1 = aim1_compute( param,changePop );
        changeAim2 = 1- aim2_compute( param,changePop );
        if currentAim1<changeAim1
           if currentAim2<=changeAim2
               if flag
                   AI(AINum+1,:) = [currentPop,currentAim1,currentAim2];
                   AINum = AINum+1;
                   flag = false;
               end
           else
            AI(AINum+1,:) = [changePop,changeAim1,changeAim2];
            AINum = AINum+1;
             if flag
                   AI(AINum+1,:) = [currentPop,currentAim1,currentAim2];
                   AINum = AINum+1;
                   flag = false;
             end
           end
        else
            if currentAim2>changeAim2
                    AI(AINum+1,:) =  [changePop,changeAim1,changeAim2];
                    AINum = AINum+1;    
            else
             AI(AINum+1,:) = [changePop,changeAim1,changeAim2];
            AINum = AINum+1;
             if flag
                    AI(AINum+1,:) = [currentPop,currentAim1,currentAim2];
                    AINum = AINum+1;
                    flag = false;
             end
            end
        end
        
    end
end
AI = AI(1:AINum,:);
N = PopNum;
FunctionValue = AI(:,PopSize+1:PopSize+2);
[FrontValue,MaxFront] = F_NDSort(FunctionValue,'all');
CrowdDistance         = F_distance(FunctionValue,FrontValue);
%
while(1)
    NoN         = numel(FrontValue,FrontValue<MaxFront);
if (NoN>N)
    MaxFront = MaxFront-1;
else
    break;
end
end

Next(1:NoN) = find(FrontValue<MaxFront);
%
Last          = find(FrontValue==MaxFront);
[~,Rank]      = sort(CrowdDistance(Last),'descend');
Next(NoN+1:N) = Last(Rank(1:N-NoN));
FunctionValue = FunctionValue(Next,:);
Population = AI(Next,1:PopSize);
end

