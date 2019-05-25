function [ensembleBeforeLater] = BOOSTING4_2018_4_24( param )
%BOOSTIN 此处显示有关此函数的摘要
%   此处显示详



intanceFront1_Pop = param.intanceFront1_Pop;
intanceFront1_PopNum = size(intanceFront1_Pop,1);
featureFront1_Pop = param.featureFront1_Pop;

pairwiseFeature = param.pairwiseFeature;
pairwiseLabel = param.pairwiseLabel;

FunctionValuef = featureFront1_Pop(:,1:2);
FunctionValuef(:,2) = 1-FunctionValuef(:,2);
[FrontValuef,~] = F_NDSort(FunctionValuef,'all');
featureFront1_Pop = featureFront1_Pop(FrontValuef==1,:);
% 
featureFront1_PopNum = size(featureFront1_Pop,1);
% 
% 
%
PopW  = zeros(featureFront1_PopNum,param.FSnumber);
% 
X = param.X;
Y = param.Y;
% 
ensembleBeforeC = zeros(1,featureFront1_PopNum);

for i=1:featureFront1_PopNum
    currentFeaturePop = logical(featureFront1_Pop(i,3:end-5));
  %
    PopW2 = zeros(intanceFront1_PopNum,param.FSnumber);
    C_NDCG = zeros(intanceFront1_PopNum,1);
    %
    for j=1:intanceFront1_PopNum
        currentIstancePop = logical(intanceFront1_Pop(j,3:end));
        TrainX = sparse(pairwiseFeature(currentIstancePop,currentFeaturePop));
        TrainY = pairwiseLabel(currentIstancePop);
%         
             model = train(TrainY,TrainX,['-s 2  -q']);
        W = model.w;
        complete_W = zeros(1,param.FSnumber);
        complete_W(currentFeaturePop==1) = W;
        NDCG =  compute_ndcgBOOST4(complete_W*X',Y,10);
        %    
        C_NDCG(j) = NDCG;
        PopW2(j,:) = complete_W;
    end
    [maxvalue,maxIndex]=max(C_NDCG);
    minvalue = min(C_NDCG);
    if maxvalue==minvalue
        maxIndex = intanceFront1_PopNum;
    end
  %%%%%%%%%%%%%%
%    
        [complete_W,NDCG,BestParamterC] =SelectParamterC(pairwiseFeature(logical(intanceFront1_Pop(maxIndex,3:end)),currentFeaturePop),pairwiseLabel(logical(intanceFront1_Pop(maxIndex,3:end))),currentFeaturePop,param,1);
      
        ensembleBeforeC(i)=BestParamterC;
%  
      PopW(i,:) =complete_W;
        featureFront1_Pop(i,2) = NDCG;
        
%     writeIntoText2(param.Xtest* PopW(i,:)','MulitiRANK',...
%         featureFront1_Pop(i,1),featureFront1_Pop(i,2),'test',param);%跑多次
    
end
ensembleBefore =featureFront1_Pop(:,2);
ensembleBeforeLater{1} = ensembleBefore';
ensembleBeforeLater{2} =  ensembleBeforeC;
% plot(featureFront1_Pop(:,1),featureFront1_Pop(:,2),'+b');                 %作图
% hold on;
% axis([0,65,0,1.5]);xlabel('F1_');ylabel('F_2');title('RANK')
%%
% {
selectedFeature = sum(PopW);
selectedFeatureIndex =(selectedFeature~=0);
%
N = 50;
Generations = 20; 
%
poplength1 = sum(selectedFeatureIndex);
%    
poplength2 = poplength1;
%
minvalue1 =  ones(1,poplength1)*(0);%%边界
maxvalue1  =  ones(1,poplength1);
Boundary1    = [maxvalue1;minvalue1];  %上下边界
minvalue2 =  ones(1,poplength2)*(0);%%边界
maxvalue2  =  ones(1,poplength2);
Boundary2    = [maxvalue2;minvalue2];  %上下边界
Coding1      =  'Binary';
Coding2      =  'Real';

%
Population1 = randi([0 1],N,poplength1);%随机初始化
Population2 = rand(N,poplength2);

[BOOST4aim1]=computeBOOST4Aim1(Population1);
[BOOST4aim2,AllW1] = computeBOOST4Aim2(Population1,Population2,PopW,selectedFeatureIndex,param);

%    

FunctionValue  =  [BOOST4aim1,BOOST4aim2];
FrontValue                   = F_NDSort(FunctionValue,'half');
CrowdDistance                = F_distance(FunctionValue,FrontValue);
for Gene = 1 : Generations
    MatingPool            = F_mating(Population1,FrontValue,CrowdDistance);
    MatingPoo2            = F_mating(Population2,FrontValue,CrowdDistance);
    Offspring1             = P_generatorBOOSTING(MatingPool,Boundary1,Coding1,N);
    Offspring2            = P_generatorBOOSTING(MatingPoo2,Boundary2,Coding2,N);
    Population1            = [Population1;Offspring1];
    Population2            = [Population2;Offspring2];
    %       
    [BOOST4aim1]=computeBOOST4Aim1(Offspring1);
    [BOOST4aim2,AllW2]  = computeBOOST4Aim2(Offspring1,Offspring2,PopW,selectedFeatureIndex,param);
    
    if Gene==1
        FunctionValue(N+1:2*N,1:2) = [BOOST4aim1,BOOST4aim2];
        AllW = [AllW1;AllW2];
    else
        FunctionValue(1:N,1:2) = FunctionValue(Next,:);
        FunctionValue(N+1:2*N,1:2) = [BOOST4aim1,BOOST4aim2];
        AllW(1:N,:) = AllW(Next,:);
        AllW(N+1:2*N,:) = AllW2;
    end
    [FrontValue,MaxFront] = F_NDSort(FunctionValue,'half');
    CrowdDistance         = F_distance(FunctionValue,FrontValue);
  
    %
    Next        = zeros(1,N);
    NoN         = numel(FrontValue,FrontValue<MaxFront);
    Next(1:NoN) = find(FrontValue<MaxFront);
    
    %
    Last          = find(FrontValue==MaxFront);
    [~,Rank]      = sort(CrowdDistance(Last),'descend');
    Next(NoN+1:N) = Last(Rank(1:N-NoN));
    
    %
    Population1    = Population1(Next,:);
    Population2    = Population2(Next,:);
    FrontValue    = FrontValue(Next);
   
    CrowdDistance = CrowdDistance(Next);
    %       
end
AllW = AllW(Next,:);
AllW = AllW(FrontValue==1,:);
afterEnsembleSelected = AllW;
afterEnsembleSelected(afterEnsembleSelected~=0)=1;
afterEnsembleSelectedAim = zeros(size(afterEnsembleSelected,1),3);
for i=1:size(afterEnsembleSelected,1)
    currentFeaturePop = logical(afterEnsembleSelected(i,:));
    currentFeatureIS = logical(intanceFront1_Pop(end,3:end));
%    
        [complete_W,NDCG,BestParamterC] =SelectParamterC(pairwiseFeature(currentFeatureIS,currentFeaturePop)...
            ,pairwiseLabel(currentFeatureIS),currentFeaturePop,param,2);
      afterEnsembleSelectedAim(i,3) = BestParamterC;
%   
    afterEnsembleSelectedAim(i,1) = sum(currentFeaturePop);
    afterEnsembleSelectedAim(i,2) = NDCG;
    writeIntoText(param.Xtest*complete_W','BOOSTING',afterEnsembleSelectedAim(i,1) ,  afterEnsembleSelectedAim(i,2) ,'test',param);%跑多次
end
afterEnsembleSelectedAim = sortrows(afterEnsembleSelectedAim);
ensembleLater =  afterEnsembleSelectedAim(:,2);
ensembleBeforeLater{3}=ensembleLater';
ensembleBeforeLater{4}=afterEnsembleSelectedAim(:,3)'; 

plot(afterEnsembleSelectedAim(:,1),afterEnsembleSelectedAim(:,2),'xr');                 %作图
hold on;
axis([0,65,0,1.5]);xlabel('F_1');ylabel('F_2');title('RANK')
end

function [BOOST4aim1]=computeBOOST4Aim1(Population)
BOOST4aim1 = sum(Population,2);
end

function [BOOST4aim2,AllW] = computeBOOST4Aim2(Population1,Population2,PopW,selectedFeatureIndex,param)
X = param.X;
Y = param.Y;
PopNum = size(Population1,1);
BOOST4aim2 = zeros(PopNum,1);
AllW = zeros(PopNum,param.FSnumber);
for i=1:PopNum
    currentPop1 =logical(Population1(i,:));
    %
    completeW =zeros(1,param.FSnumber);
    currentPop2=Population2(i,:);
    currentPop2(~currentPop1)=0;
    completeW(selectedFeatureIndex)=currentPop2;
    AllW(i,:) = completeW;
    %***
    NDCG =  compute_ndcgBOOST4(completeW*X',Y,10);
    BOOST4aim2(i) =1- NDCG;
end
%}
end
%%

function [bestW,ndcg_2,BestParamterC] =SelectParamterC(TrainX,TrainY,currentFeaturePop,param,flag)
if flag==1
    X = param.X;
    Y = param.Y;
elseif flag==2
    X = param.Xtest;
    Y = param.Ytest;
end
TrainX = sparse(TrainX);
ndcg_ =0;
ndcg_2 =0;
bestW = zeros(1,param.FSnumber);
BestParamterC=0;
for i=-10:1:10
    paramterC = 2^i;
    model = train(TrainY,TrainX,['-s 2 -c ' num2str(paramterC) ' -q']);
    W = model.w;
    complete_W = zeros(1,param.FSnumber);
    complete_W(currentFeaturePop==1) = W;
    NDCG =  compute_ndcgBOOST4(complete_W*X',Y,10);
    %   
    NDCG2=NDCG;
    %     
    if NDCG>ndcg_
        ndcg_ = NDCG;
        bestW = complete_W;
        ndcg_2 = NDCG2;
        BestParamterC=paramterC;
    end
end

end

%%
function writeIntoText(output,foldername,featureNumber,ndcgValue,DataType,param)
%   output = output + 1e-10*randn(length(output),1);  % Break ties at random
fname1 = ['Result/' param.dataset '/outputWX/' param.foldName '/' DataType '.fold' num2str(param.Fold_number)...
    '_times=' num2str(param.times) '_FN=' num2str(featureNumber) '_ndcg=' num2str(ndcgValue)];
fname2 = ['Result/' param.dataset '/' foldername '/' param.foldName '/' DataType '.fold' num2str(param.Fold_number)...
    '_times=' num2str(param.times) '_FN=' num2str(featureNumber) '_ndcg=' num2str(ndcgValue)];
save( fname1,'output','-ascii');
% Either copy the evaluation script in the current directory or
% change the line below with the correct path
if   ~strcmpi(param.dataset,'MQ08')&&~strcmpi(param.dataset,'MQ07')
    %   %跑3.0数据集
    system(['perl Eval-Score-3.0.pl ' param.dname  DataType ...
        '.txt ' fname1 ' ' fname2 '.metric 0']);
else
    %跑4.0数据集
    system(['perl Eval-Score-4.0.pl ' param.dname  DataType ...
        '.txt ' fname1 ' ' fname2 '.metric 0']);
end
end
%%
function writeIntoText2(output,foldername,featureNumber,ndcgValue,DataType,param)
%   output = output + 1e-10*randn(length(output),1);  % Break ties at random
fname1 = ['Result/' param.dataset '/outputWX/' param.foldName '/' DataType '.fold' num2str(param.Fold_number)...
    '_times=' num2str(param.times) '_FN=' num2str(featureNumber)  '_ndcg=' num2str(ndcgValue)];
fname2 = ['Result/' param.dataset '/' foldername '/' param.foldName '/' DataType '.fold' num2str(param.Fold_number)...
    '_times=' num2str(param.times) '_FN=' num2str(featureNumber)  '_ndcg=' num2str(ndcgValue)];
save( fname1,'output','-ascii');
% Either copy the evaluation script in the current directory or
% change the line below with the correct path
%   %跑3.0数据集
if ~strcmpi(param.dataset,'MQ08')&&~strcmpi(param.dataset,'MQ07')
    system(['perl Eval-Score-3.0.pl ' param.dname  DataType ...
        '.txt ' fname1 ' ' fname2 '.metric 0']);
else
    %跑4.0数据集
    system(['perl Eval-Score-4.0.pl ' param.dname  DataType ...
        '.txt ' fname1 ' ' fname2 '.metric 0']);
end
end
%%
function meanndcg = compute_ndcgBOOST4(Y,Yt,p)
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
end