function [Front1_Pop,Gene1] = featureSelectedNSGAII( param )

%-----------------------------------------------------------------------参数设置

N=param.N2;
Generations=param.Generations2;
poplength  = param.FSnumber;%对应于种群的维度，即列数
%
Population = lhsamp(N, poplength);%拉丁超立方体抽样
% Population = zeros(N,poplength);%全0初始
% Population = ones(N,poplength);%全1初始
minvalue1  =  ones(1,poplength)*(0);%%边界
maxvalue1  =  ones(1,poplength);
Boundary    = [maxvalue1;minvalue1];                          %%上下边界
Coding      =  'Binary';
%                        %%上下边界
[ aim1 ] = aim1_compute(param,Population );
[ aim2 ] =aim2_compute(param,Population );
%
FunctionValue  =  [aim1,aim2];
FrontValue                   = F_NDSort(FunctionValue,'half');
CrowdDistance                = F_distance(FunctionValue,FrontValue);
Next = 1:N;
flagNum = 0;
%% -----------------------------------------------------------------循环迭代
for Gene1 = 1 : Generations 
    oldFunctionValue = FunctionValue(Next,:);
MatingPool            = F_mating(Population(:,1:end),FrontValue,CrowdDistance);
%
k1 =0;
k2 =0;
if Gene1==1
    k1 =100;
    k2 = 100;
elseif Gene1<10
    k1=0.1+k1/exp(Gene1);
    k2=0.1+k2/exp(Gene1);
end
Offspring             = P_generatorCW(MatingPool,Boundary,Coding,N,k1,k2,param);
Population            = [Population;Offspring];
[ aim1 ] = aim1_compute( param,Offspring );
[ aim2 ] =aim2_compute(param,Offspring );
FunctionValue(1:N,:) = FunctionValue(Next,:);
FunctionValue(N+1:2*N,:) = [aim1,aim2];
%+++++++++++++
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
Population    = Population(Next,:);
FrontValue    = FrontValue(Next);
%
CrowdDistance = CrowdDistance(Next);
 %% 判断是否收敛
u = 0;
currentFunctionValue =  FunctionValue(Next,:);
    for NP = 1:N
        for m=1:2
            u = u + abs(oldFunctionValue(NP,m) - currentFunctionValue(NP,m) );
        end
    end
    %
    if u<=1
     flagNum = flagNum+1;
    else
       flagNum=0; 
    end
    %
    if flagNum>=4
        break;
    end
end
LAST_FunctionValue = FunctionValue(Next,:);
Front1_Pop(:,1:2)  = LAST_FunctionValue;
 Front1_Pop(:,2) = 1- Front1_Pop(:,2);
 Front1_Pop(:,3:poplength+2) = Population;
end
%% 超立方体抽样
function S = lhsamp(m, n)%超立方体抽样
if nargin < 1, m = 1; end
if nargin < 2, n = m; end

S = zeros(m,n);
for i = 1 : n
    S(:, i) = (rand(1, m) + (randperm(m) - 1))' / m;
end
S(S>0.5) = 1;
S(S<=0.5) = 0;
end

