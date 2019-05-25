function [Front1_Pop] = instnceSelected( param )

%-----------------------------------------------------------------------参数设置
N=param.N1;%种群大小
Generations=param.Generations1;%迭代次数
ISnumber  =  param.ISnumber;%对应于种群的维度，即列数

%  Population = lhsamp(N, ISnumber);%拉丁超立方体抽样
% Population = [zeros(N/2, ISnumber);ones(N/2, ISnumber)];%初始化种群
Population = zeros(N,ISnumber);
% end
minvalue1  =  ones(1,ISnumber)*(0);%%边界
maxvalue1  =  ones(1,ISnumber);
Boundary    = [maxvalue1;minvalue1];                          %%上下边界
Coding      =  'Binary';

%初始迭代及代理数据库设置
% 
[ istanceAim1 ] = istanceAim1Compute(param,Population );
[ istanceAim2 ] = istanceAim2Compute(param,Population );

FunctionValue  =  [istanceAim1,istanceAim2];
FrontValue                   = F_NDSort(FunctionValue,'half');
CrowdDistance                = F_distance(FunctionValue,FrontValue);
Next = 1:N;
%% -----------------------------------------------------------------循环迭代
for Gene1 = 1 : Generations
MatingPool            = F_mating(Population(:,1:end),FrontValue,CrowdDistance);
Offspring           = P_generator(MatingPool,Boundary,Coding,N);
Population            = [Population;Offspring];
[ istanceAim1 ] = istanceAim1Compute( param,Offspring );
[ istanceAim2] =istanceAim2Compute( param,Offspring );
%
FunctionValue(1:N,:) = FunctionValue(Next,:);
FunctionValue(N+1:2*N,:) = [istanceAim1,istanceAim2];
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
end
%+----------找出第一前页面解，放入数组
LAST_FunctionValue = FunctionValue(Next,:);
Front1_Pop(:,1:2) =  LAST_FunctionValue(FrontValue==1,:);
Front1_Pop(:,2) = 1 - Front1_Pop(:,2);
Front1_Pop(:,3:ISnumber+2) = Population(FrontValue==1,1:end);
Front1_Pop = sortrows(Front1_Pop);
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