function [Front1_Pop] = instnceSelected( param )

%-----------------------------------------------------------------------��������
N=param.N1;%��Ⱥ��С
Generations=param.Generations1;%��������
ISnumber  =  param.ISnumber;%��Ӧ����Ⱥ��ά�ȣ�������

%  Population = lhsamp(N, ISnumber);%���������������
% Population = [zeros(N/2, ISnumber);ones(N/2, ISnumber)];%��ʼ����Ⱥ
Population = zeros(N,ISnumber);
% end
minvalue1  =  ones(1,ISnumber)*(0);%%�߽�
maxvalue1  =  ones(1,ISnumber);
Boundary    = [maxvalue1;minvalue1];                          %%���±߽�
Coding      =  'Binary';

%��ʼ�������������ݿ�����
% 
[ istanceAim1 ] = istanceAim1Compute(param,Population );
[ istanceAim2 ] = istanceAim2Compute(param,Population );

FunctionValue  =  [istanceAim1,istanceAim2];
FrontValue                   = F_NDSort(FunctionValue,'half');
CrowdDistance                = F_distance(FunctionValue,FrontValue);
Next = 1:N;
%% -----------------------------------------------------------------ѭ������
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
%+----------�ҳ���һǰҳ��⣬��������
LAST_FunctionValue = FunctionValue(Next,:);
Front1_Pop(:,1:2) =  LAST_FunctionValue(FrontValue==1,:);
Front1_Pop(:,2) = 1 - Front1_Pop(:,2);
Front1_Pop(:,3:ISnumber+2) = Population(FrontValue==1,1:end);
Front1_Pop = sortrows(Front1_Pop);
end
 %% �����������
function S = lhsamp(m, n)%�����������
if nargin < 1, m = 1; end
if nargin < 2, n = m; end

S = zeros(m,n);
for i = 1 : n
    S(:, i) = (rand(1, m) + (randperm(m) - 1))' / m;
end
S(S>0.5) = 1;
S(S<=0.5) = 0;
end