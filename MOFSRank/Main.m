clear ; 
tic;
%%  ���ݳ�ʼ����
%�������ݼ�
addpath(genpath('./'));
 dataname_ = {'HP04'};%���ݼ�
 for dataNums =1:1%ѡ�����ݼ�

fprintf('%s',dataset);
dataset = dataname_{dataNums};
if ~exist(['Result\' dataset '\outputWX']','file')
    mkdir(['Result\' dataset],['outputWX']);
end
% if ~exist(['Result\' dataset '\MulitiRANK'],'file')
%     mkdir(['Result\' dataset],['MulitiRANK']);
% end
if ~exist(['Result\' dataset '\PopResultRANK'],'file')
    mkdir(['Result\' dataset],['PopResultRANK']);
end
for foldNums =1:1%si��ʾ�ڼ���
    Fold_number =foldNums;
    fprintf(' fold=%d',Fold_number);

    dname = ['H:\ʵ��\mulit\���ҵĴ���\matlab\dataset\' dataset '\Fold' num2str(Fold_number) '\'];
    param.dname = dname;
    [X, Y] = read_letor([dname '/train.txt']);%�������ݼ�
    [Xt,Yt] = read_letor([dname '/vali.txt']);
    [Xtest,Ytest] = read_letor([dname '/test.txt']);
        param.dataset = dataset;
        param.Fold_number=Fold_number;
        param.X = [X;Xt];
        param.Y = [Y,Yt]; 
        param.Xt = Xt;
        param.Yt = Yt;
        param.Xtest = Xtest;
        param.Ytest =  Ytest ;
    [ pairwiseFeature,pairwiseLabel,param] = getPairwise(  param.X, param.Y,param );%����������
    clear X Xt Xtest Y Yt Ytest;%����ռ䣬�ܶ��ʵ��ע�͵�
    for runtimes =3:3%s2��ʾ�ڼ���ʵ��
        param.times = runtimes;
         fprintf(' times=%d',runtimes);
        %���ݼ���ṹ��
        %% ����Ԥ����
        param.ISnumber =size(pairwiseFeature,1);
        param.FSnumber = size(pairwiseFeature,2);
        %
       param.pairwiseFeature = pairwiseFeature;
       param.pairwiseLabel  = pairwiseLabel;
        %
        [Each_Corr] = abs(getCorr( param.X ));%���������
        param.Each_Corr = Each_Corr;
        param.Person_Sum_corr = sum(Each_Corr)-1;
        [Person_Importance_Weight] = Compute_Weight( param );%����Ȩ��
        param.Person_Importance_Weight = Person_Importance_Weight;
        %
     clear pairwiseFeature pairwiseLabel Each_Corr Person_Importance_Weight;%����ռ䣬�ܶ��ʵ��ע�͵�
        %% ���ֲ������ã���������㷨����
        %ʵ��ѡ��
        param.N1=50;%��Ⱥ��С
        param.Generations1=10;%��������
        %�����㷨 %ʵ��ѡ��
% %       
%         param.paramterC=-1;
        [intanceFront1_Pop]=instnceSelected(param);%ʵ��ѡ��
        Q = intanceFront1_Pop(:,1)~=0;
        intanceFront1_Pop  = intanceFront1_Pop(Q,:);
        param.intanceFront1_Pop = intanceFront1_Pop;
        %         %����ѡ��
        param.N2=50;%��Ⱥ��С
        param.Generations2=20;%��������
        
        foldName = ['fold'  num2str(foldNums) 'Generations1_2=' num2str(param.Generations1) '_' num2str(param.Generations2)]; %'_k1=' num2str(param.k1) '_k2=' num2str(param.k2)];
        param.foldName=foldName;
        if ~exist(['Result\' dataset '\BOOSTING\' foldName],'file')
            mkdir(['Result\' dataset  '\BOOSTING\'],foldName);
        end
%         if ~exist(['Result\' dataset '\MulitiRANK\' foldName],'file')
%             mkdir(['Result\' dataset '\MulitiRANK\'],foldName);
%         end
        if ~exist(['Result\' dataset '\outputWX\' foldName],'file')
            mkdir(['Result\' dataset '\outputWX\'],foldName);
        end
        [featureFront1_Pop,featureGene]=  featureSelectedNSGAII(param);%����ѡ��
        M = featureFront1_Pop(:,1)~=0;
        featureFront1_Pop = featureFront1_Pop(M,:);
        param.featureFront1_Pop = featureFront1_Pop;
%         save(['Result\' dataset '\MulitiRANK\' foldName ...
%             '\' foldName '_times=' num2str(runtimes) '.mat'],'intanceFront1_Pop','featureFront1_Pop');
%       %  (param);
      [ensembleBeforeLater]=BOOSTING(param);%����
    end
end
 end

toc;
% end