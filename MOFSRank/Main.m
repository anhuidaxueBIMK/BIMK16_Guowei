clear ; 
tic;
%%  数据初始部分
%加载数据集
addpath(genpath('./'));
 dataname_ = {'HP04'};%数据集
 for dataNums =1:1%选择数据集

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
for foldNums =1:1%si表示第几折
    Fold_number =foldNums;
    fprintf(' fold=%d',Fold_number);

    dname = ['H:\实验\mulit\寝室的代码\matlab\dataset\' dataset '\Fold' num2str(Fold_number) '\'];
    param.dname = dname;
    [X, Y] = read_letor([dname '/train.txt']);%加载数据集
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
    [ pairwiseFeature,pairwiseLabel,param] = getPairwise(  param.X, param.Y,param );%产生样本对
    clear X Xt Xtest Y Yt Ytest;%清理空间，跑多次实验注释掉
    for runtimes =3:3%s2表示第几次实验
        param.times = runtimes;
         fprintf(' times=%d',runtimes);
        %数据加入结构体
        %% 数据预处理
        param.ISnumber =size(pairwiseFeature,1);
        param.FSnumber = size(pairwiseFeature,2);
        %
       param.pairwiseFeature = pairwiseFeature;
       param.pairwiseLabel  = pairwiseLabel;
        %
        [Each_Corr] = abs(getCorr( param.X ));%计算相关性
        param.Each_Corr = Each_Corr;
        param.Person_Sum_corr = sum(Each_Corr)-1;
        [Person_Importance_Weight] = Compute_Weight( param );%计算权重
        param.Person_Importance_Weight = Person_Importance_Weight;
        %
     clear pairwiseFeature pairwiseLabel Each_Corr Person_Importance_Weight;%清理空间，跑多次实验注释掉
        %% 部分参数设置，进入进化算法部分
        %实例选择
        param.N1=50;%种群大小
        param.Generations1=10;%迭代次数
        %进化算法 %实例选择
% %       
%         param.paramterC=-1;
        [intanceFront1_Pop]=instnceSelected(param);%实例选择
        Q = intanceFront1_Pop(:,1)~=0;
        intanceFront1_Pop  = intanceFront1_Pop(Q,:);
        param.intanceFront1_Pop = intanceFront1_Pop;
        %         %特征选择
        param.N2=50;%种群大小
        param.Generations2=20;%迭代次数
        
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
        [featureFront1_Pop,featureGene]=  featureSelectedNSGAII(param);%特征选择
        M = featureFront1_Pop(:,1)~=0;
        featureFront1_Pop = featureFront1_Pop(M,:);
        param.featureFront1_Pop = featureFront1_Pop;
%         save(['Result\' dataset '\MulitiRANK\' foldName ...
%             '\' foldName '_times=' num2str(runtimes) '.mat'],'intanceFront1_Pop','featureFront1_Pop');
%       %  (param);
      [ensembleBeforeLater]=BOOSTING(param);%集成
    end
end
 end

toc;
% end