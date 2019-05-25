function [BestParameterC] =inStartSelectedParamterC(TrainX,TrainY,param)

X = param.Xt;
Y = param.Yt; 
% X2 = param.Xtest;
% Y2 = param.Ytest; 
TrainX = sparse(TrainX);
ndcg_ =0;
% ndcg_2 =0;
% bestW = zeros(1,param.FSnumber);
for i=-15:1:20
    paramterC = 2^i;
    model = train(TrainY,TrainX,['-s 2 -c ' num2str(paramterC) ' -q']);
    W = model.w;
    NDCG =  compute_ndcg4(W*X',Y,10); 
%     NDCG2 =  compute_ndcgBOOST4(complete_W*X2',Y2,10); 
    if NDCG>ndcg_
        ndcg_ = NDCG;
        BestParameterC = paramterC;
%         bestW = complete_W;
%         ndcg_2 = NDCG2;
    end
end

end