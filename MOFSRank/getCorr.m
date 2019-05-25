function [Each_Corr] = getCorr( X )
%GETCORR 此处显示有关此函数的摘要
%   此处显示详细说明

        [Each_Corr] = corr(X);
        Each_Corr(isnan(Each_Corr))=1;
        
     
end

