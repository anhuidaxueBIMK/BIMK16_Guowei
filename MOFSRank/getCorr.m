function [Each_Corr] = getCorr( X )
%GETCORR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

        [Each_Corr] = corr(X);
        Each_Corr(isnan(Each_Corr))=1;
        
     
end

