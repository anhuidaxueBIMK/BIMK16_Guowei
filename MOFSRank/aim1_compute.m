function [ aim1_feature_num ] = aim1_compute( param ,newpopulation)
%COMPUTE_FEATURE_NUM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% A= newpopulation(:,1:param.FSnumber);
aim1_feature_num = sum(newpopulation,2);
end

