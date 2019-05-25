function [ aim1_feature_num ] = aim1_compute( param ,newpopulation)
%COMPUTE_FEATURE_NUM 此处显示有关此函数的摘要
%   此处显示详细说明
% A= newpopulation(:,1:param.FSnumber);
aim1_feature_num = sum(newpopulation,2);
end

