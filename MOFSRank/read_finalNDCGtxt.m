function [ NDCG ] = read_finalNDCGtxt( filepath )
%READ_FINALNDCGTXT 此处显示有关此函数的摘要
%   此处显示详细说明

filename = [filepath '.metric']; 
f = fopen(filename);
l = fgetl(f);
l = fgetl(f);
l = fgetl(f);
l = fgetl(f);
l = fgetl(f);
[buf, ~, ~, ind] = sscanf(l,'%s NDCG:'); 
     l(1:ind-1)=[];
      [nqid, foo1, foo2, ind] = sscanf(l,'%f'); 
      NDCG = nqid(end);
fclose(f);
end

