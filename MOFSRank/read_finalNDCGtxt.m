function [ NDCG ] = read_finalNDCGtxt( filepath )
%READ_FINALNDCGTXT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

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

