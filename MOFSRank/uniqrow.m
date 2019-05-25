function [A,flag]=uniqrow(A)
%出进来的A是按照第M+1行排序的
[m,n] = size(A);
C = zeros(m,n);
L  = unique(A(:,n-1));
flag = 0;
for i=1:length(L)
    sitepre = find(L(i)<=A(:,n-1));
    sitepre = sitepre(1);
    siterear = find(L(i)+1<=A(:,n-1));
   
%    site = site(1);
if isempty(siterear)
     B = A(sitepre:end,:);
else
     siterear = siterear(1);
     B = A(sitepre:siterear-1,:);
end
 
    B=myfun( B);
     [m2]= size(B,1);
     C(flag+1:flag+m2,:) = B;
     flag = flag+m2;
end
A = C(1:flag,:);
end
%%
function a=myfun(a)
[m,n]=size(a);
b=zeros(m,1);
for i=1:m-1
for j=i+1:m
if sum(a(i,:)==a(j,:))==n
b(j,1)=1;
end
 end
end
for i=m:-1:1
if b(i,1)==1
a(i,:)=[];
 end
end
end