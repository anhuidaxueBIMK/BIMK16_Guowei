function a=removeDuplicate(a)
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