% clear
% load('finalbankseason');

wBank=zeros(63,20,20);

for l=1:63
a=zeros(20,20);
k=0;   
for i=1:20
    j=k+1;
    a(i,:)=TotalData(j:j+19,l)';
k=j;
end;
wBank(l,:,:)=a;
end;