% clc
% clear  all
N=20;
N_MC=800;
p=0.1;q=0.00000099;

w=zeros(2*N,2*N);
wp=zeros(2*N,2*N);

 load('stockbank.mat');
 
    stock(:,:)=correlation_matrix(28,:,:);
    InvStock=inv(stock);
    InvStock= InvStock-diag(diag(InvStock));
    R=reshape(InvStock,1,400);
    InvSt=(InvStock-mean(R))./std(R);
    InvSt=InvSt-diag(diag(InvSt));
    bank(:,:)=wBank(28,:,:);
    RB=reshape(bank,1,400);
    RS=reshape(InvSt,1,400);

% [a]=generate_random_graph(2,N,3)    ;
% [b]=generate_random_graph(2,N,3)    ;
% % spy(field1)
[a]=bank;
[b]=InvSt;
% 
% [a]=generate_random_graph(2,N,3)    ;
% [b]=generate_random_graph(2,N,3)    ;
% spy(field1)

hold on
w(1:N,1:N)=a;
w(N+1:2*N,N+1:2*N)=b;
[V, D] = eig(w);             %diagonal D is eigenvalues and V is eigenvectors
Eigenvalue=diag(D);
IPR=zeros(size(D,1),1);
IPR=1./(sum(V.^4));
%  figure
% plot(Eigenvalue,IPR,'o');
% xlabel('Eigenvalue');
% ylabel('IPR');
%==========================================
for i=1:100
    for rep=1:30
        A=zeros(N,N);
        A= triu(randi([1 100], N,N));
        A = A+A';
        zero=find(A>i);
        one=find(A<=i);
        A(zero)=0;
        A(one)=1;
        wp(1:N,1:N)=a;
        wp(N+1:2*N,N+1:2*N)=b;
        wp(N+1:2*N,1:N)=A;
        wp(1:N,N+1:2*N)=A;
        
        [Vp, Dp] = eig(wp);             %diagonal D is eigenvalues and V is eigenvectors
        Eigenvaluep=diag(Dp);
        Max_EIG(rep,i)=max(Eigenvaluep);
        IPRp=zeros(size(Dp,1),1);
        IPRp=1./(sum(Vp.^4));
        Max_IPR(rep,i)=max(IPRp);
    end;
    MAX_EIG(i,1)=mean(Max_EIG(:,i));
    MAX_IPR(i,1)=mean(Max_IPR(:,i));
end;
% figure
% plot(Eigenvalue1,IPR1,'o');
% xlabel('Eigenvalue');
% ylabel('IPR');
figure;
plot( MAX_EIG,'o-');
ylabel('MAX-Eig');
xlabel('Percent of coupling');
figure
plot( MAX_IPR,'o-');
ylabel('MAX-IPR');
xlabel('Percent of coupling');

%================
% [N,edges] = histcounts(Eigenvalue1,50);
% plot(edges(1:50),N,'o-r')