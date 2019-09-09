clear;
load('stockbank.mat');
load('GDP.mat');
N=20;

W=zeros(2*N,2*N);

WP=ones(2*N,2*N)-eye(2*N,2*N);
IPR=zeros(2*N,59);
% NPRc=zeros(2*N,60);
Eigenvalue=zeros(2*N,59);
MaxEig=zeros(1,59);
for j=2:60
    WC=zeros(2*N,2*N);
    W=zeros(2*N,2*N);
    WP=ones(2*N,2*N)-eye(2*N,2*N);
    WB=zeros(N,N);
    WS=zeros(N,N);
    stock=zeros(N,N);
    bank=zeros(N,N);
    product=zeros(N,N);
    stock(:,:)=correlation_matrix(j,:,:);
    InvStock=inv(stock);
    InvStock= InvStock-diag(diag(InvStock));
        for kk=1:20
        for jj=1:20
           InvStock(kk,jj)=GDP(kk,1)* InvStock(kk,jj)*GDP(jj,1);
        end;
        end;
    R=reshape(InvStock,1,400);
    InvSt=(InvStock-mean(R))./std(R);
    InvSt=InvSt-diag(diag(InvSt));
    bank(:,:)=wBank(j,:,:);
    RB=reshape(bank,1,400);
    RS=reshape(InvSt,1,400);
    %==================================
    WB=bank;
    [VB, DB] = eig(WB);             %diagonal D is eigenvalues and V is eigenvectors
    EigenvalueB(:,j-1)=diag(DB);
    MaxEig_B(1,j-1)=max(EigenvalueB(:,j-1));
    IPRMAX_B(:,j-1)=1./(sum(VB(:,end).^4));
    IPRB(:,j-1)=1./(sum(VB.^4));
    %==================================
    WS=InvSt;
    [VS, DS] = eig(WS);             %diagonal D is eigenvalues and V is eigenvectors
    EigenvalueS(:,j-1)=diag(DS);
    MaxEig_S(1,j-1)=max(EigenvalueS(:,j-1));
    IPRMAX_S(:,j-1)=1./(sum(VS(:,end).^4));
    IPRS(:,j-1)=1./(sum(VS.^4));
    
    
    %================================zero
    W(1:N,1:N)=bank;
    W(N+1:2*N,N+1:2*N)=-InvSt;
    [V, D] = eig(W);             %diagonal D is eigenvalues and V is eigenvectors
    Eigenvalue(:,j-1)=diag(D);
    MaxEig(1,j-1)=max(Eigenvalue(:,j-1));
    IPRMAX(:,j-1)=1./(sum(V(:,end).^4));
    IPR(:,j-1)=1./(sum(V.^4));
    %=================================one
    WP(1:N,1:N)=bank;
    WP(N+1:2*N,N+1:2*N)=InvSt;
    [Vp, Dp] = eig(WP);             %diagonal D is eigenvalues and V is eigenvectors
    Eigenvaluep(:,j-1)=diag(Dp);
    MaxEig_P(1,j-1)=max(Eigenvaluep(:,j-1));
    IPRMAX_P(:,j-1)=1./(sum(Vp(:,end).^4));
    IPRp(:,j-1)=1./(sum(Vp.^4));
    %==================================coupling
    WC(1:N,1:N)=bank;
    WC(N+1:2*N,N+1:2*N)=InvSt;
    for kk=1:20
        for jj=1:20
            product(kk,jj)=(bank(:,kk)'*InvSt(:,jj))/20;
        end;
    end;
    WC(1:N,N+1:2*N)=product;
    WC(N+1:2*N,1:N)=(product)';
    
    
    
    [Vc, Dc] = eig(WC);             %diagonal D is eigenvalues and V is eigenvectors
    Eigenvaluec(:,j-1)=diag(Dc);
    MaxEig_C(1,j-1)=max(Eigenvaluec(:,j-1));
    IPRMAX_C(:,j-1)=1./(sum(Vc(:,end).^4));
    IPRc(:,j-1)=1./(sum(Vc.^4));
    
    NPRMAX(:,j-1)=(Vc(:,end)).^4;
    NPRMAX_FINAL(:,j-1)=NPRMAX(:,j-1)./mean(NPRMAX(:,j-1));
    SORT_NPR(:,j-1)=sort(NPRMAX(:,j-1));
    for ll=0:5
        MAX_SORT(ll+1,j-1)=find(NPRMAX(:,j-1)==SORT_NPR(40-ll,j-1));
    end;
    
end;



Date=['2000/05';'2000/08';'2000/11';'2001/02';'2001/05';'2001/08';'2001/11';'2002/02';'2002/05';'2002/08';'2002/11'; '2003/02';'2003/05';'2003/08';'2003/11'; '2004/02';'2004/05';'2004/08';'2004/11';'2005/02';'2005/05';'2005/08';'2005/11';'2006/02';'2006/05';'2006/08';'2006/11';'2007/02';'2007/05';'2007/08';'2007/11';'2008/02';'2008/05';'2008/08';'2008/11';'2009/02';'2009/05';'2009/08';'2009/11'; '2010/02';'2010/05';'2010/08';'2010/11';'2011/02';'2011/05';'2011/08';'2011/11';'2012/02';'2012/05';'2012/08';'2012/11';'2013/02';'2013/05';'2013/08';'2013/11';'2014/02';'2014/05';'2014/08';'2014/11'];
celldata = cellstr(Date);

figure
subplot(2,1,1);
plot(MaxEig,'.-');
xlabel('Time');
ylabel('MaxEig');
set(gca, 'XTick',1:10:59, 'xticklabel',celldata(1:10:59));
subplot(2,1,2);
plot(IPRMAX,'.-');
xlabel('Time');
ylabel('IPRMax');
set(gca, 'XTick',1:10:59, 'xticklabel',celldata(1:10:59));
%++++++++++++++++++++++++++++++
figure
subplot(2,1,1);
plot(MaxEig_P,'.-');
xlabel('Time');
ylabel('MaxEig');
set(gca, 'XTick',1:10:59, 'xticklabel',celldata(1:10:59));
subplot(2,1,2);
plot(IPRMAX_P,'.-');
xlabel('Time');
ylabel('IPRMax');
set(gca, 'XTick',1:10:59, 'xticklabel',celldata(1:10:59));
%+++++++++++++++++++
figure
subplot(2,1,1);
plot(MaxEig_C,'.-');
xlabel('Time');
ylabel('MaxEig');
set(gca, 'XTick',1:10:59, 'xticklabel',celldata(1:10:59));
subplot(2,1,2);
plot(IPRMAX_C,'.-');
xlabel('Time');
ylabel('IPRMax');
set(gca, 'XTick',1:10:59, 'xticklabel',celldata(1:10:59));
% subplot(3,1,3);
%
% xlabel('Time');
% ylabel('NPR');
%=================
figure
plot(NPRMAX_FINAL(6,:),'.-');
hold on
plot(NPRMAX_FINAL(16,:),'.-');
hold on
plot(NPRMAX_FINAL(4,:),'.-');
hold on
% plot(NPRMAX_FINAL(16,:),'.-');
plot(NPRMAX_FINAL(21,:),'.-');
% for jj=1:40
% %     S_NPR(:,j)=sort(NPRMAX(:,j));
% %     four_NPR(:,j)=NPRMAX(:,j).^4/
%  MEAN_NPR(jj,:)=four_NPR(jj,:)./mean(NPRMAX(jj,:),2);
% %  plot(NPRc(kkk,:),'.-');
% % hold on;
% end;

% ff=find(MAX_SORT==28);
% havij=['2000/2';'2000/5';'2000/8';'2000/11';'2000/2';'2000/8';'2000/2';'2000/5';'2000/8';'2000/2';'2000/5';'2000/8';'2000/2';'2000/5';'2000/8';'2000/2';'2000/5';'2000/8';'2000/2';'2000/5';'2000/8']

%  havij=['2000/2';'2000/5';'2000/8';'2000/11';'2001/2';'2001/5';'2001/8';'2001/11';'2002/8';'2000/2';'2000/5';'2000/8']


%     ]startDate = datenum('01-01-2000');
% endDate = datenum('01-01-2015');
% xData = linspace(startDate,endDate,60);
% % plot(xData, NPRMAX_FINAL,'.-');
% datetick('x','yyyy','keeplimits');
% xlabel('Time');
% ylabel('NPR');