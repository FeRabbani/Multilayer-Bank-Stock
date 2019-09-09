% load('International_exchange_29.xlsx');

% j=0;
% for i=95:157
% j=j+1;
% BankData(:,j)=Internationalexchange29(:,i);
% end;
% C = table2cell(BankData);
% BankMatrix = cell2mat(C);
% FinalBank=zeros(841,63);
% l=551;
% for k=639:667 %Netherland
%     l=l+1;
%    FinalBank(l,:) =BankMatrix(k,:);
% end;
%
% for k=30:58  %Netherland
%
%    FinalBank(k,:) =BankMatrix(k,:);
% end;
% Total=zeros(400,63);


%normalize bank data
% for i=1:841
%     
%     F(i,:)= (FinalBank(i,:)-mean(FinalBank(i,:)))./var(FinalBank(i,:));
% end;
% F(isnan(F))=0;
% FinalBank=F;
s=0;
for m=1:29:560         % select countaries we want
                        % in this step we remove 9 countries from matrix
    h=1+s;
    Total(h:h+4,:)=FinalBank(m:m+4,:);
    Total(h+5,:)=FinalBank(m+7,:);
    Total(h+6:h+7,:)=FinalBank(m+9:m+10,:);
    Total(h+8:h+9,:)=FinalBank(m+13:m+14,:);
    Total(h+10,:)=FinalBank(m+16,:);
    Total(h+11:h+19,:)=FinalBank(m+20:m+28,:);
    s=h+19;
end;


f=0;
for m=1:20:400    % in this step we move rows according to name of cuntries in stock matrix
    a=1+f;
    TotalData(a,:)= Total(m+12,:);
    TotalData(a+1,:)= Total(m+1,:);
    TotalData(a+2,:)= Total(m+6,:);
    TotalData(a+3,:)= Total(m+7,:);
    TotalData(a+4,:)= Total(m+2,:);
    TotalData(a+5,:)= Total(m+18,:);
    TotalData(a+6,:)= Total(m+8,:);
    TotalData(a+7,:)= Total(m+3,:);
    TotalData(a+8,:)= Total(m+11,:);
    TotalData(a+9,:)= Total(m+9,:);
    TotalData(a+10,:)= Total(m+14,:);
    TotalData(a+11,:)= Total(m+10,:);
    TotalData(a+12,:)= Total(m+5,:);
    TotalData(a+13,:)= Total(m+16,:);
    TotalData(a+14,:)= Total(m+17,:);
    TotalData(a+15,:)= Total(m+19,:);
    TotalData(a+16,:)= Total(m,:);
    TotalData(a+17,:)= Total(m+4,:);
    TotalData(a+18,:)= Total(m+15,:);
    TotalData(a+19,:)= Total(m+13,:);
    
    f=a+19;
end;






