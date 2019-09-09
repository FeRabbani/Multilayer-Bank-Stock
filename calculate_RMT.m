% clear;
%data=dlmread('./market_data_preformat');
%data=dlmread('./sse_preformat_fix');
%data=dlmread('./sp_preformat_fix_final');
%data=dlmread('./tehran_preformat_new_fix2');
%data=dlmread('./DAX_preformat_fix');
%data=dlmread('IPC_Mexico_preformat');
%data=dlmread('britain_FTSE_100_preformat');
%data=dlmread('HSI_HONGKONG_preformat_fix');
%data=dlmread('DOW_Jones_preformat_fix');
%data=dlmread('markets_final_final_preformat');
% load('markets_final_final_preformat');
load('STOK');

% data=markets_final_final_preformat;
data=MARKET;
% load('Stocktotal');
% data=MARKET;
% fill the blanks
for i=1:size(data,2)
  data(:,i)=fill_blanks(data(:,i));
end

%remove -1 and put 0 instead
temp=find(data==-1);
data(temp)=0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% this part is for removal of two companies from tehran dataset
%% r_row=91;
%% data_new=zeros(size(data,1),size(data,2)-1);
%% data_new(:,1:r_row-1)=data(:,1:r_row-1);
%% data_new(:,r_row:size(data,2)-1)=data(:,r_row+1:size(data,2));

%% clear data;
%% data=data_new;
%% clear data_new r_row;

%% r_row=15;
%% data_new=zeros(size(data,1),size(data,2)-1);
%% data_new(:,1:r_row-1)=data(:,1:r_row-1);
%% data_new(:,r_row:size(data,2)-1)=data(:,r_row+1:size(data,2));

%% clear data;
%% data=data_new;
%% clear data_new r_row;

data=abs(data);  % calculate return
a=log(data);
b=zeros(size(data));
b(1:size(data,1)-1,:)=data(2:size(data,1),:);
b(size(data,1),:)=0;
b=log(b);
c=b-a;
for ii=1:20
c(:,ii)=c(:,ii)*GDP(ii,1);
end;

d=normalize_signal(c);

%p_start 
% p_start=[1 3000 ]
p_start=[1:90:5387];
p_end= [365:90:5752];
correlation_matrix=zeros(length(p_start),size(d,2),size(d,2));
for i=1:length(p_start)
  temp=d(p_start(i):p_end(i),:);
  correlation_matrix(i,:,:)=calculate_correlation(temp);
end
%%clean correlation matrix ( remove Nans )
nans=find(isnan(correlation_matrix)==1);
correlation_matrix(nans)=0;
%% %clean_up
%clear p_start p_end a b c col data filter filter_num i m nans sigma2 temp d;


