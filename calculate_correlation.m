function [ correlation_matrix ] = calculate_correlation( input_data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% normalize input matrix
  input_col=size(input_data,2);
  input_data=normalize_signal(input_data);
  
  correlation_matrix=zeros(input_col,input_col);
  finites=isfinite(input_data);
  for i=1:input_col
    for j=1:input_col
      filter=finites(:,i).*finites(:,j);
      filter_num=find(filter==1);
      correlation_matrix(i,j)=sum(input_data(filter_num,i).*input_data(filter_num,j));
      correlation_matrix(i,j)=correlation_matrix(i,j)/length(filter_num);
      correlation_matrix(j,i)=correlation_matrix(i,j);
    end
  end

end

