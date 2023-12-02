%orig_data is original data in real-time monitoring phase£¬esti_data is the recovered data
%matrix, row is time, column is metric
function [p_orig_NMAEs,p_orig_COSes]=getPerformanceNC_orign(orig_data, esti_data)

    %NMAE
    J=size(orig_data,2);
    p_orig_NMAEs=zeros(1,J);
    for j=1:J
        p_orig_NMAEs(1,j)=sum(abs(esti_data(:,j)-orig_data(:,j)))/sum(abs(orig_data(:,j)));
    end
    %COS
    p_orig_COSes=zeros(1,J);
    for j=1:J
        p_orig_COSes(1,j)=sum(esti_data(:,j).*orig_data(:,j))/(sqrt(sum(esti_data(:,j).^2)) * sqrt(sum(orig_data(:,j).^2)));
    end

    