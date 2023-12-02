%% normalize for each metric
% for MMS: x = (x-min)/(max-min)
% input: MMS, is a matrix time × metric
MMS2=MMS
MMS2_min=zeros(1,7);
MMS2_max_min=zeros(1,7);
for i=1:7
    MMS2_min(1,i)=min(MMS2(:,i));
    MMS2_max_min(1,i)=max(MMS2(:,i))-min(MMS2(:,i));
    MMS2(:,i)=(MMS2(:,i)-MMS2_min(1,i))/MMS2_max_min(1,i);
end
%% normalize for each metric
% for SR: x = (x-min)/(max-min)
% input: SR, is a matrix time × metric
SR2=SR;
SR2_min=zeros(1,4);
SR2_max_min=zeros(1,4);
for i=1:4
    SR2_min(1,i)=min(SR2(:,i));
    SR2_max_min(1,i)=max(SR2(:,i))-min(SR2(:,i));
    SR2(:,i)=(SR2(:,i)-SR2_min(1,i))/SR2_max_min(1,i);
end
%% normalize for each metric
% for AQI: x = (x-min)/(max-min)
% input: AQI, is a matrix time × metric
AQI2=AQI;
AQI2_min=zeros(1,5);
AQI2_max_min=zeros(1,5);
for i=1:5
    AQI2_min(1,i)=min(AQI2(:,i));
    AQI2_max_min(1,i)=max(AQI2(:,i))-min(AQI2(:,i));
    AQI2(:,i)=(AQI2(:,i)-AQI2_min(1,i))/AQI2_max_min(1,i);
end

%% normalize for each metric
% for MAWI: x = (x-min)/(max-min)
% input: MAWI, is a matrix time × metric
MAWI2=MAWI;
MAWI2_min=zeros(1,5);
MAWI2_max_min=zeros(1,5);
for i=1:5
    MAWI2_min(1,i)=min(MAWI2(:,i));
    MAWI2_max_min(1,i)=max(MAWI2(:,i))-min(MAWI2(:,i));
    MAWI2(:,i)=(MAWI2(:,i)-MAWI2_min(1,i))/MAWI2_max_min(1,i);
end
%% normalize for each metric
% for CAIDA: x = (x-min)/(max-min)
% input: CAIDA, is a matrix time × metric
CAIDA2=CAIDA;
CAIDA2_min=zeros(1,5);
CAIDA2_max_min=zeros(1,5);
for i=1:5
    CAIDA2_min(1,i)=min(CAIDA2(:,i));
    CAIDA2_max_min(1,i)=max(CAIDA2(:,i))-min(CAIDA2(:,i));
    CAIDA2(:,i)=(CAIDA2(:,i)-CAIDA2_min(1,i))/CAIDA2_max_min(1,i);
end
%% normalize for each metric
% for SMD: x = (x-min)/(max-min)
% input: SMD, is a matrix time × metric
% data = SMD_22;
SMD2=[];
SMD = [];
for i=1:38
    if length(find(data(:,i)==0)) == length(data(:,i))
        continue;
    end
    SMD2(:,size(SMD2,2)+1) = data(:,i);
    SMD(:,size(SMD,2)+1) = data(:,i);
end
SMD2_min=zeros(1,size(SMD2,2));
SMD2_max_min=zeros(1,size(SMD2,2));
for i=1:size(SMD2,2)
    SMD2_min(1,i)=min(SMD2(:,i));
    SMD2_max_min(1,i)=max(SMD2(:,i))-min(SMD2(:,i));
    SMD2(:,i)=(SMD2(:,i)-SMD2_min(1,i))/SMD2_max_min(1,i);
end
%% creat input tenor by a parameter T(cycle length)
% T: cycle length; W: window size
% performing only self-embedding transform for the first W slices
T=50;
W=3;
index=1;
M=[];
for i=1:(T*W-T+1)
%     M(:,index,:) = MMS2(i:i+T-1,:);
%     M(:,index,:) = SR2(i:i+T-1,:);
%     M(:,index,:) = AQI2(i:i+T-1,:);
    M(:,index,:) = MAWI2(i:i+T-1,:);
%     M(:,index,:) = CAIDA2(i:i+T-1,:);
%     M(:,index,:) = SMD(i:i+T-1,:);
    index = index+1;
end
W_size=i
for i=T*W+1:T:size(MAWI,1)
%     M(:,index,:)=MMS2(i:i+T-1,:);
%     M(:,index,:) = SR2(i:i+T-1,:);
%     M(:,index,:) = AQI2(i:i+T-1,:);
    M(:,index,:) = MAWI2(i:i+T-1,:);
%     M(:,index,:) = CAIDA2(i:i+T-1,:);
%     M(:,index,:) = SMD(i:i+T-1,:);
    index=index+1;
end
%% creat input tenor by a parameter T(cycle length)
% T: cycle length; W: window size
% performing only self-embedding transform for the first W slices
T=50;
index=1;
W=3;
index=1;
for i = 1:floor(size(MAWI,1)/T)
    Mres2(:,i,:) = MAWI2((i-1)*T+1:i*T,:);
end
% for i=1:(T*W-T+1)
% %     M_mdt(:,index,:) = MMS2(i:i+T-1,:);
% %     M_mdt(:,index,:) = SR2(i:i+T-1,:);
%     M_mdt(:,index,:) = AQI2(i:i+T-1,:);
%     index = index+1;
% end
% W_size=i
% for i=T*W+1:T:5950
% %     M_mdt(:,index,:)=MMS2(i:i+T-1,:);
% %     M_mdt(:,index,:) = SR2(i:i+T-1,:);
%     M_mdt(:,index,:) = AQI2(i:i+T-1,:);
% 
%     index=index+1;
% end
