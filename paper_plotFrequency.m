
IDX = W_size+1:J;
Num4slice = sum(sum(omega(:,IDX,:),3),1);

figure
plot(1:length(IDX), Num4slice, 'O-', 'LineWidth',3,'MarkerSize',10);
grid on;
ylim([0 I*K])
% ylabel('Number of samples')
ylabel('Frequency')
xlabel('Cycle in network telemetry')
set(gca,'FontName','Times New Roman' ,'FontSize',24, 'FontWeight','bold');

figure
plot(1:length(IDX), rs(IDX), '-', 'LineWidth',3,'MarkerSize',10);
grid on;
ylabel('Estimated rank')
xlabel('Cycle in network telemetry')
set(gca,'FontName','Times New Roman' ,'FontSize',24, 'FontWeight','bold');

% figure
% plot(1:length(IDX), rss(:,IDX), '-', 'LineWidth',2,'MarkerSize',5)
% grid on;
% % ylim([0 200])
% ylabel('Estimated rank')
% xlabel('Cycle in network telemetry')
% legend('Memory usage','CPU usage','ReceiveBytes','TransmitBytes','Latency','SuccessRate','Workload')
% set(gca,'FontName','Times New Roman' ,'FontSize',24, 'FontWeight','bold');