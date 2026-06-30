clc
figure(2);

% com1 = 1;
% com2 = 80;
% plot(squeeze(S(com1,com1,:)),'r'); 
% hold on; 
% stem(squeeze(S(com2,com2,:)),'k');
% 
% plot(idxS, '-+r')
% hold on 
% % plot(idxK,'-+r')
% plot(idxS1,'-*g')
% % legend('Index', 'Index K', 'Index S1');
% xlabel('Sample Number');
% ylabel('Value');
% title('Plot of Indices');

% hold on
% plot(Features(1,:),'--r');
% plot(Features(2,:),'--g');
% plot(Features(3,:),'-*k');
% plot(Features(4,:),'-+b');
% plot(Features(5,:),'-*m');
% 
% legend;
% xlabel('Sample Number');
% ylabel('Feature Value');
% title('Plot of Features');
 
for b = 2
    extracted_signal = squeeze(reconstructed_data_all(b, 2, view_range));
    plot(t(view_range), extracted_signal, 'Color', colors(b), 'LineWidth', 1.2);
    hold on;
end
plot(t(view_range), data(2, view_range), '--r', 'LineWidth', 1);