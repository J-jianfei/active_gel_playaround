

for i = 1:length(U_data);
t = dtsave*(0:length(U_data)-1);

f = figure;
f.Position = [f.Position(1),f.Position(2),2*f.Position(3),f.Position(4)];
subplot(1,2,1);
plot(x,phi_data(i,:),'-k','LineWidth',2.0);
xticks(0:0.2:1);
yticks(0:0.2:1);
xlim([0 1]);
ylim([0 1]);
xlabel('x');
ylabel('\phi');
title('Polymer volume fraction')
set(gca,'FontSize',20);
set(gca,'LineWidth',2.0);

subplot(1,2,2)

plot(t,U_data,'--k','LineWidth',2.0)
ylim([8 16])
xlim([-0.01,0.31])
hold on;
scatter(t(i),U_data(i),50,'or','LineWidth',2);
ylabel('Gel layer velocity U');
xlabel('t');
title('Instantaneous gel velocity')
set(gca,'FontSize',20);
set(gca,'LineWidth',2.0);

frame(i) = getframe(gcf);
pause(0.1)
close(f);

end