function dotrajectory(data)

plot(data(1,3:2:end),'-')
hold on;
plot(data(2,3:2:end),'--')

legend('Chk to Bit','Bit to Chk')
xl = xlabel('Iteration'); 
yl = ylabel('Information')
set(xl,'fontsize',15);
set(yl,'fontsize',15);
