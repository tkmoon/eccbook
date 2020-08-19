function z = doexitchart(data,color)
% Make Exit Charts
% Data should be vertically concatenated exitpath arrays

if(nargin==1)
  color = 'b';
end

data = [[0; 0] data];

% plot([data(1,:), data(1,end:-2:2), data(1,3:2:end)],...
%      [data(2,:), data(2,end:-2:2), data(2,3:2:end)],[color '--'])

p1 = plot(data(1,:),data(2,:),color); hold on;
p2 = plot(data(1,end:-2:2),data(2,end:-2:2),[color '--']);
p3 = plot(data(1,3:2:end),data(2,3:2:end),color);

set(p3,'linewidth',2);
set(p2,'linewidth',2);
% xl = xlabel('Check-to-Bit Info');
% yl = ylabel('Bit-to-Check Info');
xl = xlabel('I_{C\rightarrow B}^{[i]}, I_{C\rightarrow B}^{[i+1]}');
yl = ylabel('I_{B\rightarrow C}^{[i+1]}, I_{B\rightarrow C}^{[i+1]}');

set(xl,'FontSize',15);
set(yl,'FontSize',15);
axis equal;
axis([0 1 0 1]);
%title('EXIT Chart')