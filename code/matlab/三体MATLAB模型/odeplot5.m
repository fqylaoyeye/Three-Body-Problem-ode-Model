function odeplot5
% 三星系统
close all;

G = 6.67e-11; %万有引力常数
m1 = 1/G;     %小质量星体
m2 = m1*5;    %中质量星体
m3 = m1*15;    %大质量星体

% 求解区
%-------------------------------------------------------------
u1 = [-0.01 0.1 0.1];
u2 = [10 5 -0.8];
u3 = [-0.15 0.15 0.12];
u4 = [4 8 -3.1];

u5 = ([0 0 0]-m1*u1-m2*u3)/m3;
u6 = ([0 0 0]-m1*u2-m2*u4)/m3;

x0 = [u1 u2 u3 u4 u5 u6];
% 定步长控制
tspan = [0:0.00001:10];
% 求解微分方程
[t,x] = ode45(@(t,x) guidao(t,x,m1,m2,m3),tspan,x0);

% 绘图区
%-------------------------------------------------------------
figure('Color','k','Units','normalized','Position',[0.1 0.1 0.8 0.8])
% 绘制三星
movep1 = scatter3(u1(1),u1(2),u1(3),'SizeData',400,'MarkerFaceColor','r','MarkerEdgeColor','none');
hold on
movep2 = scatter3(u3(1),u3(2),u3(3),'SizeData',800,'MarkerFaceColor','g','MarkerEdgeColor','none');
movep3 = scatter3(u5(1),u5(2),u5(3),'SizeData',3000,'MarkerFaceColor','b','MarkerEdgeColor','none');

% 添加轨迹线
h1 = animatedline('MaximumNumPoints',fix(length(tspan)),'Color','r','LineWidth',1.5);
h2 = animatedline('MaximumNumPoints',fix(length(tspan)),'Color','g','LineWidth',1.5);
h3 = animatedline('MaximumNumPoints',fix(length(tspan)),'Color','b','LineWidth',1.5);

ax = h1.Parent;
ax.Color = 'none';
ax.XAxis.Color = 'y';
ax.XAxis.FontName = 'times';
ax.YAxis.Color = 'y';
ax.ZAxis.Color = 'y';
ax.View = [22 33];
% view(2)
axis equal
axis off

for k = 1:length(t)
    if k<10
        pause(0.1)
    end
    
    addpoints(h1,x(k,1),x(k,2),x(k,3));
    movep1.XData = x(k,1);
    movep1.YData = x(k,2);
    movep1.ZData = x(k,3);

    addpoints(h2,x(k,7),x(k,8),x(k,9));
    movep2.XData = x(k,7);
    movep2.YData = x(k,8);
    movep2.ZData = x(k,9);
    
    addpoints(h3,x(k,13),x(k,14),x(k,15));
    movep3.XData = x(k,13);
    movep3.YData = x(k,14);
    movep3.ZData = x(k,15);
%     ax.View(1) = ax.View(1)+0.01;
%     ax.View(2) = ax.View(2)+0.001;
    drawnow limitrate
end

function dx=guidao(t,x,m1,m2,m3)
% 三星系统的运动微分方程
G = 6.67e-11; %万有引力常数

dx=zeros(18,1);

u1 = x(1:3);
u2 = x(4:6);
u3 = x(7:9);
u4 = x(10:12);
u5 = x(13:15);
u6 = x(16:18);

u1x = u2;
u2x = G*(m2*(u3-u1)/(norm(u3-u1))^3+m3*(u5-u1)/(norm(u5-u1))^3);
u3x = u4;
u4x = G*(m3*(u5-u3)/(norm(u5-u3))^3+m1*(u1-u3)/(norm(u1-u3))^3);
u5x = u6;
u6x = G*(m1*(u1-u5)/(norm(u1-u5))^3+m2*(u3-u5)/(norm(u3-u5))^3);

dx = [u1x;u2x;u3x;u4x;u5x;u6x];



