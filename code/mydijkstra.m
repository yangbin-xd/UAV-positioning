function [mydistance,mypath]=mydijkstra(a,sb,db);
% 输入：a—邻接矩阵(aij)是指i到j之间的距离，可以是有向的
% sb—起点的标号, db—终点的标号
% 输出：mydistance—最短路的距离, mypath—最短路的路径
n=size(a,1); visited(1:n) = 0;
distance(1:n) = inf; % 保存起点到各顶点的最短距离
distance(sb) = 0; parent(1:n) = 0;
for i = 1: n-1
temp=distance;
id1=find(visited==1); %查找已经标号的点
temp(id1)=inf; %已标号点的距离换成无穷
[t, u] = min(temp); %找标号值最小的顶点
visited(u) = 1; %标记已经标号的顶点
id2=find(visited==0); %查找未标号的顶点
for v = id2
if a(u, v) + distance(u) < distance(v)
distance(v) = distance(u) + a(u, v); %修改标号值
parent(v) = u;
end
end
end
mypath = [];
if parent(db) ~= 0 %如果存在路!
t = db; mypath = [db];
while t ~= sb
p = parent(t);
mypath = [p mypath];
t = p;
end
end
mydistance = distance(db);
return