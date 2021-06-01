% floyd.m
% 采用floyd算法计算图a中每对顶点最短路
% d是矩离矩阵
% r是路由矩阵
function [d,r]=floyd(a)
n=size(a,1);
% 初始化距离矩阵
d=a;
% 初始化路由矩阵
for i=1:n
    for j=1:n
        r(i,j)=j;
    end 
end 
r;

% Floyd算法开始
for k=1:n
    for i=1:n
        for j=1:n
            if d(i,k)+d(k,j)<d(i,j)
                d(i,j)=d(i,k)+d(k,j);
                r(i,j)=r(i,k);
            end 
        end 
    end
    k;
    d;
    r;
end

