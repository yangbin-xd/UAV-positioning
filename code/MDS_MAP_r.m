
% MDS-MAP定位精度与通信范围的关系
clear;clc;close all;tic;

N=50;          %顶点数
eta=3;         %维数
r=50;        %通信范围
de=0.5;        %测距误差
ae=2;        %测角误差
a=0.1;        %GPS配比
total=100;        %循环次数(time=total/10)

for r=40:10:100
z=1;
while z~=total+1

X=rand(N,eta)*100-50;        %生成[-50，50]的顶点
nosie1=normrnd(0,5/3,N,1);
nosie2=normrnd(0,5/3,N,1);
nosie3=normrnd(0,10/3,N,1);
noise_GPS=[nosie1,nosie2,nosie3];
Xgps=X+noise_GPS;        %GPS定位误差

num=[];        %ID集合
for num1=1:N
    num=[num;num1];
end
num=num2cell(num);
% figure(1);        %节点拓扑图    
% scatter3(X(:,1),X(:,2),X(:,3),'ko');hold on;
% text(X(:,1)+2,X(:,2)+0.5,X(:,3)+0.5,num);
% xlabel('x'),ylabel('y'),zlabel('z');
% axis([-60,60,-60,60,-60,60]);
% set(gca,'XTick',-60:20:60);
% set(gca,'YTick',-60:20:60);
% set(gca,'ZTick',-60:20:60);

% 通信范围内节点间连线
d=zeros(N);        %距离矩阵
H=zeros(N);        %邻接矩阵
for i=1:N
    for j=1:N
        if sqrt((X(i,1)-X(j,1))^2+(X(i,2)-X(j,2))^2+(X(i,3)-X(j,3))^2)<=r
            d(i,j)=sqrt((X(i,1)-X(j,1))^2+(X(i,2)-X(j,2))^2+(X(i,3)-X(j,3))^2);
            H(i,j)=1;
%             line([X(i,1),X(j,1)],[X(i,2),X(j,2)],[X(i,3),X(j,3)],'linestyle','-','color','r');
%             hold on;
        else
            d(i,j)=inf;
        end
    end
end
% legend('UAV','link');

noise=normrnd(0,de/3,N,N);       %测距误差
dn=d+noise;
dn=dn-diag(diag(dn));        %主对角线元素置零
dn=1/2*(dn+dn');        %平均测距结果

% MDS-MAP
D=zeros(N,N);
[D,R]=floyd(dn);        %Floyd算法
J=eye(N)-1/N*ones(N);        %中心化矩阵
B=-1/2*J*(D.^2)*J;        %二次中心化
try
[V,D]=eigs(B,eta,'la');        %特征值分解
X1=V(:,1:eta)*D(1:eta,1:eta).^(1/2);        %相对坐标

n=N*a;        
for i=1:n
    Ygps(i,:)=Xgps(i,:);
    Y1(i,:)=X1(i,:);
end

%坐标变换(将相对坐标X1变换到GPS定位坐标Xgps上)
Ygps_mean=Ygps-(sum(Ygps)'/n*ones(1,n))';Ygps_mean=Ygps_mean';
Y1_mean=Y1-(sum(Y1)'/n*ones(1,n))';Y1_mean=Y1_mean';
P=Ygps_mean*Y1_mean';
[U,S,V]=svd(P);        %奇异值分解
R=U*V';t=sum(Ygps)'/n-sum(Y1)'/n;
Y1=R*Y1';Y1=Y1';
s=sum(Ygps)/n-sum(Y1)/n;X1=X1';
X2=zeros(eta,N);
for i=1:N
    X2(:,i)=R*X1(:,i)+s';
end
X2=X2';        %变换后的绝对坐标X2

MAP_err=0;GPS_err=0;
for i=1:N
    MAP_err=MAP_err+sqrt((X(i,1)-X2(i,1))^2+(X(i,2)-X2(i,2))^2+(X(i,3)-X2(i,3))^2);
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
MAP_e(z)=MAP_err/N;GPS(z)=GPS_err/N;
if MAP_e(z)~=0        %保证有效跑total次
    z=z+1;
end
end
end

j=r/10-3;
MAP(:,j)=MAP_e;
end

figure(2);
boxplot(MAP,'Labels',{'40','50','60','70','80','90','100'},'Whisker',1.5);
% boxplot(MAP,'Labels',{'40','50','60','70','80','90','100','110','120','130','140','150','160','170','180','190','200'},'Whisker',1.5);
xlabel('range(m)');ylabel('error(m)');
save data1;
toc;





