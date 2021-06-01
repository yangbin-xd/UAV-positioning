
% MDS-MAP  MDS_MAP(P)  SMDS  SMDS(P)
% 比较四种算法的定位误差与通信范围的关系

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               MDS-MAP                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;tic;

N=50;          %顶点数
eta=3;         %维数
r=50;        %通信范围
de=0.5;        %测距误差
ae=2;        %测角误差
a=0.2;        %GPS配比
total=10000;        %循环次数(time=total/10)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               MDS-MAP(P)                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;

N=50;          %顶点数
eta=3;         %维数
r=50;        %通信范围
de=0.5;        %测距误差
ae=2;        %测角误差
a=0.2;        %GPS配比
total=10000;        %循环次数(time=total)

for r=40:10:100
zp=1;MAPP_e=zeros(1,total);
while zp~=total+1

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

% MDS-MAP(P)
H1=H-eye(N);        %连通矩阵(主对角线元素置零)
c=H1*ones(N,1);        %连通度
nat=(1:N);
for i=1:N-1
    nat=[nat;(1:N)];
end
adj=H.*nat;        %邻接矩阵(包含零元素，使用时需先删除)

% 选定初始点，连通度最大的点
% 建立初始地图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=find(c==max(c));
cu=adj(x,:);
cu(find(cu==0))=[];        %删除零元素
[row,col]=size(cu);        %只用col
n=col;
d1=zeros(n);
for i=1:n
    for j=1:n
        d1(i,j)=dn(cu(i),cu(j));        %局部距离矩阵
    end
end
D=zeros(n,n);
[D,R]=floyd(d1);        %Floyd算法
J=eye(n)-1/n*ones(n);        %中心化矩阵
B=-1/2*J*(D.^2)*J;        %二次中心化
try
[V,D]=eigs(B,eta,'la');        %特征值分解
X1=V(:,1:eta)*D(1:eta,1:eta).^(1/2);        %相对坐标
end

FM=cu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算公共点
% 更新连通度
% 循环选定下一个点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=1;zz=1;
while col~=N&&zz~=N        %还有未定位的点

SAM=zeros(col);        %公共节点
for i=1:col
    cu=adj(FM(i),:);
    cu(find(cu==0))=[];
    same=0;
    for j=1:col
        same=same+length(find(cu==FM(j)));
    end
    c(FM(i))=c(FM(i))-same+1;        %加上自身
    [row1,col1]=size(intersect(cu,FM));
    SAM(i,1:col1)=intersect(cu,FM);
end

% 筛选，至少有四个公共节点的簇头
[row2,col2]=size(FM);
for i=1:col2
    if SAM(i,4)==0
        c(i)=0;
    end
end

% 选择连通度最大的点
xx=find(c(FM)==max(c(FM)));
x1=FM(xx(1));
FM1=adj(x1,:);
FM1(find(FM1==0))=[];        %删除零元素
[row3,col3]=size(FM1);
n=col3;
d1=zeros(n);
for i=1:n
    for j=1:n
        d1(i,j)=d(FM1(i),FM1(j));        %局部距离矩阵
    end
end
D=zeros(n,n);
[D,R]=floyd(d1);        %Floyd算法
J=eye(n)-1/n*ones(n);        %中心化矩阵
B=-1/2*J*(D.^2)*J;        %二次中心化
try
[V,D]=eigs(B,eta,'la');        %特征值分解
X2=V(:,1:eta)*D(1:eta,1:eta).^(1/2);        %相对坐标

% 将X1与X2进行缝合
% 确定公共点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FM=unique(FM);FM1=unique(FM1);
sam=intersect(FM,FM1);

% 获取公共点在X1和X2中坐标
[tf1 loc1]=ismember(sam,FM);        %loc为公共点在X1中位置
[tf2 loc2]=ismember(sam,FM1);        %loc1为公共点在X2中位置
[row4,col4]=size(loc1);        %col4为公共点个数

%最小二乘坐标变换
mu1=0;mu2=0;
P1=X1(loc1(1),:);P2=X2(loc2(1),:);
for i=2:col4
    P1=[P1;X1(loc1(i),:)];        %公共点在X1中坐标
    P2=[P2;X2(loc2(i),:)];        %公共点在X2中坐标
end
mu1=sum(P1)/col4;mu2=sum(P2)/col4;        %求公共点的平均值
P11=P1';P22=P2';
mu1=mu1';mu2=mu2';
P11=P11-mu1*ones(1,col4);P22=P22-mu2*ones(1,col4);
[U,S,V]=svd(P11*P22');
R=U*V';        %旋转矩阵
t=(sum(P1)-sum((R*P2')'))/col4;        %平移向量
X2=R*X2'+t'*ones(1,col3);X2=X2';

%更新X1，公共点取X1与X2的平均值
fm=FM;
fm1=FM1;
FM=[FM,FM1];
FM=unique(FM);        %更新FM    
[row5,col]=size(FM);        %col5为更新后FM中元素个数

X3=zeros(col,eta);
for i=1:col
    [tf3 loc3]=ismember(fm,FM(i));
    [tf4 loc4]=ismember(fm1,FM(i));
    if max(loc3)~=0 && max(loc4)==0
        [m1,n1]=find(loc3==1);n1=n1(1);
        X3(i,:)=X1(n1(1),:);
    elseif max(loc3)==0 && max(loc4)~=0
        [m2,n2]=find(loc4==1);n2=n2(1);
        X3(i,:)=X2(n2,:);
    elseif max(loc3)~=0 && max(loc4)~=0
        [m1,n1]=find(loc3==1);n1=n1(1);
        [m2,n2]=find(loc4==1);n2=n2(1);
        X3(i,:)=(X1(n1,:)+X2(n2,:))/2;
    end
end
X1=X3;
end

c=H*ones(N,1);
z=z+1;
x(1)=x1;
zz=zz+1;
end

if zz==50
    X1=zeros(N,eta);
end

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

MAPP_err=0;GPS_err=0;
for i=1:N
    MAPP_err=MAPP_err+sqrt((X(i,1)-X2(i,1))^2+(X(i,2)-X2(i,2))^2+(X(i,3)-X2(i,3))^2);
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
MAPP_e(zp)=MAPP_err/N;GPS(zp)=GPS_err/N;
if MAPP_e(zp)~=0        %保证有效跑total次
    zp=zp+1;
end
end

j=r/10-3;
MAPP(:,j)=MAPP_e;
end

figure(2);
boxplot(MAPP,'Labels',{'40','50','60','70','80','90','100'},'Whisker',1.5);
xlabel('range(m)');ylabel('error(m)');
save data2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  SMDS(P)                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;

N=50;          %顶点数
eta=3;         %维数
r=50;        %通信范围
de=0.5;        %测距误差
ae=2;        %测角误差
a=0.2;        %GPS配比
total=10000;        %循环次数(time=total/2)

for r=40:10:100
zp=1;
while zp~=total+1

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
H=H-eye(N);        %邻接矩阵
c=H*ones(N,1);        %连通度
d=d-diag(diag(d));        %主对角线元素置零（噪声）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 选定初始点，选择连通度最大的点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=find(c==max(c));
CU=[x(1)];

LN=zeros(N);        %簇矩阵
for i=1:N
    z=find(H(i,:)==max(H(i,:)));
    [row,column]=size(z);
    LN(i,1:column)=z;
end

FM=LN(x(1),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 对选定的簇进行SMDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FM(find(FM==0))=[];        %删除零元素
[row1,column1]=size(FM);
v=zeros(column1,eta);

for i=1:column1
    v(i,:)=X(FM(i),:)-X(x(1),:);        %生成向量（利用实际坐标）
    vn(i,:)=v(i,:)/norm(v(i,:))*dn(FM(i),x(1));
end

K=zeros(column1);
for i=1:column1
    for j=1:column1
        K(i,j)=dot(v(i,:),v(j,:));        %生成矩阵
    end
end
theta=zeros(column1);
for i=1:column1
    for j=1:column1
        cos_theta(i,j)=K(i,j)/norm(v(i,:))/norm(v(j,:));
        theta(i,j)=acos(cos_theta(i,j));
    end
end
theta=abs(theta)/2/pi*360;
angle_noise=normrnd(0,ae/3,column1,column1);
theta=theta+angle_noise;        %测角误差
theta=theta/360*2*pi;

%生成有误差的K
for i=1:column1
    for j=1:column1
        K(i,j)=norm(vn(i,:))*norm(vn(j,:))*cos(theta(i,j));
    end
end

% nystrom SMDS
p=column1-1;        %方阵
A=K(1:p,1:p);
T=K(1:p,p+1:column1);
[V,D]=eigs(A,eta,'la');
V_A=V(:,1:eta)*D(1:eta,1:eta).^(1/2);
% V_T=(sqrt(D^(-1))*V_A'*T)';
V_T=(pinv(V_A)*T)';
V_AT=[V_A;V_T];
XX1=V_AT;
X1=[0,0,0;V_AT];        % X1为初始拓扑图

% 循环选定下一个点
% 更新连通度
% 计算公共点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz=1;zzz=1;
while column1~=N&&zzz~=N        %还有未定位的点

% for zz=1:5

SAM=zeros(column1);        %公共节点
for i=1:column1
    CLU=LN(FM(i),:);
    CLU(find(CLU==0))=[];
    same=0;
    for j=1:column1
        same=same+length(find(CLU==FM(j)));
    end
    c(FM(i))=c(FM(i))-same;
    [row2,column2]=size(intersect(CLU,FM));
    SAM(i,1:column2)=intersect(CLU,FM);
end

% 筛选，至少有四个公共节点的簇头，即SAM中不全零(待补充)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row10,column10]=size(FM);
for i=1:column10
    if SAM(i,3)==0
        c(i)=0;
    end
end

% 选择连通度最大的点
xx=find(c(FM)==max(c(FM)));
x_1=FM(xx(1));
% q(qq)=x_1;
CU=[CU,x_1];
FM1=LN(x_1,:);

FM1(find(FM1==0))=[];        %删除零元素
[row3,column3]=size(FM1);
v=zeros(column3,eta);

for i=1:column3
    v(i,:)=X(FM1(i),:)-X(x_1,:);        %生成向量（利用实际坐标）
    vn(i,:)=v(i,:)/norm(v(i,:))*dn(FM1(i),x_1);
end

K=zeros(column3);
for i=1:column3
    for j=1:column3
        K(i,j)=dot(v(i,:),v(j,:));        %生成矩阵
    end
end
cos_theta=zeros(column3);
theta=zeros(column3);
for i=1:column3
    for j=1:column3
        cos_theta(i,j)=K(i,j)/norm(v(i,:))/norm(v(j,:));
        theta(i,j)=acos(cos_theta(i,j));
    end
end
theta=abs(theta)/2/pi*360;
angle_noise=normrnd(0,ae/3,column3,column3);
theta=theta+angle_noise;        %测角误差
theta=theta/360*2*pi;

%生成有误差的K
for i=1:column3
    for j=1:column3
        K(i,j)=norm(vn(i,:))*norm(vn(j,:))*cos(theta(i,j));
    end
end

% nystrom SMDS
p=column3-1;        %方阵
A=K(1:p,1:p);
T=K(1:p,p+1:column3);
try
[V,D]=eigs(A,eta,'la');
V_A=V(:,1:eta)*D(1:eta,1:eta).^(1/2);
% V_T=(sqrt(D^(-1))*V_A'*T)';
V_T=(pinv(V_A)*T)';
V_AT=[V_A;V_T];
XX2=V_AT;
X2=[0,0,0;V_AT];        % X2为缝合图
% scatter(X2(:,2),X2(:,1),'b');

% 将X1与X2进行缝合
% 确定公共点
z=SAM(xx,:);
z(find(z==0))=[]; 

% 计算公共元素坐标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FM=[x(1),FM];
FM=unique(FM);
[mm,nn]=find(FM==x(1));

[row9,column9]=size(XX1);
XXX1=zeros(row9+1,column9);
XXX1(1:nn-1,:)=XX1(1:nn-1,:);
XXX1(nn,:)=[0,0,0];
XXX1(nn+1:row9+1,:)=XX1(nn:row9,:);

FM1=[x_1,FM1];
FM1=unique(FM1);
[mm1,nn1]=find(FM1==x_1);

[row10,column10]=size(XX2);
XXX2=zeros(row10+1,column10);
XXX2(1:nn1-1,:)=XX2(1:nn1-1,:);
XXX2(nn1,:)=[0,0,0];
XXX2(nn1+1:row10+1,:)=XX2(nn1:row10,:);

X2=XXX2;
if zz==1
    X1=XXX1;
else
end
  
sam=intersect(FM,FM1);

[tf loc]=ismember(sam,FM);
loc(find(loc==0))=[];
[row4,column4]=size(loc);

[tf1 loc1]=ismember(sam,FM1);
loc1(find(loc1==0))=[];
[row5,column5]=size(loc1);

Xp=zeros(column5,eta);
Xq=zeros(column5,eta);
for i=1:column5
    Xp(i,:)=X1(loc(i),:);
    Xq(i,:)=X2(loc1(i),:);
end

n=column5;
Xp_mean=Xp-(sum(Xp)'/n*ones(1,n))';Xp_mean=Xp_mean';
Xq_mean=Xq-(sum(Xq)'/n*ones(1,n))';Xq_mean=Xq_mean';
P=Xp_mean*Xq_mean';
[U,S,V]=svd(P);
R=U*V';t=sum(Xp)'/n-sum(Xq)'/n;
Xq=R*Xq';
Xq=Xq';
s=sum(Xp)/n-sum(Xq)/n;
[vc,vc1]=size(X2);
X2=X2';
X3=zeros(eta,n);
for i=1:vc
    X3(:,i)=R*X2(:,i)+s';
end
X3=X3';
end


% 更新X1
fm1=FM;
fm2=FM1;
FM=[FM,FM1];        
FM=unique(FM);        %保存已定位的点
[row6,column1]=size(FM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XX1=zeros(column1,eta);
for i=1:column1
    [tf2 loc2]=ismember(fm1,FM(i));
    if max(loc2)~=0
        [m1,n1]=find(loc2==1);
        loc2(find(loc2==0))=[];
        XX1(i,:)=X1(n1,:);
    else
        [tf2 loc2]=ismember(fm2,FM(i));
        [m1,n1]=find(loc2==1);
        loc2(find(loc2==0))=[];
        XX1(i,:)=X3(n1,:);
    end
end
X1=XX1;

c=H*ones(N,1);
zz=zz+1;
x(1)=x_1;
zzz=zzz+1;
end

if zz==50
    X1=zeros(N,eta);
end

n=N*a;        %20%的GPS配置
for i=1:n
    Ygps(i,:)=Xgps(i,:);
    Y1(i,:)=X1(i,:);
end

Ygps_mean=Ygps-(sum(Ygps)'/n*ones(1,n))';Ygps_mean=Ygps_mean';
Y1_mean=Y1-(sum(Y1)'/n*ones(1,n))';Y1_mean=Y1_mean';
P=Ygps_mean*Y1_mean';
[U,S,V]=svd(P);
R=U*V';t=sum(Ygps)'/n-sum(Y1)'/n;
Y1=R*Y1';
Y1=Y1';
s=sum(Ygps)/n-sum(Y1)/n;
[vc,vc1]=size(X1);
X1=X1';
X2=zeros(eta,n);
for i=1:vc
    X2(:,i)=R*X1(:,i)+s';
end
X2=X2';

SMDSP_err=0;GPS_err=0;
for i=1:N
    SMDSP_err=SMDSP_err+sqrt((X(i,1)-X2(i,1))^2+(X(i,2)-X2(i,2))^2+(X(i,3)-X2(i,3))^2);
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
SMDSP_e(zp)=SMDSP_err/N;GPS(zp)=GPS_err/N;
if SMDSP_e(zp)<10        %保证有效跑total次
    zp=zp+1;
end
end

j=r/10-3;
SMDSP(:,j)=SMDSP_e;
end

figure(2);
SMDSP=abs(SMDSP);
boxplot(SMDSP,'Labels',{'40','50','60','70','80','90','100'},'Whisker',1.5);
% boxplot(MAP,'Labels',{'40','50','60','70','80','90','100','110','120','130','140','150','160','170','180','190','200'},'Whisker',1.5);
xlabel('range(m)');ylabel('error(m)');
save data3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                    SMDS                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;tic;

N=50;          %顶点数
eta=3;         %维数
r=50;        %通信范围
de=0.5;        %测距误差
ae=2;        %测角误差
a=0.2;        %GPS配比
total=100;        %循环次数(2/217)

for r=40:10:100
zp=1;
while zp~=total+1

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

% SMDS

% SMDS生成K
M=N*(N-1)/2;        % K的大小
for i=1:N-1
    a=N*(i-1)-i*(i-1)/2+1;
    for j=1:N-i
        v(a+j-1,:)=X(i+j,:)-X(i,:);
        if sqrt((X(i+j,1)-X(i,1))^2+(X(i+j,2)-X(i,2))^2+(X(i+j,3)-X(i,3))^2)<=r
            v(a+j-1,:)=X(i+j,:)-X(i,:);
            vn(a+j-1,:)=v(a+j-1,:)/norm(v(a+j-1,:))*dn(i,i+j);        %测距误差
        else
            v(a+j-1,:)=[0,0,0];        %模拟联通缺失
            vn(a+j-1,:)=[0,0,0];
        end
    end
end

K=v(1:M,:)*v(1:M,:)';
%对K模拟测距和测角误差
for i=1:M
    for j=1:M
        cos_theta(i,j)=K(i,j)/norm(v(i,:))/norm(v(j,:));
        theta(i,j)=acos(cos_theta(i,j));
    end
end
theta=abs(theta)/2/pi*360;
angle_noise=normrnd(0,ae/3,M,M);
theta=theta+angle_noise;        %测角误差
theta=theta/360*2*pi;

%生成有误差的K
for i=1:M
    for j=1:M
        K(i,j)=norm(vn(i,:))*norm(vn(j,:))*cos(theta(i,j));
        if abs(K(i,j))<1e10
            K(i,j)=K(i,j);
        else
            K(i,j)=0;
        end    
    end
end

p=N-1;
A=K(1:p,1:p);
T=K(1:p,p+1:M);
[V,D]=eigs(A,eta,'la');
V_A=V(:,1:eta)*D(1:eta,1:eta).^(1/2);
V_T=(pinv(V_A)*T)';
V_AT=[V_A;V_T];
v1=V_AT(1:M,:);
% X1=[0,0,0;v1(1:N-1,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dijkstra算法寻最短路径
%[d1,r]=floyd(d);        %选定1为root,k为leaf

X1=zeros(N,eta);
for k=2:N
    [myd,P]=mydijkstra(d,1,k);
    [f1,f]=size(P);
    for f=1:f-1
        i=P(f);j=P(f+1);
        if i<j
            z=(i-1)*N-i*(i-1)/2+1+j-i-1;
            X1(k,:)=X1(k,:)+v1(z,:);
        elseif i>j
            z=(j-1)*N-j*(j-1)/2+1+i-j-1;        %i与j互换
            X1(k,:)=X1(k,:)-v1(z,:);
        elseif i==j
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GPS配比
n=N*0.2;        
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

% figure(2);
% scatter3(X(:,1),X(:,2),X(:,3),'ko');hold on;
% scatter3(X2(:,1),X2(:,2),X2(:,3),'r+');
% for i=1:N
%     line([X2(i,1),X(i,1)],[X2(i,2),X(i,2)],[X2(i,3),X(i,3)],'linestyle','-','color','r');
%     hold on;
% end
% legend('true','SMDS');
% text(X(:,1)+2,X(:,2)+0.5,X(:,3),num);
% xlabel('x'),ylabel('y'),zlabel('z');
% axis([-60,60,-60,60,-60,60]);
% set(gca,'XTick',-60:20:60);
% set(gca,'YTick',-60:20:60);
% set(gca,'ZTick',-60:20:60);

SMDS_err=0;GPS_err=0;
for i=1:N
    SMDS_err=SMDS_err+sqrt((X(i,1)-X2(i,1))^2+(X(i,2)-X2(i,2))^2+(X(i,3)-X2(i,3))^2);
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
SMDS_e(zp)=SMDS_err/N;
GPS(zp)=GPS_err/N;

if SMDS_e(zp)<10        %保证有效跑total次
    zp=zp+1;
end
end

j=r/10-3;
SMDS(:,j)=SMDS_e;
end

figure(2);
boxplot(SMDS,'Labels',{'40','50','60','70','80','90','100'},'Whisker',1.5);
% boxplot(MAP,'Labels',{'40','50','60','70','80','90','100','110','120','130','140','150','160','170','180','190','200'},'Whisker',1.5);
xlabel('range(m)');ylabel('error(m)');
save data4r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc;



