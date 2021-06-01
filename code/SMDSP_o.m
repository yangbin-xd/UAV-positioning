
% SMDS(P)定位精度与通信范围的关系
clear;clc;close all;tic;

N=50;          %顶点数
eta=3;         %维数
r=80;        %通信范围
de=0.5;        %测距误差
ae=2;        %测角误差
a=0.2;        %GPS配比
total=100;        %循环次数(100/17)

for ae=0:4:20
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

j=ae/4+1;j=round(j);
SMDSP(:,j)=SMDSP_e;
end

SMDSP=real(SMDSP);
figure(2);
boxplot(SMDSP,'Labels',{'0','4','8','12','16','20'},'Whisker',2,'sym',' ');
xlabel('3\sigma_{\theta}({\circ})');ylabel('error(m)');
save data3o;
toc;
