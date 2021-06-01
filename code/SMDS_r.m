
% SMDS定位精度与通信范围的关系
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
toc;



