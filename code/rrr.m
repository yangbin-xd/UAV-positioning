
% MDS-MAP  MDS_MAP(P)  SMDS  SMDS(P)
% �Ƚ������㷨�Ķ�λ�����ͨ�ŷ�Χ�Ĺ�ϵ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               MDS-MAP                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;tic;

N=50;          %������
eta=3;         %ά��
r=50;        %ͨ�ŷ�Χ
de=0.5;        %������
ae=2;        %������
a=0.2;        %GPS���
total=10000;        %ѭ������(time=total/10)

for r=40:10:100
z=1;
while z~=total+1

X=rand(N,eta)*100-50;        %����[-50��50]�Ķ���
nosie1=normrnd(0,5/3,N,1);
nosie2=normrnd(0,5/3,N,1);
nosie3=normrnd(0,10/3,N,1);
noise_GPS=[nosie1,nosie2,nosie3];
Xgps=X+noise_GPS;        %GPS��λ���

num=[];        %ID����
for num1=1:N
    num=[num;num1];
end
num=num2cell(num);
% figure(1);        %�ڵ�����ͼ    
% scatter3(X(:,1),X(:,2),X(:,3),'ko');hold on;
% text(X(:,1)+2,X(:,2)+0.5,X(:,3)+0.5,num);
% xlabel('x'),ylabel('y'),zlabel('z');
% axis([-60,60,-60,60,-60,60]);
% set(gca,'XTick',-60:20:60);
% set(gca,'YTick',-60:20:60);
% set(gca,'ZTick',-60:20:60);

% ͨ�ŷ�Χ�ڽڵ������
d=zeros(N);        %�������
H=zeros(N);        %�ڽӾ���
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

noise=normrnd(0,de/3,N,N);       %������
dn=d+noise;
dn=dn-diag(diag(dn));        %���Խ���Ԫ������
dn=1/2*(dn+dn');        %ƽ�������

% MDS-MAP
D=zeros(N,N);
[D,R]=floyd(dn);        %Floyd�㷨
J=eye(N)-1/N*ones(N);        %���Ļ�����
B=-1/2*J*(D.^2)*J;        %�������Ļ�
try
[V,D]=eigs(B,eta,'la');        %����ֵ�ֽ�
X1=V(:,1:eta)*D(1:eta,1:eta).^(1/2);        %�������

n=N*a;        
for i=1:n
    Ygps(i,:)=Xgps(i,:);
    Y1(i,:)=X1(i,:);
end

%����任(���������X1�任��GPS��λ����Xgps��)
Ygps_mean=Ygps-(sum(Ygps)'/n*ones(1,n))';Ygps_mean=Ygps_mean';
Y1_mean=Y1-(sum(Y1)'/n*ones(1,n))';Y1_mean=Y1_mean';
P=Ygps_mean*Y1_mean';
[U,S,V]=svd(P);        %����ֵ�ֽ�
R=U*V';t=sum(Ygps)'/n-sum(Y1)'/n;
Y1=R*Y1';Y1=Y1';
s=sum(Ygps)/n-sum(Y1)/n;X1=X1';
X2=zeros(eta,N);
for i=1:N
    X2(:,i)=R*X1(:,i)+s';
end
X2=X2';        %�任��ľ�������X2

MAP_err=0;GPS_err=0;
for i=1:N
    MAP_err=MAP_err+sqrt((X(i,1)-X2(i,1))^2+(X(i,2)-X2(i,2))^2+(X(i,3)-X2(i,3))^2);
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
MAP_e(z)=MAP_err/N;GPS(z)=GPS_err/N;
if MAP_e(z)~=0        %��֤��Ч��total��
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

N=50;          %������
eta=3;         %ά��
r=50;        %ͨ�ŷ�Χ
de=0.5;        %������
ae=2;        %������
a=0.2;        %GPS���
total=10000;        %ѭ������(time=total)

for r=40:10:100
zp=1;MAPP_e=zeros(1,total);
while zp~=total+1

X=rand(N,eta)*100-50;        %����[-50��50]�Ķ���
nosie1=normrnd(0,5/3,N,1);
nosie2=normrnd(0,5/3,N,1);
nosie3=normrnd(0,10/3,N,1);
noise_GPS=[nosie1,nosie2,nosie3];
Xgps=X+noise_GPS;        %GPS��λ���

num=[];        %ID����
for num1=1:N
    num=[num;num1];
end
num=num2cell(num);
% figure(1);        %�ڵ�����ͼ    
% scatter3(X(:,1),X(:,2),X(:,3),'ko');hold on;
% text(X(:,1)+2,X(:,2)+0.5,X(:,3)+0.5,num);
% xlabel('x'),ylabel('y'),zlabel('z');
% axis([-60,60,-60,60,-60,60]);
% set(gca,'XTick',-60:20:60);
% set(gca,'YTick',-60:20:60);
% set(gca,'ZTick',-60:20:60);

% ͨ�ŷ�Χ�ڽڵ������
d=zeros(N);        %�������
H=zeros(N);        %�ڽӾ���
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

noise=normrnd(0,de/3,N,N);       %������
dn=d+noise;
dn=dn-diag(diag(dn));        %���Խ���Ԫ������
dn=1/2*(dn+dn');        %ƽ�������

% MDS-MAP(P)
H1=H-eye(N);        %��ͨ����(���Խ���Ԫ������)
c=H1*ones(N,1);        %��ͨ��
nat=(1:N);
for i=1:N-1
    nat=[nat;(1:N)];
end
adj=H.*nat;        %�ڽӾ���(������Ԫ�أ�ʹ��ʱ����ɾ��)

% ѡ����ʼ�㣬��ͨ�����ĵ�
% ������ʼ��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=find(c==max(c));
cu=adj(x,:);
cu(find(cu==0))=[];        %ɾ����Ԫ��
[row,col]=size(cu);        %ֻ��col
n=col;
d1=zeros(n);
for i=1:n
    for j=1:n
        d1(i,j)=dn(cu(i),cu(j));        %�ֲ��������
    end
end
D=zeros(n,n);
[D,R]=floyd(d1);        %Floyd�㷨
J=eye(n)-1/n*ones(n);        %���Ļ�����
B=-1/2*J*(D.^2)*J;        %�������Ļ�
try
[V,D]=eigs(B,eta,'la');        %����ֵ�ֽ�
X1=V(:,1:eta)*D(1:eta,1:eta).^(1/2);        %�������
end

FM=cu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���㹫����
% ������ͨ��
% ѭ��ѡ����һ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=1;zz=1;
while col~=N&&zz~=N        %����δ��λ�ĵ�

SAM=zeros(col);        %�����ڵ�
for i=1:col
    cu=adj(FM(i),:);
    cu(find(cu==0))=[];
    same=0;
    for j=1:col
        same=same+length(find(cu==FM(j)));
    end
    c(FM(i))=c(FM(i))-same+1;        %��������
    [row1,col1]=size(intersect(cu,FM));
    SAM(i,1:col1)=intersect(cu,FM);
end

% ɸѡ���������ĸ������ڵ�Ĵ�ͷ
[row2,col2]=size(FM);
for i=1:col2
    if SAM(i,4)==0
        c(i)=0;
    end
end

% ѡ����ͨ�����ĵ�
xx=find(c(FM)==max(c(FM)));
x1=FM(xx(1));
FM1=adj(x1,:);
FM1(find(FM1==0))=[];        %ɾ����Ԫ��
[row3,col3]=size(FM1);
n=col3;
d1=zeros(n);
for i=1:n
    for j=1:n
        d1(i,j)=d(FM1(i),FM1(j));        %�ֲ��������
    end
end
D=zeros(n,n);
[D,R]=floyd(d1);        %Floyd�㷨
J=eye(n)-1/n*ones(n);        %���Ļ�����
B=-1/2*J*(D.^2)*J;        %�������Ļ�
try
[V,D]=eigs(B,eta,'la');        %����ֵ�ֽ�
X2=V(:,1:eta)*D(1:eta,1:eta).^(1/2);        %�������

% ��X1��X2���з��
% ȷ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FM=unique(FM);FM1=unique(FM1);
sam=intersect(FM,FM1);

% ��ȡ��������X1��X2������
[tf1 loc1]=ismember(sam,FM);        %locΪ��������X1��λ��
[tf2 loc2]=ismember(sam,FM1);        %loc1Ϊ��������X2��λ��
[row4,col4]=size(loc1);        %col4Ϊ���������

%��С��������任
mu1=0;mu2=0;
P1=X1(loc1(1),:);P2=X2(loc2(1),:);
for i=2:col4
    P1=[P1;X1(loc1(i),:)];        %��������X1������
    P2=[P2;X2(loc2(i),:)];        %��������X2������
end
mu1=sum(P1)/col4;mu2=sum(P2)/col4;        %�󹫹����ƽ��ֵ
P11=P1';P22=P2';
mu1=mu1';mu2=mu2';
P11=P11-mu1*ones(1,col4);P22=P22-mu2*ones(1,col4);
[U,S,V]=svd(P11*P22');
R=U*V';        %��ת����
t=(sum(P1)-sum((R*P2')'))/col4;        %ƽ������
X2=R*X2'+t'*ones(1,col3);X2=X2';

%����X1��������ȡX1��X2��ƽ��ֵ
fm=FM;
fm1=FM1;
FM=[FM,FM1];
FM=unique(FM);        %����FM    
[row5,col]=size(FM);        %col5Ϊ���º�FM��Ԫ�ظ���

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

%����任(���������X1�任��GPS��λ����Xgps��)
Ygps_mean=Ygps-(sum(Ygps)'/n*ones(1,n))';Ygps_mean=Ygps_mean';
Y1_mean=Y1-(sum(Y1)'/n*ones(1,n))';Y1_mean=Y1_mean';
P=Ygps_mean*Y1_mean';
[U,S,V]=svd(P);        %����ֵ�ֽ�
R=U*V';t=sum(Ygps)'/n-sum(Y1)'/n;
Y1=R*Y1';Y1=Y1';
s=sum(Ygps)/n-sum(Y1)/n;X1=X1';
X2=zeros(eta,N);
for i=1:N
    X2(:,i)=R*X1(:,i)+s';
end
X2=X2';        %�任��ľ�������X2

MAPP_err=0;GPS_err=0;
for i=1:N
    MAPP_err=MAPP_err+sqrt((X(i,1)-X2(i,1))^2+(X(i,2)-X2(i,2))^2+(X(i,3)-X2(i,3))^2);
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
MAPP_e(zp)=MAPP_err/N;GPS(zp)=GPS_err/N;
if MAPP_e(zp)~=0        %��֤��Ч��total��
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

N=50;          %������
eta=3;         %ά��
r=50;        %ͨ�ŷ�Χ
de=0.5;        %������
ae=2;        %������
a=0.2;        %GPS���
total=10000;        %ѭ������(time=total/2)

for r=40:10:100
zp=1;
while zp~=total+1

X=rand(N,eta)*100-50;        %����[-50��50]�Ķ���
nosie1=normrnd(0,5/3,N,1);
nosie2=normrnd(0,5/3,N,1);
nosie3=normrnd(0,10/3,N,1);
noise_GPS=[nosie1,nosie2,nosie3];
Xgps=X+noise_GPS;        %GPS��λ���

num=[];        %ID����
for num1=1:N
    num=[num;num1];
end
num=num2cell(num);
% figure(1);        %�ڵ�����ͼ    
% scatter3(X(:,1),X(:,2),X(:,3),'ko');hold on;
% text(X(:,1)+2,X(:,2)+0.5,X(:,3)+0.5,num);
% xlabel('x'),ylabel('y'),zlabel('z');
% axis([-60,60,-60,60,-60,60]);
% set(gca,'XTick',-60:20:60);
% set(gca,'YTick',-60:20:60);
% set(gca,'ZTick',-60:20:60);

% ͨ�ŷ�Χ�ڽڵ������
d=zeros(N);        %�������
H=zeros(N);        %�ڽӾ���
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

noise=normrnd(0,de/3,N,N);       %������
dn=d+noise;
dn=dn-diag(diag(dn));        %���Խ���Ԫ������
dn=1/2*(dn+dn');        %ƽ�������

% MDS-MAP
H=H-eye(N);        %�ڽӾ���
c=H*ones(N,1);        %��ͨ��
d=d-diag(diag(d));        %���Խ���Ԫ�����㣨������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ѡ����ʼ�㣬ѡ����ͨ�����ĵ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=find(c==max(c));
CU=[x(1)];

LN=zeros(N);        %�ؾ���
for i=1:N
    z=find(H(i,:)==max(H(i,:)));
    [row,column]=size(z);
    LN(i,1:column)=z;
end

FM=LN(x(1),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ��ѡ���Ĵؽ���SMDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FM(find(FM==0))=[];        %ɾ����Ԫ��
[row1,column1]=size(FM);
v=zeros(column1,eta);

for i=1:column1
    v(i,:)=X(FM(i),:)-X(x(1),:);        %��������������ʵ�����꣩
    vn(i,:)=v(i,:)/norm(v(i,:))*dn(FM(i),x(1));
end

K=zeros(column1);
for i=1:column1
    for j=1:column1
        K(i,j)=dot(v(i,:),v(j,:));        %���ɾ���
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
theta=theta+angle_noise;        %������
theta=theta/360*2*pi;

%����������K
for i=1:column1
    for j=1:column1
        K(i,j)=norm(vn(i,:))*norm(vn(j,:))*cos(theta(i,j));
    end
end

% nystrom SMDS
p=column1-1;        %����
A=K(1:p,1:p);
T=K(1:p,p+1:column1);
[V,D]=eigs(A,eta,'la');
V_A=V(:,1:eta)*D(1:eta,1:eta).^(1/2);
% V_T=(sqrt(D^(-1))*V_A'*T)';
V_T=(pinv(V_A)*T)';
V_AT=[V_A;V_T];
XX1=V_AT;
X1=[0,0,0;V_AT];        % X1Ϊ��ʼ����ͼ

% ѭ��ѡ����һ����
% ������ͨ��
% ���㹫����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz=1;zzz=1;
while column1~=N&&zzz~=N        %����δ��λ�ĵ�

% for zz=1:5

SAM=zeros(column1);        %�����ڵ�
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

% ɸѡ���������ĸ������ڵ�Ĵ�ͷ����SAM�в�ȫ��(������)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row10,column10]=size(FM);
for i=1:column10
    if SAM(i,3)==0
        c(i)=0;
    end
end

% ѡ����ͨ�����ĵ�
xx=find(c(FM)==max(c(FM)));
x_1=FM(xx(1));
% q(qq)=x_1;
CU=[CU,x_1];
FM1=LN(x_1,:);

FM1(find(FM1==0))=[];        %ɾ����Ԫ��
[row3,column3]=size(FM1);
v=zeros(column3,eta);

for i=1:column3
    v(i,:)=X(FM1(i),:)-X(x_1,:);        %��������������ʵ�����꣩
    vn(i,:)=v(i,:)/norm(v(i,:))*dn(FM1(i),x_1);
end

K=zeros(column3);
for i=1:column3
    for j=1:column3
        K(i,j)=dot(v(i,:),v(j,:));        %���ɾ���
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
theta=theta+angle_noise;        %������
theta=theta/360*2*pi;

%����������K
for i=1:column3
    for j=1:column3
        K(i,j)=norm(vn(i,:))*norm(vn(j,:))*cos(theta(i,j));
    end
end

% nystrom SMDS
p=column3-1;        %����
A=K(1:p,1:p);
T=K(1:p,p+1:column3);
try
[V,D]=eigs(A,eta,'la');
V_A=V(:,1:eta)*D(1:eta,1:eta).^(1/2);
% V_T=(sqrt(D^(-1))*V_A'*T)';
V_T=(pinv(V_A)*T)';
V_AT=[V_A;V_T];
XX2=V_AT;
X2=[0,0,0;V_AT];        % X2Ϊ���ͼ
% scatter(X2(:,2),X2(:,1),'b');

% ��X1��X2���з��
% ȷ��������
z=SAM(xx,:);
z(find(z==0))=[]; 

% ���㹫��Ԫ������
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


% ����X1
fm1=FM;
fm2=FM1;
FM=[FM,FM1];        
FM=unique(FM);        %�����Ѷ�λ�ĵ�
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

n=N*a;        %20%��GPS����
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
if SMDSP_e(zp)<10        %��֤��Ч��total��
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

N=50;          %������
eta=3;         %ά��
r=50;        %ͨ�ŷ�Χ
de=0.5;        %������
ae=2;        %������
a=0.2;        %GPS���
total=100;        %ѭ������(2/217)

for r=40:10:100
zp=1;
while zp~=total+1

X=rand(N,eta)*100-50;        %����[-50��50]�Ķ���
nosie1=normrnd(0,5/3,N,1);
nosie2=normrnd(0,5/3,N,1);
nosie3=normrnd(0,10/3,N,1);
noise_GPS=[nosie1,nosie2,nosie3];
Xgps=X+noise_GPS;        %GPS��λ���

num=[];        %ID����
for num1=1:N
    num=[num;num1];
end
num=num2cell(num);
% figure(1);        %�ڵ�����ͼ    
% scatter3(X(:,1),X(:,2),X(:,3),'ko');hold on;
% text(X(:,1)+2,X(:,2)+0.5,X(:,3)+0.5,num);
% xlabel('x'),ylabel('y'),zlabel('z');
% axis([-60,60,-60,60,-60,60]);
% set(gca,'XTick',-60:20:60);
% set(gca,'YTick',-60:20:60);
% set(gca,'ZTick',-60:20:60);

% ͨ�ŷ�Χ�ڽڵ������
d=zeros(N);        %�������
H=zeros(N);        %�ڽӾ���
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

noise=normrnd(0,de/3,N,N);       %������
dn=d+noise;
dn=dn-diag(diag(dn));        %���Խ���Ԫ������
dn=1/2*(dn+dn');        %ƽ�������

% SMDS

% SMDS����K
M=N*(N-1)/2;        % K�Ĵ�С
for i=1:N-1
    a=N*(i-1)-i*(i-1)/2+1;
    for j=1:N-i
        v(a+j-1,:)=X(i+j,:)-X(i,:);
        if sqrt((X(i+j,1)-X(i,1))^2+(X(i+j,2)-X(i,2))^2+(X(i+j,3)-X(i,3))^2)<=r
            v(a+j-1,:)=X(i+j,:)-X(i,:);
            vn(a+j-1,:)=v(a+j-1,:)/norm(v(a+j-1,:))*dn(i,i+j);        %������
        else
            v(a+j-1,:)=[0,0,0];        %ģ����ͨȱʧ
            vn(a+j-1,:)=[0,0,0];
        end
    end
end

K=v(1:M,:)*v(1:M,:)';
%��Kģ����Ͳ�����
for i=1:M
    for j=1:M
        cos_theta(i,j)=K(i,j)/norm(v(i,:))/norm(v(j,:));
        theta(i,j)=acos(cos_theta(i,j));
    end
end
theta=abs(theta)/2/pi*360;
angle_noise=normrnd(0,ae/3,M,M);
theta=theta+angle_noise;        %������
theta=theta/360*2*pi;

%����������K
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

% Dijkstra�㷨Ѱ���·��
%[d1,r]=floyd(d);        %ѡ��1Ϊroot,kΪleaf

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
            z=(j-1)*N-j*(j-1)/2+1+i-j-1;        %i��j����
            X1(k,:)=X1(k,:)-v1(z,:);
        elseif i==j
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GPS���
n=N*0.2;        
for i=1:n
    Ygps(i,:)=Xgps(i,:);
    Y1(i,:)=X1(i,:);
end

%����任(���������X1�任��GPS��λ����Xgps��)
Ygps_mean=Ygps-(sum(Ygps)'/n*ones(1,n))';Ygps_mean=Ygps_mean';
Y1_mean=Y1-(sum(Y1)'/n*ones(1,n))';Y1_mean=Y1_mean';
P=Ygps_mean*Y1_mean';
[U,S,V]=svd(P);        %����ֵ�ֽ�
R=U*V';t=sum(Ygps)'/n-sum(Y1)'/n;
Y1=R*Y1';Y1=Y1';
s=sum(Ygps)/n-sum(Y1)/n;X1=X1';
X2=zeros(eta,N);
for i=1:N
    X2(:,i)=R*X1(:,i)+s';
end
X2=X2';        %�任��ľ�������X2

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

if SMDS_e(zp)<10        %��֤��Ч��total��
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



