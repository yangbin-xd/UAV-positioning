
% MDS-MAP(P)��λ������ͨ�ŷ�Χ�Ĺ�ϵ
clear;clc;close all;tic;

N=50;          %������
eta=3;         %ά��
r=50;        %ͨ�ŷ�Χ
de=0.5;        %������
ae=2;        %������
a=0.1;        %GPS���
total=10;        %ѭ������(time=total)

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
toc;




