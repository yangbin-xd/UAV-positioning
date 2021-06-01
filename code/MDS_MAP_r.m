
% MDS-MAP��λ������ͨ�ŷ�Χ�Ĺ�ϵ
clear;clc;close all;tic;

N=50;          %������
eta=3;         %ά��
r=50;        %ͨ�ŷ�Χ
de=0.5;        %������
ae=2;        %������
a=0.1;        %GPS���
total=100;        %ѭ������(time=total/10)

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
toc;





