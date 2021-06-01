

% GPS error
clear;clc;
N=50;eta=3;

for j=1:10000
X=rand(N,eta)*100-50;        %生成[-50，50]的顶点
nosie1=normrnd(0,5/3,N,1);
nosie2=normrnd(0,5/3,N,1);
nosie3=normrnd(0,10/3,N,1);
noise_GPS=[nosie1,nosie2,nosie3];
Xgps=X+noise_GPS;        %GPS定位误差

GPS_err=0;
for i=1:N
    GPS_err=GPS_err+sqrt((X(i,1)-Xgps(i,1))^2+(X(i,2)-Xgps(i,2))^2+(X(i,3)-Xgps(i,3))^2);
end
j=round(j);
GPS(j)=GPS_err/N;
end

GPS1=[GPS',GPS',GPS',GPS',GPS',GPS'];
boxplot(GPS1,'Labels',{'0','4','8','12','16','20'},'Whisker',2,'sym',' ');
save gpsa;
