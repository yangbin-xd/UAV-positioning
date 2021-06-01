
clear;clc;close all;

load data1rr;load data2rr;load data3rr;load data4rr;
MAP=real(MAP);MAPP=real(MAPP);SMDSP=real(SMDSP);SMDS=real(SMDS);
Jan_O=MAP(:,1);Feb_O=MAP(:,2);Mar_O=MAP(:,3);Apr_O=MAP(:,4);
May_O=MAP(:,5);Jun_O=MAP(:,6);Jul_O=MAP(:,7);

% Boxplot for the observed temperature from January to December
Temp_O = [Jan_O, Feb_O, Mar_O, Apr_O, May_O, Jun_O, Jul_O];
position_O = 1.3:1:7.3; 
% Define position for 12 Month_O boxplots 
box_O = boxplot(Temp_O,'colors','r','positions',position_O,'width',0.13,'sym',' ','whisker',2);hold on;
set(gca,'XTickLabel',{' '})  % Erase xlabels  
hold on  % Keep the Month_O boxplots on figure overlap the Month_S boxplots  
% Boxplot for the simulated temperature from January to December

load gpsa;GPS2=[GPS1,GPS1(:,1)];
position_A = 1.15:1:7.15;  % Define position for 12 Month_S boxplots 
box_S = boxplot(GPS2,'colors','k','positions',position_A,'width',0.13,'sym',' ','whisker',2);

Jan_S=MAPP(:,1);Feb_S=MAPP(:,2);Mar_S=MAPP(:,3);Apr_S=MAPP(:,4);
May_S=MAPP(:,5);Jun_S=MAPP(:,6);Jul_S=MAPP(:,7);
Temp_S = [Jan_S, Feb_S, Mar_S, Apr_S, May_S, Jun_S, Jul_S];
position_S = 1.45:1:7.45;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_S,'colors','y','positions',position_S,'width',0.13,'sym',' ','whisker',2);

Jan_Z=SMDSP(:,1);Feb_Z=SMDSP(:,2);Mar_Z=SMDSP(:,3);Apr_Z=SMDSP(:,4);
May_Z=SMDSP(:,5);Jun_Z=SMDSP(:,6);Jul_Z=SMDSP(:,7);
Temp_Z = [Jan_Z, Feb_Z, Mar_Z, Apr_Z, May_Z, Jun_Z, Jul_Z];
position_Z = 1.75:1:7.75;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_Z,'colors','g','positions',position_Z,'width',0.13,'sym',' ','whisker',2);  

Jan_T=SMDS(:,1);Feb_T=SMDS(:,2);Mar_T=SMDS(:,3);Apr_T=SMDS(:,4);
May_T=SMDS(:,5);Jun_T=SMDS(:,6);Jul_T=SMDS(:,7);
Temp_T = [Jan_T, Feb_T, Mar_T, Apr_T, May_T, Jun_T, Jul_T];
position_T = 1.6:1:7.6;  % Define position for 12 Month_S boxplots 
box_T = boxplot(Temp_T,'colors','b','positions',position_T,'width',0.13,'sym',' ','whisker',2);  

dd=[11;12];
Jan_D=dd;Feb_D=dd;Mar_D=dd;Apr_D=dd;May_D=dd;
Jun_D=dd;Jul_D=dd;
Temp_D = [Jan_D, Feb_D, Mar_D, Apr_D, May_D, Jun_D, Jul_D];
position_D = 1.45:1:7.45;  % Define position for 12 Month_S boxplots 
box_D = boxplot(Temp_D,'colors','k','positions',position_D,'width',0.13,'sym',' ','whisker',2,'Labels',{'40','50','60','70','80','90','100'});

% legend('MDS-MAP','MDS-MAP(P)','SMDS','SMDS(P)','GPS');
box_vars=findall(gca,'Tag','Box');
hLegend = legend(box_vars([32,39,23,11,16]), {'GPS','MDS-MAP','MDS-MAP(P)','SMDS-Ny','SMDS(P)-Ny-PM'});

% xlabel('r(m)','FontSize',15);
xlabel('$r$(m)','interpreter','latex','fontsize',15);% M
ylabel('Error(m)');
set(gca,'XLim',[0.8 7]);%X轴的数据显示范围
set(gca,'YLim',[0 5]);%X轴的数据显示范围
