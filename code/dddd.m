
clear;clc;close all;

load data1d;load data2d;load data3d;load data4d;
MAP=real(MAP);MAPP=real(MAPP);SMDSP=real(SMDSP);SMDS=real(SMDS);
% Jan_O=MAP(:,1);Feb_O=MAP(:,2);Mar_O=MAP(:,3);Apr_O=MAP(:,4);May_O=MAP(:,5);
% Jun_O=MAP(:,6);Jul_O=MAP(:,7);Aug_0=MAP(:,8);Sep_0=MAP(:,9);Oct_0=MAP(:,10);Nov_0=MAP(:,11);
Jan_O=MAP(:,1);Mar_O=MAP(:,3);May_O=MAP(:,5);Jul_O=MAP(:,7);Sep_0=MAP(:,9);Nov_0=MAP(:,11);

% Boxplot for the observed temperature from January to December
Temp_O = [Jan_O, Mar_O, May_O, Jul_O, Sep_0, Nov_0];
position_O = 1.3:1:6.3; 
% Define position for 12 Month_O boxplots 
box_O = boxplot(Temp_O,'colors','r','positions',position_O,'width',0.13,'sym',' ','whisker',2);hold on;
set(gca,'XTickLabel',{' '})  % Erase xlabels  
hold on  % Keep the Month_O boxplots on figure overlap the Month_S boxplots  
% Boxplot for the simulated temperature from January to December

load gpsa;
position_A = 1.15:1:6.15;  % Define position for 12 Month_S boxplots 
box_S = boxplot(GPS1,'colors','k','positions',position_A,'width',0.13,'sym',' ','whisker',2);

% Jan_S=MAPP(:,1);Feb_S=MAPP(:,2);Mar_S=MAPP(:,3);Apr_S=MAPP(:,4);May_S=MAPP(:,5);
% Jun_S=MAPP(:,6);Jul_S=MAPP(:,7);Aug_S=MAPP(:,8);Sep_S=MAPP(:,9);Oct_S=MAPP(:,10);Nov_S=MAPP(:,11);
Jan_S=MAPP(:,1);Mar_S=MAPP(:,2);May_S=MAPP(:,3);Jul_S=MAPP(:,4);Sep_S=MAPP(:,5);Nov_S=MAPP(:,6);
Temp_S = [Jan_S, Mar_S, May_S, Jul_S, Sep_S, Nov_S];
position_S = 1.45:1:6.45;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_S,'colors','y','positions',position_S,'width',0.13,'sym',' ','whisker',2);

% Jan_Z=SMDSP(:,1);Feb_Z=SMDSP(:,2);Mar_Z=SMDSP(:,3);Apr_Z=SMDSP(:,4);May_Z=SMDSP(:,5);
% Jun_Z=SMDSP(:,6);Jul_Z=SMDSP(:,7);Aug_Z=SMDSP(:,8);Sep_Z=SMDSP(:,9);Oct_Z=SMDSP(:,10);Nov_Z=SMDSP(:,11);
Jan_Z=SMDSP(:,1);Mar_Z=SMDSP(:,3);May_Z=SMDSP(:,5);Jul_Z=SMDSP(:,7);Sep_Z=SMDSP(:,9);Nov_Z=SMDSP(:,11);
Temp_Z = [Jan_Z,  Mar_Z,  May_Z,  Jul_Z,  Sep_Z,  Nov_Z];
position_Z = 1.75:1:6.75;  % Define position for 12 Month_S boxplots 
box_Z = boxplot(Temp_Z,'colors','g','positions',position_Z,'width',0.13,'sym',' ','whisker',2);  

% Jan_T=SMDS(:,1);Feb_T=SMDS(:,2);Mar_T=SMDS(:,3);Apr_T=SMDS(:,4);May_T=SMDS(:,5);
% Jun_T=SMDS(:,6);Jul_T=SMDS(:,7);Aug_T=SMDS(:,8);Sep_T=SMDS(:,9);Oct_T=SMDS(:,10);Nov_T=SMDS(:,11);
Jan_T=SMDS(:,1);Mar_T=SMDS(:,3);May_T=SMDS(:,5);Jul_T=SMDS(:,7);Sep_T=SMDS(:,9);Nov_T=SMDS(:,11);
Temp_T = [Jan_T, Mar_T, May_T, Jul_T, Sep_T,  Nov_T];
position_T = 1.6:1:6.6;  % Define position for 12 Month_S boxplots 
box_T = boxplot(Temp_T,'colors','b','positions',position_T,'width',0.13,'sym',' ','whisker',2);  

dd=[11;12];
% Jan_D=dd;Feb_D=dd;Mar_D=dd;Apr_D=dd;May_D=dd;
Jun_D=dd;Jul_D=dd;Aug_D=dd;Sep_D=dd;Oct_D=dd;Nov_D=dd;
Temp_D = [Jun_D, Jul_D, Aug_D, Sep_D, Oct_D, Nov_D];
position_D = 1.45:1:6.45;  % Define position for 12 Month_S boxplots 
box_D = boxplot(Temp_D,'colors','k','positions',position_D,'width',0.13,'sym',' ','whisker',2,'Labels',{'0','2','4','6','8','10'});

% legend('MDS-MAP','MDS-MAP(P)','SMDS','SMDS(P)','GPS');
box_vars=findall(gca,'Tag','Box');
hLegend = legend(box_vars([1,36,21,11,16]), {'GPS','MDS-MAP','MDS-MAP(P)','SMDS-Ny','SMDS(P)-Ny-PM',});

% xlabel('r(m)','FontSize',15);
xlabel('3\sigma_{d}(m)');ylabel('Error(m)');
xlabel('$3\sigma_{d}$(m)','FontSize',16);
xlabel('3$\sigma_{d}$(m)','interpreter','latex','fontsize',15);% M
set(gca,'XLim',[0.8 7]);%X轴的数据显示范围
set(gca,'YLim',[0 5]);%X轴的数据显示范围
