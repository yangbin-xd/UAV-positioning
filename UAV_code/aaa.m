
clear;clc;close all;

load data1a;load data2a;load data3a;load data4a;
MAP=real(MAP);MAPP=real(MAPP);SMDSP=real(SMDSP);SMDS=real(SMDS);
Jan_O=MAP(:,1);Feb_O=MAP(:,2);Mar_O=MAP(:,3);Apr_O=MAP(:,4);May_O=MAP(:,5);
Jun_O=MAP(:,6);Jul_O=MAP(:,7);Aug_0=MAP(:,8);Sep_0=MAP(:,9);Oct_0=MAP(:,10);

% Boxplot for the observed temperature from January to December
%Temp_O = [Jan_O, Feb_O, Mar_O, Apr_O, May_O, Jun_O, Jul_O, Aug_0, Sep_0, Oct_0];
Temp_O = [Feb_O,  Apr_O,Jun_O, Aug_0,  Oct_0];
position_O = 1.3:1:5.3; 
% Define position for 12 Month_O boxplots 
box_O = boxplot(Temp_O,'colors','r','positions',position_O,'width',0.13,'sym',' ','whisker',2);hold on; 
set(gca,'XTickLabel',{' '})  % Erase xlabels  
hold on  % Keep the Month_O boxplots on figure overlap the Month_S boxplots  
% Boxplot for the simulated temperature from January to December

Jan_Z=SMDSP(:,1);Feb_Z=SMDSP(:,2);Mar_Z=SMDSP(:,3);Apr_Z=SMDSP(:,4);May_Z=SMDSP(:,5);
Jun_Z=SMDSP(:,6);Jul_Z=SMDSP(:,7);Aug_Z=SMDSP(:,8);Sep_Z=SMDSP(:,9);Oct_Z=SMDSP(:,10);
%Temp_Z = [Jan_Z, Feb_Z, Mar_Z, Apr_Z, May_Z, Jun_Z, Jul_Z, Aug_Z, Sep_Z, Oct_Z];
Temp_Z = [ Feb_Z,Apr_Z,Jun_Z,  Aug_Z, Oct_Z];
position_Z = 1.75:1:5.75;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_Z,'colors','g','positions',position_Z,'width',0.13,'sym',' ','whisker',2,'Labels',{'0.2','0.4','0.6','0.8','1'});  

dd=[11;12];
Jan_D=dd;Feb_D=dd;Mar_D=dd;Apr_D=dd;May_D=dd;
Jun_D=dd;Jul_D=dd;Aug_D=dd;Sep_D=dd;Oct_D=dd;
%Temp_D = [Jan_D, Feb_D, Mar_D, Apr_D, May_D, Jun_D, Jul_D, Aug_D, Sep_D, Oct_D];
Temp_D = [ Feb_D, Apr_D,Jun_D, Aug_D,  Oct_D];
position_D = 1.15:1:5.15;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_D,'colors','k','positions',position_D,'width',0.13,'sym',' ','whisker',2,'Labels',{'0.2','0.4','0.6','0.8','1'});

load datagps;GPS=GPS1(:,1:5);
position_A = 1.15:1:5.15;  % Define position for 12 Month_S boxplots 
box_S = boxplot(GPS,'colors','k','positions',position_A,'width',0.13,'sym',' ','whisker',2,'Labels',{'0.2','0.4','0.6','0.8','1'});

Jan_T=SMDS(:,1);Feb_T=SMDS(:,2);Mar_T=SMDS(:,3);Apr_T=SMDS(:,4);May_T=SMDS(:,5);
Jun_T=SMDS(:,6);Jul_T=SMDS(:,7);Aug_T=SMDS(:,8);Sep_T=SMDS(:,9);Oct_T=SMDS(:,10);
%Temp_T = [Jan_T, Feb_T, Mar_T, Apr_T, May_T, Jun_T, Jul_T, Aug_T, Sep_T, Oct_T];
Temp_T = [ Feb_T, Apr_T, Jun_T,Aug_T, Oct_T];
position_T = 1.6:1:5.6;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_T,'colors','b','positions',position_T,'width',0.13,'sym',' ','whisker',2,'Labels',{'0.2','0.4','0.6','0.8','1'});  

Jan_S=MAPP(:,1);Feb_S=MAPP(:,2);Mar_S=MAPP(:,3);Apr_S=MAPP(:,4);May_S=MAPP(:,5);
Jun_S=MAPP(:,6);Jul_S=MAPP(:,7);Aug_S=MAPP(:,8);Sep_S=MAPP(:,9);Oct_S=MAPP(:,10);
%Temp_S = [Jan_S, Feb_S, Mar_S, Apr_S, May_S, Jun_S, Jul_S, Aug_S, Sep_S, Oct_S];
Temp_S = [ Feb_S, Apr_S, Jun_S, Aug_S, Oct_S];
position_S = 1.45:1:5.45;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_S,'colors','y','positions',position_S,'width',0.13,'sym',' ','whisker',2,'Labels',{'0.2','0.4','0.6','0.8','1'});

% legend('MDS-MAP','MDS-MAP(P)','SMDS','SMDS(P)','GPS');
box_vars=findall(gca,'Tag','Box');
hLegend = legend(box_vars([20,30,1,10,22]), {'GPS','MDS-MAP','MDS-MAP(P)','SMDS-Ny','SMDS(P)-Ny-PM'});

xlabel('\rho','FontSize',16);ylabel('Error(m)');
set(gca,'XLim',[0.8 6]);%X轴的数据显示范围
set(gca,'YLim',[0 5]);%X轴的数据显示范围

