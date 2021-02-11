%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code for COVID-19 pandemic prediction with the SIRV model proposed
% in the following work:
% Modeling the Effect of Population-Wide Vaccination on the Evolution of COVID-19 epidemic in Canada
% by Intissar Harizi, Soulaimane Berkane, and Abdelhamid Tayebi
% Soulaimane Berkane, Jan 15, 2021
% Contact: soulaimane.berkane@uqo.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
%% Selection of the population
population='Canada'; 
vaccine_rate_1000=[0.5 1 2]; % possible vaccination rates per 1,000 population
alpha=[0.95 0.6]; %possible vaccine efficacies
startDate= datenum('17-Jul-2020'); %start day of the data
endDate=datenum('08-Jan-2021'); %end day of the data, start day of the prediction
predictionEndDate=datenum('31-Dec-2021'); %end day of the prediction
horizon=predictionEndDate-endDate; % prediction horizon
tmonths=startDate:30:(predictionEndDate);
%% Model parameters
%beta: rate of infection
%gamma:  rate of removal
%alpha: vaccine efficacy
%mu: mortality rate
%N: total number of population
%I0: initial infected sub-population at startDate
%R0: initial removed sub-population at startDate
if strcmp(population,'Canada')
    beta = 9.76e-2; 
    gamma =7.81e-2;
    mu=1.64e-2;
    N =38e6;
    I0 = 4141;
    R0 = 105528;
elseif strcmp(population,'Ontario')
    beta = 10.67e-2; 
    gamma =8.89e-2;
    mu=1.4e-2;
    N =14.7e6;
    I0 = 1366;
    R0 = 35908;
elseif strcmp(population,'Quebec')
    beta = 10.81e-2; 
    gamma =9.06e-2;
    mu=1.96e-2;
    N =8.6e6;
    I0 = 1556;
    R0 = 55586;
elseif strcmp(population,'British Columbia')
    beta = 11.25e-2; 
    gamma =8.83e-2;
    mu=1.59e-2;
    N =5.15e6;
    I0 = 207;
    R0 = 2991;
elseif strcmp(population,'Alberta')
    beta = 8.36e-2; 
    gamma =6.34e-2;
    mu=1.14e-2;
    N =4.4e6;
    I0 = 859;
    R0 = 8360;
elseif strcmp(population,'Saskatchewan')
    beta = 7.10e-2; 
    gamma =4.93e-2;
    mu=1.25e-2;
    N =1.13e6;
    I0 = 128;
    R0 = 808;
elseif strcmp(population,'Manitoba')
    beta = 18.03e-2; 
    gamma =14.26e-2;
    mu=3.35e-2;
    N =1.38e6;
    I0 = 11;
    R0 = 325;
end
%% Simulation of the model
vaccine_rate=vaccine_rate_1000*N/1000; % possible vaccination rates
I0=I0/N; R0=R0/N;V0=0; % state variables are normalized
dt=1/24;% time interval between two integration samples (1 hour)
t0=endDate-startDate+1;
tf=t0+horizon;
t=startDate:dt:predictionEndDate;
vaccine_start_day=t0;
% Simulation of the model without vaccination
[S,I,R,~] = SIRV_model(0,beta,gamma,I0,R0,V0,tf-1,dt,zeros(1,(tf-1)/dt));
for i=2:length(I)
   newCases(i-1)=(1/dt)*(I(i)+R(i)-I(i-1)-R(i-1));
   newDeaths(i-1)=(1/dt)*mu*(R(i)-R(i-1));
end
sz=0.5;
figure(1);hold on;plot(t(1:end-1),I*N,'LineWidth',sz,'DisplayName','Without Vaccination');grid on;box on;
ax1=get(gcf,'CurrentAxes');
figure(2);hold on;plot(t(1:end-1),R*N,'LineWidth',sz,'DisplayName','Without Vaccination');grid on;box on;
ax2=get(gcf,'CurrentAxes');
figure(3);hold on;plot(t(1:end-2),newCases*N,'LineWidth',sz,'DisplayName','Without Vaccination');grid on;box on;
ax3=get(gcf,'CurrentAxes');
figure(4);hold on;plot(t(1:end-2),newDeaths*N,'LineWidth',sz,'DisplayName','Without Vaccination');grid on;box on;
ax4=get(gcf,'CurrentAxes');
% Simulation of the model with vaccination at different rates
vaccine_end_day=tf;
for k=1:length(alpha)
    sz=0.5;
    set(ax1,'ColorOrderIndex',2)
    set(ax2,'ColorOrderIndex',2)
    set(ax3,'ColorOrderIndex',2)
    set(ax4,'ColorOrderIndex',2)
    for n=1:length(vaccine_rate)
            u_vacc=vaccine_pulse(vaccine_rate(n),vaccine_start_day,vaccine_end_day,tf,dt,N);
            [S_vacc,I_vacc,R_vacc,~] = SIRV_model(alpha(k),beta,gamma,I0,R0,V0,tf-1,dt,u_vacc); 
            if any(S_vacc==0)
                vaccine_end_day=floor(find(S_vacc==0,1)*dt);
            end
            for i=2:length(I_vacc)
               newCases_vacc(i-1)=24*(I_vacc(i)+R_vacc(i)-I_vacc(i-1)-R_vacc(i-1));
               newDeaths_vacc(i-1)=24*mu*(R_vacc(i)-R_vacc(i-1));
            end
            if k==1   
                txt = [num2str(vaccine_rate_1000(n)),' vaccine per 1,000 (95\% efficacy)'];
                sz=sz+0.5;
                figure(1);plot(t(1:end-1),I_vacc*N,'LineWidth',sz,'DisplayName',txt);
                figure(2);plot(t(1:end-1),R_vacc*N,'LineWidth',sz,'DisplayName',txt);
                figure(3);plot(t(1:end-2),newCases_vacc*N,'LineWidth',sz,'DisplayName',txt);
                figure(4);plot(t(1:end-2),newDeaths_vacc*N,'LineWidth',sz,'DisplayName',txt);
            else
                txt = [num2str(vaccine_rate_1000(n)),' vaccine per 1,000 (60\% efficacy)'];
                sz=sz+0.5;
                figure(1);plot(t(1:end-1),I_vacc*N,'--','LineWidth',sz,'DisplayName',txt);
                figure(2);plot(t(1:end-1),R_vacc*N,'--','LineWidth',sz,'DisplayName',txt);
                figure(3);plot(t(1:end-2),newCases_vacc*N,'--','LineWidth',sz,'DisplayName',txt);
                figure(4);plot(t(1:end-2),newDeaths_vacc*N,'--','LineWidth',sz,'DisplayName',txt);
            end
    end
end

figure(1)
ylim([0 max(I)*N])
line=plot([endDate endDate],[0 max(I)*N],'k--','LineWidth',2);set(get(get(line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('Infected (active cases) $I(t)$','Interpreter','latex','fontsize',14)
ax=gca;
ax.XTick=tmonths;
xtickangle(90)
datetick('x','mmm-yy')
set(gca,'TickLabelInterpreter','latex')
title(population,'Interpreter','latex','fontsize',14)


figure(2)
ylim([0 max(R)*N])
line=plot([endDate endDate],[0 max(R)*N],'k--','LineWidth',2);set(get(get(line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('Removed (closed cases) $R(t)$','Interpreter','latex','fontsize',14)
ax=gca;
ax.XTick=tmonths;
xtickangle(90)
datetick('x','mmm-yy')
set(gca,'TickLabelInterpreter','latex')
title(population,'Interpreter','latex','fontsize',14)
legend show;
legend(gca, 'Location','best','Interpreter','latex','fontsize',14)


figure(3)
ylim([0 max(newCases)*N])
line=plot([endDate endDate],[0 max(newCases)*N],'k--','LineWidth',2);set(get(get(line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('Daily new infections (new cases)','Interpreter','latex','fontsize',14)
ax=gca;
ax.XTick=tmonths;
xtickangle(90)
datetick('x','mmm-yy')
set(gca,'TickLabelInterpreter','latex')
title(population,'Interpreter','latex','fontsize',14)

figure(4)
ylim([0 max(newDeaths)*N])
line=plot([endDate endDate],[0 max(newCases)*N],'k--','LineWidth',2);set(get(get(line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('Daily deaths','Interpreter','latex','fontsize',14)
ax=gca;
ax.XTick=tmonths;
xtickangle(90)
datetick('x','mmm-yy')
set(gca,'TickLabelInterpreter','latex')
title(population,'Interpreter','latex','fontsize',14)