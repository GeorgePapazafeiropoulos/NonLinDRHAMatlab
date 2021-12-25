%% SDOF system, dynamic analysis with direct integration (Chopra,2020)

%% Statement of the problem
% * This example comes from the 'Dynamics of Structures - Theory and
% Applications to Earthquake Engineering' by A. Chopra (edition 2020),
% Chapter 7, Figure 7.4.2.
% * Consider a linearly elastic system with weight w, natural vibration
% period Tn = 0.5 sec, and no damping. The excitation selected is the El
% Centro ground motion. Based on the time variation of the elastic
% resisting force fS, the peak value of this force fo is given by
% fo/mg=1.37.
% * Consider the response of an elastoplastic system having the same mass
% and initial stiffness as the linearly elastic system, with normalized
% strength fybar =0.125. The yield strength of this system is fy=0.125*fo,
% where fo=1.37*m*g (Fig. 7.4.1); therefore, fy=0.125*(1.37*m*g)=0.171*m*g.
% Calculate the following for the elastoplastic system: (a) deformation;
% (b) resisting force and acceleration; (c) time intervals of yielding; (d)
% force–deformation relation.
%
%% Initialization of structural input data
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=1;
%%
% Set the lumped mass equal to unity. This is permitted since the various
% response measures are normalized.
m=1;
%%
% Set the lateral stiffness.
Tn=0.5;
k_hi=m*(2*pi/Tn)^2;
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
inflvec=ones(nDOFs,1);
%% Load earthquake data
% Earthquake acceleration time history of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
D=load('elcentro.dat');
dt=D(2,1)-D(1,1);
xgtt=9.81*D(:,2);
t=dt*(0:(numel(xgtt)-1));
%%
% Post-yield stiffness
k_lo=0.001*k_hi;
%%
% Yield displacement
g=9.81;
fy=0.125*(1.37*m*g);
uy=fy/k_hi;
%%
% Set the critical damping ratio
ksi=0.0;
%%
% Initial displacement
u0=zeros(nDOFs,1);
%%
% Initial velocity
ut0=zeros(nDOFs,1);

%% Dynamic Response History Analysis (DRHA) with direct integration 
% Perform DRHA analysis for nonlinear structure
[u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy,ksi,u0,ut0);

%% Plot results
% Plot the following: (a) deformation; (b) resisting force and
% acceleration; (c) time intervals of yielding; (d) force–deformation
% relation.
FigHandle=figure('Name','Elastoplastic response','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 700, 750]);
%1st subplot
subplot(4,1,1)
plot(t,100*u(1,:),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t,zeros(size(t)),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold off
xlim([0,10])
ylim([-5,5])
set(gca,'XTick', (0:10));
set(gca,'YTick', (-5:2.5:5));
xlabel('Time (sec)','FontSize',10);
ylabel('u (cm)','FontSize',10);
title('(a)','FontSize',10)
grid on
%2nd subplot
subplot(4,1,2)
plot(t,Fs/(m*g),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t,fy/(m*g)*ones(size(t)),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
plot(t,-fy/(m*g)*ones(size(t)),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
text(9,0.22,num2str(fy/(m*g)))
text(9,-0.22,num2str(-fy/(m*g)))
hold off
xlim([0,10])
ylim([-0.3,0.3])
set(gca,'XTick', (0:10));
set(gca,'YTick', [-0.3,-0.15,0,0.15,0.3]);
xlabel('Time (sec)','FontSize',10);
ylabel('fs/mg','FontSize',10);
title('(b)','FontSize',10)
grid on
%3rd subplot
subplot(4,1,3)
plot(t,(abs(Fs-fy)<0.1)-(abs(Fs+fy)<0.1),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t,zeros(size(t)),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold off
xlim([0,10])
ylim([-1,1])
set(gca,'XTick', (0:10));
set(gca,'YTick', [-1,0,1]);
xlabel('Time (sec)','FontSize',10);
ylabel('Yield','FontSize',10);
title('(c)','FontSize',10)
grid on
%4th subplot
subplot(4,1,4)
plot(100*u(1,:),Fs/(m*g),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot([-5,5],zeros(size([-5,5])),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
plot(zeros(size([-0.3,0.3])),[-0.3,0.3],'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold off
xlim([-5,5])
ylim([-0.3,0.3])
set(gca,'XTick', (-5:5));
set(gca,'YTick', (-0.3:0.15:0.3));
xlabel('u (cm)','FontSize',10);
ylabel('fs/mg','FontSize',10);
title('(d)','FontSize',10)
grid on

%%
% Verify each subplot [(a)-(d)] with the Figure 7.4.2 of 'Dynamics of
% Structures - Theory and Applications to Earthquake Engineering' by A.
% Chopra.
%
% <<Chopra742.png>>
%

%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
