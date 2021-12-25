%% Energy dissipation of SDOF system (Chopra,2020)

%% Statement of the problem
% * This example comes from the 'Dynamics of Structures - Theory and
% Applications to Earthquake Engineering' by A. Chopra (edition 2020),
% Chapter 7, Figure 7.9.1.
% * Calculate the variation of the various energy quantities with time for
% two SDOF systems subjected to the El Centro ground motion. The results
% presented are for a linearly elastic system with natural period Tn=0.5
% sec and damping ratio ksi=0.05, and for an elastoplastic system with the
% same properties in the elastic range and normalized strength fybar=0.25.
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
% Yield displacement for elastic response
uy1=1e10;
%%
% Set the critical damping ratio
ksi=0.05;
%%
% Initial displacement
u0=zeros(nDOFs,1);
%%
% Initial velocity
ut0=zeros(nDOFs,1);
%% Dynamic Response History Analysis (DRHA) with direct integration 
% Perform DRHA analysis for linear structure
[u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy1,ksi,u0,ut0);

%% Plot results
% Plot the damping energy and strain energy of the linearly elastic SDOF
% system.
FigHandle=figure('Name','Elastic response','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 700, 750]);
plot(t',cumsum(Ed)*1e4,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t',cumsum(Ed)*1e4+Es*1e4,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold off
xlim([0,30])
ylim([0,8000])
set(gca,'XTick', (0:10:30));
set(gca,'YTick', (0:2000:8000));
xlabel('Time (sec)','FontSize',10);
ylabel('Energy/umit mass (cm/s)^2','FontSize',10);
title('(a)','FontSize',10)
grid on
%%
% Verify with the Figure 7.9.1(a) of 'Dynamics of Structures - Theory and
% Applications to Earthquake Engineering' by A. Chopra.
%
% <<Chopra791a.png>>
%
%% Dynamic Response History Analysis (DRHA) with direct integration 
% normalized yield strength
fybar=0.25;
%%
% Yield displacement for nonlinear response
uy2=fybar*max(abs(u));
%%
% Perform DRHA analysis for nonlinear structure
[u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy2,ksi,u0,ut0);

%% Plot results
% Plot the damping energy and strain energy of the linearly elastic SDOF
% system.
FigHandle=figure('Name','Elastoplastic response','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 700, 750]);
plot(t',cumsum(Ed)*1e4,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t',cumsum(Ed)*1e4+cumsum(Ey)*1e4-Es*1e4,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
plot(t',cumsum(Ed)*1e4+cumsum(Ey)*1e4,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold off
xlim([0,30])
ylim([0,8000])
set(gca,'XTick', (0:10:30));
set(gca,'YTick', (0:2000:8000));
xlabel('Time (sec)','FontSize',10);
ylabel('Energy/umit mass (cm/s)^2','FontSize',10);
title('(b)','FontSize',10)
grid on
%%
% Verify with the Figure 7.9.1(b) of 'Dynamics of Structures - Theory and
% Applications to Earthquake Engineering' by A. Chopra.
%
% <<Chopra791b.png>>
%

%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
