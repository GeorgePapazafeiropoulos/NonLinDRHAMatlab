%% Two-storey shear frame, dynamic analysis with direct integration (NPTEL)

%% Statement of the problem
% * This example comes from the Introduction to Earthquake Engineering -
% Web course of NPTEL (National Programme on Technology Enhanced Learning),
% Chapter 7, Non-linear Seismic Response of Structures, Example 7.3.
% * A two-story building is modeled as 2-DOF system and rigid floors as
% shown in the following figure. Determine the floor displacement responses
% due to El-Centro, 1940 earthquake ground motion. Take the inter-story
% stiffness, k =197.392 × 10^3 N/m, the floor mass, m = 2500 kg. The
% columns of the building are having elasto-plastic behavior with yield
% displacement of 0.05m.
%
% <<NPTEL79.png>>
%
%% Initialization of structural input data
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=2;
%%
% Set the lateral stiffness of each storey in N/m.
k_hi=197.392e3*[1;2];
%%
% Set the lumped mass at each floor in kg.
m=2500*[1;2];
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
k_lo=0.01*k_hi;
%%
% Yield displacement
uy=[0.05;0.05];
%%
% Set the critical damping ratio
ksi=0.02;
%%
% Initial displacement
u0=zeros(nDOFs,1);
%%
% Initial velocity
ut0=zeros(nDOFs,1);
%%
% Time integration algorithm
AlgID='U0-V0-Opt';
%%
% Minimum absolute value of the eigenvalues of the amplification matrix
rinf=1;
%%
% Maximum tolerance of convergence of the Full Newton Raphson method for
% numerical computation of acceleration
maxtol=0.01;
%%
% Maximum number of iterations per increment
jmax=200;
%%
% Infinitesimal acceleration for the calculation of the derivetive required
% for the convergence of the Newton-Raphson iteration
dak=eps;

%% Dynamic Response History Analysis (DRHA) with direct integration 
% Perform DRHA analysis for nonlinear structure
[u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy,ksi,u0,ut0,AlgID,rinf,maxtol,jmax,dak);

%% Plot results
% Plot the displacement time history of the second storey (1st storey
% according to the numbering convention of the NLDRHA_MDOF_1D function).
FigHandle=figure('Name','Roof displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 600, 250]);
plot(t,u(1,:),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-0.15,0.15])
xlabel('Time (sec)','FontSize',10);
ylabel('U1 (m)','FontSize',10);
title(['Second storey, maxU1=',num2str(max(abs(u(1,:)))),' m'],...
    'FontSize',10)
%%
% Verify with the top subplot of Figure 7.10 of Chapter 7, Non-linear
% Seismic Response of Structures, Example 7.3 of NPTEL
%
% <<NPTEL7101.png>>
%
%%
% Plot the force time history of the second storey (1st storey
% according to the numbering convention of the NLDRHA_MDOF_1D function).
FigHandle=figure('Name','Roof force','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 600, 250]);
plot(t,Fs(1,:)/1000,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-15,15])
xlabel('Time (sec)','FontSize',10);
ylabel('F1 (kN)','FontSize',10);
title(['Second storey, maxU1=',num2str(max(abs(Fs(1,:)/1000))),' kN'],...
    'FontSize',10)
%%
% Verify with the middle subplot of Figure 7.10 of Chapter 7, Non-linear
% Seismic Response of Structures, Example 7.3 of NPTEL
%
% <<NPTEL7102.png>>
%
%%
% Plot the force versus relative displacement (interstorey drift) between
% the second and first storeys (1st and 2nd storeys respectively according
% to the numbering convention of the NLDRHA_MDOF_1D function).
FigHandle=figure('Name','Roof force','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 600, 250]);
plot(u(1,:)-u(2,:),Fs(1,:)/1000,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([-0.1,0.1])
ylim([-15,15])
xlabel('U1-U2 (m)','FontSize',10);
ylabel('F1 (kN)','FontSize',10);
title(['Interstorey drift, max(U1-U2)=',num2str(max(u(1,:)-u(2,:))),' m'],...
    'FontSize',10)
%%
% Verify with the bottom subplot of Figure 7.10 of Chapter 7, Non-linear
% Seismic Response of Structures, Example 7.3 of NPTEL
%
% <<NPTEL7103.png>>
%
%%
% Plot the displacement time history of the first storey (2nd storey
% according to the numbering convention of the NLDRHA_MDOF_1D function).
FigHandle=figure('Name','Bottom floor displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 600, 250]);
plot(t,u(2,:),'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-0.1,0.1])
xlabel('Time (sec)','FontSize',10);
ylabel('U2 (m)','FontSize',10);
title(['Bottom floor displacement, maxU2=',num2str(max(u(2,:))),' m'],...
    'FontSize',10)
%%
% Verify with the top subplot of Figure 7.11 of Chapter 7, Non-linear
% Seismic Response of Structures, Example 7.3 of NPTEL
%
% <<NPTEL7111.png>>
%
%%
% Plot the force time history of the first storey (2nd storey
% according to the numbering convention of the NLDRHA_MDOF_1D function).
FigHandle=figure('Name','Bottom floor force','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 600, 250]);
plot(t,(Fs(1,:)+Fs(2,:))/1000,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-25,25])
xlabel('Time (sec)','FontSize',10);
ylabel('F2 (kN)','FontSize',10);
title(['Bottom floor force, maxF2=',num2str(max((Fs(1,:)+Fs(2,:))/1000)),' kN'],...
    'FontSize',10)
%%
% Verify with the middle subplot of Figure 7.11 of Chapter 7, Non-linear
% Seismic Response of Structures, Example 7.3 of NPTEL
%
% <<NPTEL7112.png>>
%
%%
% Plot the bottom floor force versus its displacement (i.e. 2nd storey
% according to the numbering convention of the NLDRHA_MDOF_1D function).
FigHandle=figure('Name','Bottom floor force','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 600, 250]);
plot(u(2,:),(Fs(1,:)+Fs(2,:))/1000,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([-0.1,0.1])
ylim([-25,25])
xlabel('U1 (m)','FontSize',10);
ylabel('F1 (kN)','FontSize',10);
title(['Displacement, maxU1=',num2str(max(u(2,:))),' m'],...
    'FontSize',10)
%%
% Verify with the bottom subplot of Figure 7.11 of Chapter 7, Non-linear
% Seismic Response of Structures, Example 7.3 of NPTEL
%
% <<NPTEL7113.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
