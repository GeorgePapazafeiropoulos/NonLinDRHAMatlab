%% Single-storey shear frame, dynamic analysis with direct integration (Rajasekaran, 2009)

%% Statement of the problem
% * This example comes from the 'Structural dynamics of earthquake
% engineering Theory and application using MATHEMATICA and MATLAB' book by
% S. Rajasekaran (edition 2009), Example 7.7.
% * A shear frame structure shown in Fig. 7.4a is subjected to time varying
% force shown in Fig. 7.21. Evaluate the elastic and elasto-plastic
% response of the structure. 
% * Input parameters: k=1897251N/m, m=43848kg, c=7767.7Ns/m, Rm=66825.6N,
% dt=0.001
% * NOTE: For the damping coefficient, the value of the Matlab code program
% 7.10 is valid, c=7767.7Ns/m, and not the one given in the problem
% statement c=34605.4Ns/m.
%
% <<Rajasekaran721.png>>
%
%% Initialization of structural input data
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=1;
%%
% Set the lateral stiffness.
k_hi=1897251;
%%
% Set the lumped mass.
m=43848;
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
inflvec=ones(nDOFs,1);
%% Load earthquake data
% Acceleration time history imposed at the base of the structure
dt=0.001;
xgtt=-[133651/m/0.5*(0:dt:0.5),133651/m/0.5*(1-((0.5+dt):dt:1)),...
    zeros(size((1+dt):dt:3))]';
t=0:dt:3;
%%
% Post-yield stiffness
k_lo=0.001*k_hi;
%%
% Yield displacement (two values, a practically infinity value for the
% linear structure and a finite value for the nonlinear structure)
uy1=1e10;
Rm=66825.6;
uy2=Rm/k_hi;
%%
% Set the critical damping ratio
ccrit=2*sqrt(k_hi*m);
c=7767.7;
ksi=c/ccrit;
%%
% Initial displacement
u0=zeros(nDOFs,1);
%%
% Initial velocity
ut0=zeros(nDOFs,1);

%% Dynamic Response History Analysis (DRHA) with direct integration
% Perform DRHA analysis for the non-yielding linear structure (uy1->inf)
[u1,ut1,utt1,Fs1,Ey1,Es1,Ed1,jiter1] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy1,ksi,u0,ut0);
%%
% Perform DRHA analysis for the yielding nonlinear structure (uy2=Rm/k_hi)
[u2,ut2,utt2,Fs2,Ey2,Es2,Ed2,jiter2] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy2,ksi,u0,ut0);

%% Plot results
% Plot the elastic displacement response.
FigHandle=figure('Name','Elastic displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 350]);
plot(t,u1,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,3])
ylim([-0.1,0.15])
xlabel('Time (sec)','FontSize',10);
ylabel('U (m)','FontSize',10);
title(['Elastic displacement, maxU=',num2str(max(abs(u1))),' m'],...
    'FontSize',10)
%%
% Verify with Figure 7.22(a) of 'Structural dynamics of earthquake
% engineering Theory and application using MATHEMATICA and MATLAB', by S.
% Rajasekaran
%
% <<Rajasekaran722a.png>>
%
%%
% Plot the elastic displacement response.
FigHandle=figure('Name','Elastoplastic displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 350]);
plot(t,u2,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,3])
ylim([0,0.4])
xlabel('Time (sec)','FontSize',10);
ylabel('U (m)','FontSize',10);
title(['Elastoplastic displacement, maxU=',num2str(max(abs(u2))),' m'],...
    'FontSize',10)
%%
% Verify with Figure 7.22(b) of 'Structural dynamics of earthquake
% engineering Theory and application using MATHEMATICA and MATLAB', by S.
% Rajasekaran
%
% <<Rajasekaran722b.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
