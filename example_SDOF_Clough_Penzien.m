%% Elastoplastic SDOF frame, dynamic analysis with direct integration (Clough & Penzien, 2003)

%% Statement of the problem
% * This example comes from the 'Dynamics of Structures' book by
% Clough & Penzien (edition 2003), Example E7-2.
% * The response of the elastoplastic SDOF frame shown in Fig. E7­3 to the
% loading history indicated is required.
% * Input parameters: k=5kips/in, m=0.1kips*sec^2/in, c=0.2kips*sec/in,
% uy=1.2in, dt=0.1
%
% <<CloughPenzienE73.png>>
%
%% Initialization of structural input data
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=1;
%%
% Set the lateral stiffness.
k_hi=5;
%%
% Set the lumped mass.
m=0.1;
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
inflvec=ones(nDOFs,1);
%% Load earthquake data
% Acceleration time history imposed at the base of the structure
dt=0.1;
xgtt=-[0;4;7;8;7;5;3;2;1;0.5;0;0;0]/m;
t=0:dt:1.2;
%%
% Post-yield stiffness
k_lo=0.001*k_hi;
%%
% Yield displacement (two values, a practically infinity value for the
% linear structure and a finite value for the nonlinear structure)
uy1=1e10;
uy2=1.2;
%%
% Set the critical damping ratio
ccrit=2*sqrt(k_hi*m);
c=0.2;
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
% Perform DRHA analysis for the yielding nonlinear structure
[u2,ut2,utt2,Fs2,Ey2,Es2,Ed2,jiter2] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy2,ksi,u0,ut0);

%% Plot results
% Plot the displacement response.
FigHandle=figure('Name','Displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 350]);
plot(t,u1,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t,u2,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[1 0 0],'markeredgecolor',[1 0 0])
plot(t,-xgtt'*m/k_hi,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 1],'markeredgecolor',[0 0 1])
hold off
grid on
xlim([0,1.2])
ylim([-1.5,3])
xlabel('Time (sec)','FontSize',10);
ylabel('U (in)','FontSize',10);
legend({'Elastic response','Elastoplastic response','Static displacement'},...
    'Location','Southwest')

%%
% Verify with Figure E7-4 of 'Dynamics of Structures', by Clough & Penzien
%
% <<CloughPenzienE74.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
