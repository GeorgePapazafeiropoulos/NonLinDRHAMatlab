function [u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLDRHA_MDOF_1D(dt,xgtt,m,...
    inflvec,k_hi,k_lo,uy,ksi,u0,ut0,varargin)
%
% Non Linear Dynamic Response History Analysis of a bilinear kinematic MDOF
% shear building with elastic damping
% 
% [U,UT,UTT,FS,EY,ES,ED,JITER] = NLDRHA_MDOF_ENERGY(DT,XGTT,M,INFLVEC,...
%     K_HI,K_LO,UY,KSI,U0,UT0,ALGID,RINF,MAXTOL,JMAX,DAK)
%
% Description
%     General nonlinear direct time integration of equation of motion of a
%     bilinear elastoplastic hysteretic shear building, or any other MDOF
%     system with elastic damping, with lumped masses connected with
%     stiffness and damping elements in series (1D).
%     The General Single Step Single Solve (GSSSS) family of algorithms
%     published by X.Zhou & K.K.Tamma (2004) is employed for direct time
%     integration of the equations of motion of the Multiple Degree of
%     Freedom (MDOF) dynamic problem. Selection among 9 algorithms, all
%     designed according to the above journal article, can be made in this
%     routine. These algorithms encompass the scope of Linear Multi-Step
%     (LMS) methods and are limited by the Dahlquist barrier theorem
%     (Dahlquist,1963).
%
% Input parameters
%     DT [double(1 x 1)] is the size of the time step
%     XGTT [double(:NumSteps x 1)] is the acceleration time history which
%         is imposed at the lumped masses of teh structure, according to
%         the influence vector INFLVEC.
%     M [double(:ndof x 1)] is the lumped masses of the structure. Define
%         the lumped masses from the top to the bottom, excluding the fixed
%         dof at the base
%     INFLVEC [double(:ndof x 1)] is the influence vector. It determines
%         the dofs at which the acceleration prescribed in xgtt will be
%         imposed.
%     K_HI ([ndof x 1]): is the initial stiffness of the system before
%         its first yield, i.e. the high stiffness. Give the stiffness of
%         each storey from top to bottom.
%     K_LO ([ndof x 1]): is the post-yield stiffness of the system,
%         i.e. the low stiffness. Give the stiffness of each storey from
%         top to bottom.
%     UY ([ndof x 1]): is the yield limit of the stiffness elements of
%         the structure. The element is considered to yield, if the
%         interstorey drift between degrees of freedom i and i+1 exceeds
%         UY(i). Give the yield limit of each storey from top to bottom.
%     KSI [double(1 x 1)] is the ratio of critical viscous damping of the
%         system, assumed to be unique for all damping elements of the
%         structure.
%     ALGID [char(1 x :inf)] is the algorithm to be used for the time
%         integration. It can be one of the following strings for superior
%         optimally designed algorithms:
%             'generalized a-method': The generalized a-method (Chung &
%             Hulbert, 1993)
%             'HHT a-method': The Hilber-Hughes-Taylor method (Hilber,
%             Hughes & Taylor, 1977)
%             'WBZ': The Wood–Bossak–Zienkiewicz method (Wood, Bossak &
%             Zienkiewicz, 1980)
%             'U0-V0-Opt': Optimal numerical dissipation and dispersion
%             zero order displacement zero order velocity algorithm
%             'U0-V0-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) zero order displacement zero order
%             velocity algorithm
%             'U0-V0-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) zero order displacement zero order
%             velocity algorithm
%             'U0-V1-Opt': Optimal numerical dissipation and dispersion
%             zero order displacement first order velocity algorithm
%             'U0-V1-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) zero order displacement first order
%             velocity algorithm
%             'U0-V1-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) zero order displacement first order
%             velocity algorithm
%             'U1-V0-Opt': Optimal numerical dissipation and dispersion
%             first order displacement zero order velocity algorithm
%             'U1-V0-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) first order displacement zero order
%             velocity algorithm
%             'U1-V0-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) first order displacement zero order
%             velocity algorithm
%             'Newmark ACA': Newmark Average Constant Acceleration method
%             'Newmark LA': Newmark Linear Acceleration method
%             'Newmark BA': Newmark Backward Acceleration method
%             'Fox-Goodwin': Fox-Goodwin formula
%     U0 [double(:ndof x 1)] is the initial displacement. Give the initial
%         displacement of each storey from top to bottom.
%     UT0 [double(:ndof x 1)] is the initial velocity. Give the initial
%         velocity of each storey from top to bottom.
%     RINF [double(1 x 1)] is the minimum absolute value of the eigenvalues
%         of the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004).
%     MAXTOL [double(1 x 1)] is the maximum tolerance of convergence of the
%         Full Newton Raphson method for numerical computation of
%         acceleration.
%     JMAX [double(1 x 1)] is the maximum number of iterations per
%         increment. If JMAX=0 then iterations are not performed and the
%         MAXTOL parameter is not taken into account.
%     DAK [double(:ndof x 1)] is the infinitesimal acceleration for the
%         calculation of the derivetive required for the convergence of the
%         Newton-Raphson iteration.
%
% Output parameters
%     U [double(:ndof x :NumSteps)] is the time-history of displacement
%     UT [double(:ndof x :NumSteps)] is the time-history of velocity
%     UTT [double(:ndof x :NumSteps)] is the time-history of acceleration
%     FS [double(:ndof x :NumSteps)] is the time-history of the internal
%         force of the structure analysed.
%     EY [double(:ndof x :NumSteps)] is the time history of the sum of the
%         energy dissipated by yielding during each time step and the
%         recoverable strain energy of the system (incremental).
%         cumsum(EY(i,:))-ES(i,:) gives the time history of the
%         total energy dissipated by yielding from the start of the dynamic
%         analysis.
%     ES [double(:ndof x :NumSteps)] is the time-history of the recoverable
%         strain energy of the system (total and not incremental).
%     ED [double(:ndof x :NumSteps)] is the time-history of the energy
%         dissipated by viscoelastic damping during each time step
%         (incremental). cumsum(ED(i,:)) gives the time history of the
%         total energy dissipated at dof i from the start of the dynamic
%         analysis.
%     JITER [double(1 x :NumSteps)] is the iterations that are made per
%         increment until convergence
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin<10
    error('Not enough input arguments');
end
% Checks for dt
if ~isa(dt,'double')
    error('dt must be double');
end
if isempty(dt)
    error('dt is empty');
end
if any(~isfinite(dt))
    error('Inf or NaN values in dt');
end
if ~isempty(dt) && numel(dt)~=1
    error('dt must be scalar');
end
% Checks for xgtt
if ~isa(xgtt,'double')
    error('xgtt must be double');
end
if isempty(xgtt)
    error('xgtt is empty');
end
if any(~isfinite(xgtt))
    error('Inf or NaN values in xgtt');
end
if ~isempty(xgtt) && size(xgtt,2)~=1
    error('xgtt must be column vector');
end
% Checks for m
if ~isa(m,'double')
    error('m must be double');
end
if isempty(m)
    error('m is empty');
end
if any(~isfinite(m))
    error('Inf or NaN values in m');
end
if ~isempty(m) && size(m,2)~=1
    error('m must be column vector');
end
% Checks for inflvec
if ~isa(inflvec,'double')
    error('inflvec must be double');
end
if isempty(inflvec)
    error('inflvec is empty');
end
if any(~isfinite(inflvec))
    error('Inf or NaN values in inflvec');
end
if ~isempty(inflvec) && size(inflvec,2)~=1
    error('inflvec must be column vector');
end
% Checks for k_hi
if ~isa(k_hi,'double')
    error('k_hi must be double');
end
if isempty(k_hi)
    error('k_hi is empty');
end
if any(~isfinite(k_hi))
    error('Inf or NaN values in k_hi');
end
if ~isempty(k_hi) && size(k_hi,2)~=1
    error('k_hi must be column vector');
end
% Checks for k_lo
if ~isa(k_lo,'double')
    error('k_lo must be double');
end
if isempty(k_lo)
    error('k_lo is empty');
end
if any(~isfinite(k_lo))
    error('Inf or NaN values in k_lo');
end
if ~isempty(k_lo) && size(k_lo,2)~=1
    error('k_lo must be column vector');
end
% Checks for uy
if ~isa(uy,'double')
    error('uy must be double');
end
if isempty(uy)
    error('uy is empty');
end
if any(~isfinite(uy))
    error('Inf or NaN values in uy');
end
if ~isempty(uy) && size(uy,2)~=1
    error('uy must be column vector');
end
% Checks for ksi
if ~isa(ksi,'double')
    error('ksi must be double');
end
if isempty(ksi)
    error('ksi is empty');
end
if any(~isfinite(ksi))
    error('Inf or NaN values in ksi');
end
if ~isempty(ksi) && numel(ksi)~=1
    error('ksi must be scalar');
end
% Checks for u0
if ~isa(u0,'double')
    error('u0 must be double');
end
if isempty(u0)
    error('u0 is empty');
end
if any(~isfinite(u0))
    error('Inf or NaN values in u0');
end
if ~isempty(u0) && size(u0,2)~=1
    error('u0 must be column vector');
end
% Checks for ut0
if ~isa(ut0,'double')
    error('ut0 must be double');
end
if isempty(ut0)
    error('ut0 is empty');
end
if any(~isfinite(ut0))
    error('Inf or NaN values in ut0');
end
if ~isempty(ut0) && size(ut0,2)~=1
    error('ut0 must be column vector');
end
% Set default values for undefined input arguments
optargs = {'U0-V0-Opt' 1 0.01 200 eps};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[AlgID,rinf,maxtol,jmax,dak] = optargs{:};
% Checks for AlgID
if ~isa(AlgID,'char')
    error('AlgID must be of char type');
end
if isempty(AlgID)
    error('AlgID is empty');
end
if ~isempty(AlgID) && size(AlgID,1)~=1
    error('AlgID must be a row vector');
end
% Checks for rinf
if ~isa(rinf,'double')
    error('rinf must be double');
end
if isempty(rinf)
    error('rinf is empty');
end
if any(~isfinite(rinf))
    error('Inf or NaN values in rinf');
end
if ~isempty(rinf) && numel(rinf)~=1
    error('rinf must be scalar');
end
% Checks for maxtol
if ~isa(maxtol,'double')
    error('maxtol must be double');
end
if isempty(maxtol)
    error('maxtol is empty');
end
if any(~isfinite(maxtol))
    error('Inf or NaN values in maxtol');
end
if ~isempty(maxtol) && numel(maxtol)~=1
    error('maxtol must be scalar');
end
% Checks for jmax
if ~isa(jmax,'double')
    error('jmax must be double');
end
if isempty(jmax)
    error('jmax is empty');
end
if any(~isfinite(jmax))
    error('Inf or NaN values in jmax');
end
if ~isempty(jmax) && numel(jmax)~=1
    error('jmax must be scalar');
end
% Checks for dak
if ~isa(dak,'double')
    error('dak must be double');
end
if isempty(dak)
    error('dak is empty');
end
if any(~isfinite(dak))
    error('Inf or NaN values in dak');
end
if ~isempty(dak) && numel(dak)~=1
    error('dak must be scalar');
end

% Set integration constants
if all(size(AlgID)==[1,14])
    % define integration constants explicitly
    w1=AlgID(1);
    w2=AlgID(2);
    w3=AlgID(3);
    W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
    % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
    % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
    W1L1=AlgID(4);
    W2L2=AlgID(5);
    W3L3=AlgID(6);
    W1L4=AlgID(7);
    W2L5=AlgID(8);
    W1L6=AlgID(9);
    l1=AlgID(10);
    l2=AlgID(11);
    l3=AlgID(12);
    l4=AlgID(13);
    l5=AlgID(14);
else
    switch AlgID
        case 'U0-V0-Opt'
            % zero-order displacement & velocity overshooting behavior and
            % optimal numerical dissipation and dispersion
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % mid-point rule a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf); % suggested
            w2=15*(3-4*rinf)/(1-4*rinf); % suggested
            w3=-35*(1-rinf)/(1-4*rinf); % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/2/(1+rinf)^2;
            W1L4=1/(1+rinf);
            W2L5=1/(1+rinf)^2; % suggested
            W1L6=(3-rinf)/2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V0-CA'
            % zero-order displacement & velocity overshooting behavior and
            % continuous acceleration
            % rinf must belong to [1/3 1]
            if rinf<1/3
                rinf=1/3;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/3');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-5*rinf)/(3-7*rinf); % suggested
            w2=15*(1-13*rinf)/(3-7*rinf); % suggested
            w3=140*rinf/(3-7*rinf); % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(1+3*rinf)/2/(1+rinf);
            W2L2=(1+3*rinf)/4/(1+rinf);
            W3L3=(1+3*rinf)/4/(1+rinf)^2;
            W1L4=(1+3*rinf)/2/(1+rinf);
            W2L5=(1+3*rinf)/2/(1+rinf)^2; % suggested
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V0-DA'
            % zero-order displacement & velocity overshooting behavior and
            % discontinuous acceleration
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15; % suggested
            w2=45; % suggested
            w3=-35; % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/2/(1+rinf);
            W1L4=1;
            W2L5=1/(1+rinf); % suggested
            W1L6=(3+rinf)/2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V1-Opt'
            % zero-order displacement & first-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % This is the generalized a-method (Chung & Hulbert, 1993)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf);
            w2=15*(3-4*rinf)/(1-4*rinf);
            w3=-35*(1-rinf)/(1-4*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/(1+rinf)^3;
            W1L4=1/(1+rinf);
            W2L5=(3-rinf)/2/(1+rinf)^2;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'generalized a-method'
            % zero-order displacement & first-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % This is the generalized a-method (Chung & Hulbert, 1993)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf);
            w2=15*(3-4*rinf)/(1-4*rinf);
            w3=-35*(1-rinf)/(1-4*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/(1+rinf)^3;
            W1L4=1/(1+rinf);
            W2L5=(3-rinf)/2/(1+rinf)^2;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U0-V1-CA'
            % zero-order displacement & first-order velocity overshooting
            % behavior and continuous acceleration
            % This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
            % Taylor, 1977)
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(2-3*rinf);
            w2=15*(2-5*rinf)/(2-3*rinf);
            w3=-35*(1-3*rinf)/2/(2-3*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=2*rinf/(1+rinf);
            W2L2=rinf/(1+rinf);
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=2*rinf/(1+rinf);
            W2L5=rinf*(3-rinf)/(1+rinf)^2;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'HHT a-method'
            % zero-order displacement & first-order velocity overshooting
            % behavior and continuous acceleration
            % This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
            % Taylor, 1977)
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(2-3*rinf);
            w2=15*(2-5*rinf)/(2-3*rinf);
            w3=-35*(1-3*rinf)/2/(2-3*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=2*rinf/(1+rinf);
            W2L2=rinf/(1+rinf);
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=2*rinf/(1+rinf);
            W2L5=rinf*(3-rinf)/(1+rinf)^2;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U0-V1-DA'
            % zero-order displacement & first-order velocity overshooting
            % behavior and discontinuous acceleration
            % This is the Wood–Bossak–Zienkiewicz method (Wood, Bossak &
            % Zienkiewicz, 1980)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/(1+rinf)^2;
            W1L4=1;
            W2L5=(3-rinf)/2/(1+rinf);
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'WBZ'
            % zero-order displacement & first-order velocity overshooting
            % behavior and discontinuous acceleration
            % This is the Wood–Bossak–Zienkiewicz method (Wood, Bossak &
            % Zienkiewicz, 1980)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/(1+rinf)^2;
            W1L4=1;
            W2L5=(3-rinf)/2/(1+rinf);
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U1-V0-Opt'
            % first-order displacement & zero-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % mid-point rule a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-30*(3-8*rinf+6*rinf^2)/(9-22*rinf+19*rinf^2);
            w2=15*(25-74*rinf+53*rinf^2)/2/(9-22*rinf+19*rinf^2);
            w3=-35*(3-10*rinf+7*rinf^2)/(9-22*rinf+19*rinf^2);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(3-rinf)/2/(1+rinf);
            W2L2=1/(1+rinf)^2;
            W3L3=1/(1+rinf)^3;
            W1L4=(3-rinf)/2/(1+rinf);
            W2L5=2/(1+rinf)^3;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U1-V0-CA'
            % first-order displacement & zero-order velocity overshooting
            % behavior and continuous acceleration
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-60*(2-8*rinf+7*rinf^2)/(11-48*rinf+41*rinf^2);
            w2=15*(37-140*rinf+127*rinf^2)/2/(11-48*rinf+41*rinf^2);
            w3=-35*(5-18*rinf+17*rinf^2)/(11-48*rinf+41*rinf^2);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(1+3*rinf)/2/(1+rinf);
            W2L2=2*rinf/(1+rinf)^2;
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=(1+3*rinf)/2/(1+rinf);
            W2L5=4*rinf/(1+rinf)^3;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U1-V0-DA'
            % first-order displacement & zero-order velocity overshooting behavior
            % and discontinuous acceleration
            % This is the Newmark average acceleration a-form algorithm
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-30*(3-4*rinf)/(9-11*rinf);
            w2=15*(25-37*rinf)/2/(9-11*rinf);
            w3=-35*(3-5*rinf)/(9-11*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(3+rinf)/2/(1+rinf);
            W2L2=1/(1+rinf);
            W3L3=1/(1+rinf)^2;
            W1L4=(3+rinf)/2/(1+rinf);
            W2L5=2/(1+rinf)^2;
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'Newmark ACA'
            % Newmark Average Constant Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=0.25;
            W3L3=0.25;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=0.25;
            l4=1;
            l5=0.5;
        case 'Newmark LA'
            % Newmark Linear Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/6;
            W3L3=1/6;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=1/6;
            l4=1;
            l5=0.5;
        case 'Newmark BA'
            % Newmark Backward Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=0.5;
            W3L3=0.5;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=0.5;
            l4=1;
            l5=0.5;
        case 'Fox-Goodwin'
            % Fox-Goodwin formula
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/12;
            W3L3=1/12;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=1/12;
            l4=1;
            l5=0.5;
        otherwise
            error('No appropriate algorithm specified.');
    end
end

%% Calculation
% number of analysis increments
NumSteps=length(xgtt);
% number of dofs
ndof=length(m);
% Initialize output
u=zeros(ndof,NumSteps);
ut=zeros(ndof,NumSteps);
utt=zeros(ndof,NumSteps);
Fs=zeros(ndof,NumSteps);
Ey=zeros(ndof,NumSteps);
Es=zeros(ndof,NumSteps);
Ed=zeros(ndof,NumSteps);
jiter=zeros(1,NumSteps);
% set initial values of displacement, velocity, acceleration (u0,ut0 and
% utt0 respectively) at n=0.
u(:,1)=u0;
ut(:,1)=ut0;
% construct lumped mass matrix
M=diag(m,0);
% calculation for first increment
k=k_hi;
d=zeros(ndof,1);
[FKC0,K0,C0,k,d] = Bilinear_Kin_1D(u0,ut0,k_hi,k_lo,uy,m,ksi,k,d);
utt0=-xgtt(1)*inflvec-M\FKC0;
utt(:,1)=utt0;
Fs(:,1)=FKC0;
% initial assignments
FKCn=FKC0;
Kn=K0;
Cn=C0;
un=u0;
utn=ut0;
uttn=utt0;
% integration increments n
for n=1:NumSteps-1
    % effective force
    Feffn1k=-FKCn...
        -Kn*(W1L1*dt*utn+W2L2*dt^2*uttn)...
        -W1L4*dt*Cn*uttn...
        +M*(((1-W1)*xgtt(n)+W1*xgtt(n+1))*inflvec-uttn);
    % effective mass
    Meffn=W1L6*M+W2L5*dt*Cn+W3L3*dt^2*Kn;
    % initial estimate of da
    dan=Meffn\Feffn1k;
    
    % start iteration number k
    j=1;
    % set initial quotient of variation of da equal to maxtol
    quda=maxtol;
    % full Newton-Raphson iterations k
    while max(abs(quda))>=maxtol && j<=jmax
        % iteration k+1 of increment n+1
        
        % displacement, velocity, acceleration, internal force, stiffness
        % and damping for uttn+dan
        % calculate the residual Rn1k at uttn+dan
        % update kinematic quantities
        un1k=un+l1*utn*dt+l2*uttn*dt^2+l3*dan*dt^2;
        utn1k=utn+l4*uttn*dt+l5*dan*dt;
        uttn1k=uttn+dan;
        % force due to stiffness and damping
        [FKCn1k,Kn1k,Cn1k,~,~]=Bilinear_Kin_1D(un1k,utn1k,k_hi,k_lo,uy,m,ksi,k,d);
        % effective force
        Feffn1k=-FKCn1k...
            -Kn1k*(W1L1*dt*utn1k+W2L2*dt^2*uttn1k)...
            -Cn1k*W1L4*dt*uttn1k...
            +M*((1-W1)*xgtt(n)+W1*xgtt(n+1)-uttn1k);
        % effective mass
        Meffn1k=Kn1k*W3L3*dt^2+Cn1k*W2L5*dt+M*W1L6;
        % residual
        Rn1k=Feffn1k-Meffn1k*dan;

        % displacement, velocity, acceleration, internal force, stiffness
        % and damping for uttn+(dan+dak)
        % calculate the derivative at uttn+dan as:
        % dR/da=(dRn1k-Rn1k)/(uttn+(dan+dak)-(uttn+dan))=(dRn1k-Rn1k)/dak
        % update kinematic quantities
        dun1k=un+l1*utn*dt+l2*uttn*dt^2+l3*(dan+dak)*dt^2;
        dutn1k=utn+l4*uttn*dt+l5*(dan+dak)*dt;
        duttn1k=uttn+(dan+dak);

        % force due to stiffness and damping
        [dFKCn1k,dKn1k,dCn1k,~,~]=Bilinear_Kin_1D(dun1k,dutn1k,k_hi,k_lo,uy,m,ksi,k,d);
        % effective force
        dFeffn1k=-dFKCn1k...
            -dKn1k*(W1L1*dt*dutn1k+W2L2*dt^2*duttn1k)...
            -dCn1k*W1L4*dt*duttn1k...
            +M*((1-W1)*xgtt(n)+W1*xgtt(n+1)-duttn1k);
        % effective mass
        dMeffn1k=dKn1k*W3L3*dt^2+dCn1k*W2L5*dt+M*W1L6;
        % residual
        dRn1k=dFeffn1k-dMeffn1k*duttn1k;
        
        % Full Newton-Raphson update:
        % da_new=da-Rn1k/(dR/da)=da*(1-Rn1k/(dRn1k/dak)/da)
        % (to be checked for while loop termination)
        quda=(Rn1k./(dRn1k-Rn1k).*dak)./dan;
        % test if derivative becomes zero
        a=isinf(quda);
        if any(a)
            break
            %quda=zeros(size(quda));
        end
        % update da
        dan=(1-quda).*dan;
        % update iteration number
        j=j+1;
    end
    
    % displacement and its derivatives after iteration k+1 of increment
    % n+1
    un1k=un+l1*utn*dt+l2*uttn*dt^2+l3*dan*dt^2;
    utn1k=utn+l4*uttn*dt+l5*dan*dt;
    uttn1k=uttn+dan;
    % internal force, stiffness and damping after iteration k+1 of
    % increment n+1
    [FKCn1k,Kn1k,Cn1k,k,d] = Bilinear_Kin_1D(un1k,utn1k,k_hi,k_lo,uy,m,ksi,k,d);
    
    % assignments to output parameters
    u(:,n+1)=-un1k;
    ut(:,n+1)=-utn1k;
    utt(:,n+1)=-uttn1k;
    Fs(:,n+1)=-FKCn1k;
    Ey(:,n+1)=-(cumsum(FKCn1k-Cn1k*utn1k)+cumsum(FKCn-Cn*utn))/2.*diff([un1k-un;0]);
    Es(:,n+1)=cumsum(FKCn1k-Cn1k*utn1k).^2./k_hi/2;
    Ed(:,n+1)=-(cumsum(Cn1k*utn1k)+cumsum(Cn*utn))/2.*diff([un1k-un;0]);
    jiter(n+1)=j-1;
    % assignments for next increment
    FKCn=FKCn1k;
    Kn=Kn1k;
    Cn=Cn1k;
    un=un1k;
    utn=utn1k;
    uttn=uttn1k;
end

end
