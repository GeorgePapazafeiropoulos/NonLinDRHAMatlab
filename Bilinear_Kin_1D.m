function [f,K,C,kstat,d] = Bilinear_Kin_1D(u,ut,k_hi,k_lo,uy,m,ksi,kstat,d)
%
% Bilinear elastoplastic hysteretic model with elastic viscous damping
%
% [F,K,C,KSTAT,D] = BILINEAR_KIN_1D(U,UT,K_HI,K_LO,UY,M,KSI,KSTAT,D)
%
% Description
%     Define the internal force vector, tangent stiffness matrix and
%     tangent damping matrix of a bilinear elastoplastic hysteretic
%     structure with elastic damping as a function of displacement and
%     velocity.
%     The MDOF structure modeled with this function consists of lumped
%     masses connected with stiffness and damping elements in series, along
%     a single dimension (i.e. one-dimensional MDOF system, hence the
%     function is named as '1D'). Each lumped mass has one degree of
%     freedom. The first degree of freedom is at the top of the structure
%     and the last at its fixed base. However, the last degree of freedom
%     is not included in the input arguments of the function, i.e. not
%     contained in ndof, as it is always fixed. The nonlinear stiffness
%     is virtually of the bilinear type, where an initial stiffness and a
%     post-yield stiffness are defined. The unloading or reloading curve of
%     this model are parallel to the initial loading curve, and a
%     hysteresis loop is created by continuously loading and unloading the
%     structure above its yield limit. This behavior can be viewed as
%     hardening of the kinematic type.
%     An appropriate reference for this function definition is Hughes,
%     Pister & Taylor (1979): "Implicit-explicit finite elements in
%     nonlinear transient analysis". This function should be defined in
%     accordance with equations (3.1), (3.2) and (3.3) of this paper. This
%     representation has as special cases nonlinear elasticity and a class
%     of nonlinear “rate-type” viscoelastic materials. Tangent stiffness
%     and tangent damping matrices are the "consistent" linearized
%     operators associated to f in the sense of [Hughes & Pister,
%     "Consistent linearization in mechanics of solids", Computers and
%     Structures, 8 (1978) 391-397].
%
% Input parameters
%     U [double(:ndof x 1)] is the absolute displacement (with respect to
%         ground). ndof is the number of degrees of freedom of the model.
%     UT [double(:ndof x 1)] is the absolute velocity (with respect to
%         ground).
%     K_HI [double(:ndof x 1)] is the initial stiffness of the system
%         before its first yield, i.e. the high stiffness. Give the
%         stiffness of each storey from top to bottom.
%     K_LO [double(:ndof x 1)] is the post-yield stiffness of the system,
%         i.e. the low stiffness. Give the stiffness of each storey from
%         top to bottom.
%     UY [double(:ndof x 1)] is the yield limit of the stiffness elements
%         of the structure. The element is considered to yield, if the
%         interstorey drift between degrees of freedom i and i+1 exceeds
%         UY(i). Give the yield limit of each storey from top to bottom.
%     M ([ndof x ndof]) is the lumped mass matrix of the system. Give
%         the lumped mass of each storey from top to bottom.
%     KSI (scalar): ratio of critical viscous damping of the system,
%         assumed to be unique for all damping elements of the structure.
%     KSTAT [double(:ndof x 1)] is the stiffness vector which takes into
%         account any plastic response of the structure. It is used to
%         record the status of the structure so that it is known before the
%         next application of this function at a next (time) step.
%         Initialize by setting KSTAT=K_HI.
%     D [double(:ndof x 1)] is the equilibrium displacement vector which
%         takes into account any plastic response of the structure. It is
%         used to record the status of the structure so that it is known
%         before the next application of this function at a next (time)
%         step. Initialize by setting D=zeros(ndof,1).
%
% Output parameters
%     F [double(:ndof x 1)]: internal force vector of the structure (sum of
%         forces due to stiffness and damping) at displacement u and
%         velocity UT
%     K ([ndof x ndof]): tangent stiffness matrix (nonlinear function
%         of displacement u and velocity UT). It is equivalent to the
%         derivative d(F)/d(U)
%     C ([ndof x ndof]): tangent damping matrix (nonlinear function
%         of displacement u and velocity UT). It is equivalent to the
%         derivative d(F)/d(U)
%     KSTAT [double(:ndof x 1)] is the stiffness vector which takes into
%         account any plastic response of the structure. It is used to
%         record the status of the structure so that it is known before the
%         next application of this function at a next (time) step.
%     D [double(:ndof x 1)] is the equilibrium displacement vector which
%         takes into account any plastic response of the structure. It is
%         used to record the status of the structure so that it is known
%         before the next application of this function at a next (time)
%         step.
%
% Examples
%     % 1) SDOF system
%     u=0:0.2:4;
%     ut=0.001*ones(1,numel(u));
%     u=[u,u(end:-1:1)];
%     ut=[ut,-ut];
%     u=[u,-u];
%     ut=[ut,ut(end:-1:1)];
%     u=[u u];
%     ut=[ut ut];
%     k_hi=1000;
%     k_lo=1;
%     uy=2;
%     M=1;
%     ksi=0.05;
%     k=k_hi;
%     d=0;
%     f=zeros(1,numel(u));
%     for i=1:numel(u)
%         [f(i),K,C,k,d] = Bilinear_Kin_1D(u(i),ut(i),k_hi,k_lo,uy,M,...
%             ksi,k,d);
%     end
%     figure()
%     plot(u,f)
%     % 2) System with 2 DOFs
%     u=0:0.2:4;
%     ut=0.001*ones(1,numel(u));
%     u=[u;u/2];
%     ut=[ut;ut/2];
%     u=[u,u(:,end:-1:1)];
%     ut=[ut,-ut];
%     u=[u,-u];
%     ut=[ut,ut(:,end:-1:1)];
%     u=[u u];
%     ut=[ut ut];
%     k_hi=[1000;2000];
%     k_lo=[1;1];
%     uy=[1;1];
%     M=eye(2);
%     ksi=0.05;
%     k=k_hi;
%     d=zeros(2,1);
%     f=zeros(2,size(u,2));
%     for i=1:size(u,2)
%         [f(:,i),K,C,k,d] = Bilinear_Kin_1D(u(:,i),ut(:,i),k_hi,k_lo,...
%             uy,M,ksi,k,d);
%     end
%     figure()
%     plot(u(1,:),f(1,:))
%     figure()
%     plot(u(2,:),f(2,:))
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

ndofs=size(u,1);
% Elastic tangent stiffness matrix
K=zeros(ndofs+1);
K=K-diag(k_hi,-1)-diag(k_hi,1);
K(1:end-1,1:end-1)=K(1:end-1,1:end-1)+diag(k_hi,0);
K(2:end,2:end)=K(2:end,2:end)+diag(k_hi,0);
K(end,:)=[];
K(:,end)=[];
% Mass matrix
M=diag(m);
% Elastic tangent damping matrix
[Eigvec,Eigval] = eig(K,M);
C=zeros(ndofs);
for i=1:ndofs
    c=2*ksi*sqrt(Eigval(i,i))*(((M*Eigvec(:,i))*Eigvec(:,i)')*M);
    C=C+c;
end
% loop over dofs
u=[u;0];
ut=[ut;0];
fK=zeros(ndofs,1);
fK_bak=zeros(ndofs+1,1);
for i=1:ndofs
    % force from stiffness (not damping) of the current storey
    fK(i)=kstat(i)*(u(i)-u(i+1)-d(i));
    % eq.(46) in ...
    fy=k_lo(i)*(u(i)-u(i+1))+(k_hi(i)-k_lo(i))*(uy(i).*sign(ut(i)-ut(i+1)));
    % check for yielding or load reversal
    if kstat(i)==k_hi(i) && ut(i)-ut(i+1)>0 && fK(i)>fy
        % check for yielding
        % the system has just exceeded its positive yield force level
        kstat(i)=k_lo(i);
        d(i)=(1-k_hi(i)/k_lo(i))*uy(i);
    elseif kstat(i)==k_hi(i) && ut(i)-ut(i+1)<0 && fK(i)<fy
        % check for yielding
        % the system has just exceeded its negative yield force level
        kstat(i)=k_lo(i);
        d(i)=(k_hi(i)/k_lo(i)-1)*uy(i);
    elseif kstat(i)==k_lo(i) && fK(i)*(ut(i)-ut(i+1))<0
        % check for load reversal
        % the system reloads from negative ultimate displacement or unloads
        % from positive ultimate displacement
        kstat(i)=k_hi(i);
        d(i)=(u(i)-u(i+1))-k_lo(i)/k_hi(i)*(u(i)-u(i+1)-d(i));
    end
    fK_bak(i+1)=kstat(i)*(u(i)-u(i+1)-d(i));
end
u(end)=[];
ut(end)=[];
% Update the elastic tangent stiffness matrix
K=zeros(ndofs+1);
K=K-diag(kstat,-1)-diag(kstat,1);
K(1:end-1,1:end-1)=K(1:end-1,1:end-1)+diag(kstat,0);
K(2:end,2:end)=K(2:end,2:end)+diag(kstat,0);
K(end,:)=[];
K(:,end)=[];
% internal force due to stiffness and damping
f = diff(fK_bak) + C*ut;

end