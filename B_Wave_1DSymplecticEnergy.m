clc; clear; close all;
% 1D Wave eqn, ASrinivasan, 19Dec21, Energy Evolution
% Use this for the MOLE library
% This code compares the performance for the TVD, PEFRL, RK4 & RRK 
% time integration schemes, with Mimetic 4th order spatial scheme. 

addpath('../mole_MATLAB')

NElem = 100; 
xmin = -5; xmax = 5; 
dh = (xmax-xmin)/NElem; 
CFL = 1.15; 

% init val
xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
xNod = [xmin:dh:xmax]';
muu = 0.; sg = 0.15;
u0 = 1/sqrt(sg*2*pi)*exp(-(xGrid-muu).^2/2/sg);
v0 = zeros(size(xGrid,1),1);
u0 = [u0;v0]; 
tEnd = 100; 

    
%%%%%%%%%%%
%%% ForestRuth
    Ordr = 4;      
    dt = 0.6*dh; 
    [uFRuth, tFRuth, eFRuth] = FRuth(NElem,dh,dt,tEnd,Ordr,u0);  
    energyFRuth = EnergyCalc(eFRuth);
    
%%% PEFRL
    Ordr = 4; 
    dt = CFL*dh; 
    [uPEFRL, tPEFRL, ePEFRL] = PEFRL(NElem,dh,dt,tEnd,Ordr,u0);  
    energyPEFRL = EnergyCalc(ePEFRL);

%%% RRK
    Ordr = 4; 
    [uRRK, tRRK,gamRRK, eRRK] = RRK(1,NElem,dh,dt,tEnd,Ordr,u0);  
    energyRRK = EnergyCalc(eRRK);
    uRRKEnd = InterpEnd(uRRK, tEnd+dt, tRRK);
    
    [uRK4, tRK4,gamRK4,eRK4] = RRK(0,NElem,dh,dt,tEnd,Ordr,u0);    
    energyRK4 = EnergyCalc(eRK4);

%%% TVD4, Shu & Osher
    Ordr = 4; 
    [uTVD, tTVD,eTVD] = TVD(NElem,dh,dt,tEnd,Ordr,u0);  
    energyTVD = EnergyCalc(eTVD);


f1 = figure;
plot(tFRuth(1:end-1,1), abs(energyFRuth),'-k','LineWidth', 2)
hold on
plot(tPEFRL(1:end-1,1), abs(energyPEFRL),'-*b','LineWidth', 1)
hold on
plot(tRK4(1:end-1,1), abs(energyRK4),'-.r','LineWidth', 2)
hold on
plot(tRRK(1:end-1,1), abs(energyRRK),'-+k','LineWidth', 2)
set(gca, 'YScale', 'log'); 
legend( 'Forest Ruth', 'PEFRL', 'RK4', 'RRK4', ...    
    'location', 'best'); 
set(gca,'FontSize',10)
grid on
xlabel('Time (s)')
ylabel('Energy Norm, abs(E^n - E^0)')
ylim([1E-16 1E2])
title({'1D Wave Eqn u_{tt} + \nabla\cdot(\nabla u) = 0, 4th order Mimetic', ...
    'E = ||V||^2 + ||\nabla U||^2'})


f2 = figure;
plot(xGrid, uFRuth(1:NElem+2,end), '-k', 'LineWidth', 2);
hold on
plot(xGrid, uPEFRL(1:NElem+2,end), '-*b', 'LineWidth', 1);
hold on
plot(xGrid, uRK4(1:NElem+2,end), '-.r', 'LineWidth', 2);
hold on
plot(xGrid, uRRKEnd(1:NElem+2,1), '-+k', 'LineWidth', 2);
hold on
legend('Forest Ruth', 'PEFRL', 'RK4', 'RRK4')   %, 'RK4'
set(gca,'FontSize',10)
ylim([-0.2 1.2])
xlabel('Spatial Domain x')
ylabel('Solution u at tEnd')
title('1D Wave Eqn u_{tt} + \nabla\cdot(\nabla u) = 0, 4th order Mimetic')
movegui(f2, 'south')

f3 = figure;
plot(xGrid, uRK4(1:NElem+2,end), '-.r', 'LineWidth', 2);
hold on
plot(xGrid, uRRKEnd(1:NElem+2,1), '-+k', 'LineWidth', 2);
hold on
legend('RK4', 'RRK4')   %, 'RK4'
set(gca,'FontSize',10)
ylim([-0.2 1.2])
xlabel('Spatial Domain x')
ylabel('Solution u at tEnd')
title('1D Wave Eqn u_{tt} + \nabla\cdot(\nabla u) = 0, 4th order Mimetic, CFL = 1.22')
movegui(f3, 'southwest')

f4 = figure;
plot(tRK4(1:end-1,1), abs(energyRK4),'-.r','LineWidth', 2)
hold on
plot(tRRK(1:end-1,1), abs(energyRRK),'-+k','LineWidth', 2)
set(gca, 'YScale', 'log'); 
legend( 'RK4', 'RRK4', ...    
    'location', 'best'); 
set(gca,'FontSize',10)
grid on
xlabel('Time (s)')
ylabel('Energy Norm, abs(E^n - E^0)')
ylim([1E-16 1E2])
title({'1D Wave Eqn u_{tt} + \nabla\cdot(\nabla u) = 0, 4th order Mimetic', ...
    'E = ||V||^2 + ||\nabla U||^2, CFL = 1.22'})





function [uFRuth, tFRuth, eFRuth] = FRuth(NElem,dh,dt,tEnd,Ordr,u0)
    
% Forest Ruth algorithm
     
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);
     
    uFRuth(:,1) = u0; % initial value
    y = u0; 
    tFRuth(1,1) = 0; % initial time    
    DU = G*y(1:NElem+2,1); 
    eFRuth(1,1) = dh*(norm(y(NElem+3:end,1),2)^2 + ...
        norm(DU,2)^2); 
    
    iCount = 2; 
    t = dt; 

    while t <= tEnd+dt
        UU = y(1:NElem+2,1);
        VV = y(NElem+3:end,1);
        Tht = 1/(2 - 2^(1/3));
        
        % Step 1
        U1 = UU + Tht * dt/2 * VV;
        V1 = VV + Tht * dt * D*G*U1;
        
        % Step 2        
        U2 = U1 + (1-Tht) * dt/2 * V1;
        V2 = V1 + (1-2*Tht) * dt * D*G*U2;

        % Step 3    
        U3 = U2 + (1-Tht) * dt/2 * V2;
        V3 = V2 + Tht * dt * D*G*U3;
        
        % Step 4    
        U4 = U3 + Tht * dt/2 * V3;   

        uNew = [U4;V3];                  
        
        uFRuth(:,iCount) = uNew;   
        y = uFRuth(:,iCount);
        tFRuth(iCount,1) = t;
        eFRuth(iCount,1) = dh*(norm(V3,2)^2 + norm(G*U4,2)^2); 
       
        iCount = iCount+1; 
        t = t + dt;           
        
    end   

end

function [uPERFL, tPERFL, ePERFL] = PEFRL(NElem,dh,dt,tEnd,Ordr,u0)
    
% Position Extended Forest Ruth algorithm
     
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);

     
    uPERFL(:,1) = u0; % initial value
    y = u0; 
    tPERFL(1,1) = 0; % initial time    
    DU = G*y(1:NElem+2,1); 
    ePERFL(1,1) = dh*(norm(y(NElem+3:end,1),2)^2 + ...
        norm(DU,2)^2); 
    
    iCount = 2; 
    t = dt; 

    while t <= tEnd+dt
        UU = y(1:NElem+2,1);
        VV = y(NElem+3:end,1);
        Zeta = 0.1786178958448091E+00;
        Lmb = -0.2123418310626054E+00;
        Xi = -0.6626458266981849E-1; 
        
        % Step 1
        U1 = UU + Zeta * dt * VV;
        V1 = VV + (1-2*Lmb) * dt/2 * D*G*U1;
        
        % Step 2        
        U2 = U1 + Xi * dt * V1;
        V2 = V1 + Lmb * dt * D*G*U2;

        % Step 3    
        U3 = U2 + (1-2*(Xi + Zeta)) * dt * V2;
        V3 = V2 + Lmb * dt * D*G*U3;
        
        % Step 4    
        U4 = U3 + Xi * dt * V3;   
        V4 = V3 + (1-2*Lmb)*dt/2 * D*G*U4;
        
        U5 = U4 + Zeta*dt*V4; 

        uNew = [U5;V4];                  
        
        uPERFL(:,iCount) = uNew;   
        y = uPERFL(:,iCount);
        tPERFL(iCount,1) = t;
        ePERFL(iCount,1) = dh*(norm(V3,2)^2 + norm(G*U4,2)^2); 
       
        iCount = iCount+1; 
        t = t + dt;           
        
    end   

end

function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tEnd,Ordr,u0)
    
% Relaxation RK4
     
    Z1 = zeros(NElem+2,NElem+2);
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);
    II = speye(NElem+2,NElem+2); 
    
    AMat = [Z1,II;
            D*G,Z1]; 
    
    AMat = sparse(AMat); 
               
    C(1,1) = 0;
    C(2,1) = 1/2;
    C(3,1) = 1/2;
    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;    
   

    B(1,1) = 1/6; 
    B(2,1) = 1/3;
    B(3,1) = 1/3;
    B(4,1) = 1/6;

     
    uRRK(:,1) = u0; % initial value
    y = u0; 
    tRRK(1,1) = 0; % initial time    
    gamRRK(1,1) = 1;
    E0 = y(1:NElem+2,1)'*G'*G*y(1:NElem+2,1) + ...
            y(NElem+3:end,1)'*y(NElem+3:end,1); 
    eRRK(1,1) = E0; 
    

    iCount = 2; 
    t = dt; 

    while t <= tEnd+dt
        z1 = y; 
        [k1] = AMat*z1; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = AMat*z2;

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = AMat*z3 ;
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = AMat*z4 ;        
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4;    
    
        switch RKFlag
            case 1
                U = y(1:NElem+2,1); V = y(NElem+3:end,1);
                dU = BKsum(1:NElem+2,1); dV = BKsum(NElem+3:end,1);
                
                AA = U'*G'*G*U + V'*V;
                BB = 2*(U'*G'*G*dU + V'*dV);
                CC = dU'*G'*G*dU + dV'*dV;
                
                gam = (E0 - AA - BB)/(dt* CC);
            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;                  
        
        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
 
       eRRK(iCount,1) = uNew(1:NElem+2,1)'*G'*G*uNew(1:NElem+2,1) + ...
            uNew(NElem+3:end,1)'*uNew(NElem+3:end,1);  
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end

function [uTVD, tTVD, eTVD] = TVD(NElem,dh,dt,tEnd,Ordr,u0)
    
% TVD, Shu & Osher
     
    Z1 = zeros(NElem+2,NElem+2);
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);
    II = speye(NElem+2,NElem+2); 
    
    AMat = [Z1,II;
            D*G,Z1]; 
    
    AMat = sparse(AMat); 
               
    C(1,1) = 0;
    C(2,1) = 1/2;
    C(3,1) = 1/2;
    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = -1/10;   A(3,2) = 3/5;
    A(4,1) = 0;   A(4,2) = 1/6;     A(4,3) = 5/6;    
   

    B(1,1) = 1/6; 
    B(2,1) = 7/18;
    B(3,1) = 5/18;
    B(4,1) = 1/6;

     
    uTVD(:,1) = u0; % initial value
    y = u0; 
    tTVD(1,1) = 0; % initial time    
    
    E0 = y(1:NElem+2,1)'*G'*G*y(1:NElem+2,1) + ...
            y(NElem+3:end,1)'*y(NElem+3:end,1); 
    eTVD(1,1) = E0; 
    

    iCount = 2; 
    t = dt; 

    while t <= tEnd+dt
        z1 = y; 
        [k1] = AMat*z1; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = AMat*z2;

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = AMat*z3 ;
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = AMat*z4 ;        
       
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4;    
    
        
        uNew = y + dt*BKsum;                  
        
        uTVD(:,iCount) = uNew;   
        y = uTVD(:,iCount);
        tTVD(iCount,1) = t;
        
 
       eTVD(iCount,1) = uNew(1:NElem+2,1)'*G'*G*uNew(1:NElem+2,1) + ...
            uNew(NElem+3:end,1)'*uNew(NElem+3:end,1);  
        
        iCount = iCount+1; 
        t = t + dt;           
        
    end   

end

function energyOut = EnergyCalc(uRRK)
    
   
    for i = 1:size(uRRK,1)-1
       energyOut(i,1) = (uRRK(i,1) - uRRK(1,1)); 
    end
    
end


function uEnd = InterpEnd(uWIC, tEnd, tWIC)
    % Ending value for uRRK, if tEnd > tWICRRK(end)
    if tWIC(end,1) ~= tEnd
       u2 = uWIC(:,end);  
       u1 = uWIC(:,end-1);
       t2 = tWIC(end,1);
       t1 = tWIC(end-1,1);
       uEnd = (u2-u1)*(tEnd-t1)/(t2-t1) + u1;    
    else 
       uEnd = uWIC(:,end);
    end

end
