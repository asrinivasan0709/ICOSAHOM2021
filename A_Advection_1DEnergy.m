clc; clear; close all;
% 1D advection eqn, to check energy conservation
% ASrinivasan, 24DEC21
% Use this version for MOLE. 

addpath('../mole_MATLAB')

NElem = 100;
xmin = -5; xmax = 5; 

dh = (xmax-xmin)/NElem; 
tEnd = 50;  
CFL = 1;
dt =  CFL*dh; 

% init val
xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
xNod = [xmin:dh:xmax]';
u0 = zeros(size(xGrid,1),1);
for i = 1:size(xGrid,1)
    if abs(xGrid(i,1)) <= 0.5
       u0(i,1) = cos(pi*xGrid(i,1)).^2;  
    end
end

%%%%%%%%%%%
%%% RRK
    Ordr = 4; 
    [uRRK, tRRK,gamRRK, eRRK] = RRK(1,NElem,dh,dt,tEnd,Ordr,u0);    
    energyRRK = EnergyCalc(eRRK,dh);
    
    [uRK4, tRK4,gamRK4, eRK4] = RRK(0,NElem,dh,dt,tEnd,Ordr,u0);    
    energyRK4 = EnergyCalc(eRK4,dh);

    dt = 0.5*dh; 
    Ordr = 6; 
    [uRRK6, tRRK6,gamRRK6, eRRK6] = RRK6(1,NElem,dh,dt,tEnd,Ordr,u0);    
    energyRRK6 = EnergyCalc(eRRK6,dh);
    
    [uRK6, tRK6,gamRK6, eRK6] = RRK6(0,NElem,dh,dt,tEnd,Ordr,u0);    
    energyRK6 = EnergyCalc(eRK6,dh);

% Heun - 3rd order time, 4th space
    dt = CFL*dh; 
    Ordr = 4; 
    [uRHn, tRHn,gamRHn, eRHn] = RHeun(1,NElem,dh,dt,tEnd,Ordr,u0);    
    energyRHn = EnergyCalc(eRHn,dh);
    
    [uHn, tHn,gamHn, eHn] = RHeun(0,NElem,dh,dt,tEnd,Ordr,u0);    
    energyHn = EnergyCalc(eHn,dh);
    
    


f1 = figure;
plot(tRRK6, abs(energyRRK6),'-ok','LineWidth', 2)
hold on
plot(tRK6, abs(energyRK6), '-.k','LineWidth', 3)
hold on
plot(tRRK, abs(energyRRK),'-sb','LineWidth', 2)
hold on
plot(tRK4, abs(energyRK4), '-.b','LineWidth', 3)
hold on
plot(tRHn, abs(energyRHn),'-+r','LineWidth', 2)
hold on
plot(tHn, abs(energyHn), '-.r','LineWidth', 3)
set(gca, 'YScale', 'log'); 
legend('M6-RRK6', 'M6-RK6', 'M4-RRK4', 'M4-RK4', ...
    'M4-RRK3', 'M4-RK3', ...
    'location', 'best'); 
set(gca,'FontSize',10)
ylim([1E-15 1E5])
xlim([0 tEnd])
xlabel('Time (s)')
ylabel('Energy Norm, abs(||u^n||^2 - ||u^0||^2)')
title({'1D Advection Eqn u_t + \nabla\cdot(u) = 0', 'Mimetic RK Scheme'})


f2 = figure;
plot(xGrid,uRK4(:,end), '-.b', 'linewidth', 2);
hold on
plot(xGrid,uRRK(:,end), 'k', 'linewidth', 2);
hold on
plot(xGrid, u0, '--k', 'linewidth', 1);
legend('MIM4-RK4', 'MIM4-RRK4', 'Exact');
set(gca,'FontSize',10)
xlabel('Spatial Domain x')
ylabel('Solution u')
ylim([-0.2 1.2]);
title({'1D Advection Eqn u_t + \nabla\cdot(u) = 0', '4th Mimetic & 4th order RK Scheme'})


f3 = figure;
plot(xGrid,uHn(:,end), '-.b', 'linewidth', 2);
hold on
plot(xGrid,uRHn(:,end), 'k', 'linewidth', 2);
hold on
plot(xGrid, u0, '--k', 'linewidth', 1);
legend('MIM4-RK3', 'MIM4-RRK3', 'Exact');
set(gca,'FontSize',10)
xlabel('Spatial Domain x')
ylabel('Solution u')
ylim([-0.2 1.2]);
title({'1D Advection Eqn u_t + \nabla\cdot(u) = 0', '4th Mimetic & 3rd order RK Scheme'})



f4 = figure;
iOut = 1000;
plot(xGrid,uRK6(:,end), '-.b', 'linewidth', 1);
hold on
plot(xGrid,uRRK6(:,end), 'k', 'linewidth', 2);
hold on
plot(xGrid, u0, '--k', 'linewidth', 1);
legend('MIM6-RK6', 'MIM6-RRK6', 'Exact');
set(gca,'FontSize',10)
xlabel('Spatial Domain x')
ylabel('Solution u')
ylim([-0.2 1.2]);
title({'1D Advection Eqn u_t + \nabla\cdot(u) = 0', '6th Mimetic & 6th order RK Scheme'})

f5 = figure;
plot(tRK6, abs(energyRK6), '-.k','LineWidth', 3)
hold on
plot(tRK4, abs(energyRK4), '-.b','LineWidth', 3)
hold on
plot(tHn, abs(energyHn), '-.r','LineWidth', 3)
set(gca, 'YScale', 'log'); 
legend('MIM6-RK6',  'MIM4-RK4', ...
    'MIM4-RK3', ...
    'location', 'best'); 
set(gca,'FontSize',10)
xlim([0 tEnd])
xlabel('Time (s)')
ylabel('Energy Norm, abs(||u^n||^2 - ||u^0||^2)')
title({'1D Advection Eqn u_t + \nabla\cdot(u) = 0', 'Mimetic RK Scheme'})




f6 = figure;
plot(tRRK6, (energyRRK6),'-ok','LineWidth', 2)
hold on
plot(tRRK, (energyRRK),'-sb','LineWidth', 2)
hold on
plot(tRHn, (energyRHn),'-+r','LineWidth', 2)
legend('MIM6-RRK6',  'MIM4-RRK4', ...
    'MIM4-RRK3', ...
    'location', 'best'); 
set(gca,'FontSize',10)
xlim([0 tEnd])
xlabel('Time (s)')
ylabel('Energy Norm, (||u^n||^2 - ||u^0||^2)')
title({'1D Advection Eqn u_t + \nabla\cdot(u) = 0', 'Mimetic Relaxation RK Scheme'})


    
function [uRHeun, tRHeun, gamRHeun, eRHn] = RHeun(RKFlag,NElem,dh,dt,tEnd,Ordr,u0)
% Heun's 3rd order method
    
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);
    IntD = interpDMat(Ordr, NElem);     
    
    % Convert IntD to periodic
    
    IntD(2,:) = zeros(1,NElem+2); 
    IntD(end-1,:) = zeros(1,NElem+2); 
    
    IntD(2,2:4) = IntD(3,3:5); IntD(2,end-1) = IntD(3,2); 
    IntD(end-1,:) = fliplr(IntD(2,:));    
    
    % End conversion IntD to periodic

                
    C(1,1) = 0;
    C(2,1) = 1/3;     
    C(3,1) = 2/3;

    A(2,1) = 1/3;  
    A(3,1) = 0; A(3,2) = 2/3;

    B(1,1) = 1/4; 
    B(2,1) = 0;
    B(3,1) = 3/4; 

     
    uRHeun(:,1) = u0; % initial value
    y = u0; 
    tRHeun(1,1) = 0; % initial time    
    gamRHeun(1,1) = 1;
    E0 = y'*G'*G*y; 
    eRHn(1,1) = E0; 
     
    
    iCount = 2; 
    t = dt; 
    
    while t <= tEnd+dt
        z1 = y;          
        [k1] = - D*IntD*z1; %
        % Convert to periodic boundaries
        k1(1,1) = k1(2:3,1)'*IntD(3,4:5)' + ...
                k1(end-2:end-1)'*IntD(3,2:3)'; 
        k1(end,end) = k1(1,1);     

        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;        
        [k2] = - D*IntD*z2 ; %
        % Convert to periodic boundaries
        k2(1,1) = k2(2:3,1)'*IntD(3,4:5)' + ...
                k2(end-2:end-1)'*IntD(3,2:3)'; 
        k2(end,end) = k2(1,1);     


        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2);             
        [k3] = - D*IntD*z3 ; %
        % Convert to periodic boundaries
        k3(1,1) = k3(2:3,1)'*IntD(3,4:5)' + ...
                k3(end-2:end-1)'*IntD(3,2:3)'; 
        k3(end,end) = k3(1,1);     

              
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3;
              
    
        switch RKFlag
            case 1                
                U = y; 
                dU = BKsum; 
                
                AA = U'*G'*G*U ;
                BB = 2*(U'*G'*G*dU);
                CC = dU'*G'*G*dU;
                
                gam = (E0 - AA - BB)/(dt* CC); 

            case 0
                gam = 1;
        end
                    
        uNew = y  + gam*dt*BKsum;       
        
        
        uRHeun(:,iCount) = uNew;   
        y = uRHeun(:,iCount);
        tRHeun(iCount,1) = t;
        gamRHeun(iCount,1) = gam; 
        eRHn(iCount,1) = uNew'*G'*G*uNew; 
        
        iCount = iCount+1; 
        t = t + gam*dt;      
        
        
    end
    

end

function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tEnd,Ordr,u0)
    
% Relaxation RK4    
    
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);
    IntD = interpDMat(Ordr, NElem);     
    
%     % Convert IntD to periodic
%     
    IntD(2,:) = zeros(1,NElem+2); 
    IntD(end-1,:) = zeros(1,NElem+2); 
    
    IntD(2,2:4) = IntD(3,3:5); IntD(2,end-1) = IntD(3,2); 
    IntD(end-1,:) = fliplr(IntD(2,:));    
%     
    % End conversion IntD to periodic
    
               
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
    E0 = y'*G'*G*y; 
    eRRK(1,1) = E0;
     
    iCount = 2; 
    t = dt; 

    while t <= tEnd+dt
        z1 = y; 
        [k1] = - D*IntD*z1; %        
        % Convert to periodic boundaries
        k1(1,1) = k1(2:3,1)'*IntD(3,4:5)' + ...
                k1(end-2:end-1)'*IntD(3,2:3)'; 
        k1(end,end) = k1(1,1);     
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = - D*IntD*z2;
        % Convert to periodic boundaries
        k2(1,1) = k2(2:3,1)'*IntD(3,4:5)' + ...
                k2(end-2:end-1)'*IntD(3,2:3)'; 
        k2(end,end) = k2(1,1);     


        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = - D*IntD*z3 ;
        % Convert to periodic boundaries
        k3(1,1) = k3(2:3,1)'*IntD(3,4:5)' + ...
                k3(end-2:end-1)'*IntD(3,2:3)'; 
        k3(end,end) = k3(1,1);     

        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = - D*IntD*z4 ;
        % Convert to periodic boundaries
        k4(1,1) = k4(2:3,1)'*IntD(3,4:5)' + ...
                k4(end-2:end-1)'*IntD(3,2:3)'; 
        k4(end,end) = k4(1,1);     

        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4; % + B(5,1)*k5;   
    
       
        switch RKFlag
            case 1
%                 gam = Nr/Dr; 
                U = y; 
                dU = BKsum; 
                
                AA = U'*G'*G*U ;
                BB = 2*(U'*G'*G*dU);
                CC = dU'*G'*G*dU;
                
                gam = (E0 - AA - BB)/(dt* CC); 
            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;       
        
        
        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        eRRK(iCount,1) = uNew'*G'*G*uNew; 
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end

function [uRRK6, tRRK6, gamRRK6, eRRK6] = RRK6(RKFlag,NElem,dh,dt,tEnd,Ordr,u0)
    
    % 6th order Runge Kutta Relaxation Verner
    D = div(Ordr,NElem,dh);
    G = grad(Ordr,NElem,dh);
    IntD = interpDMat(Ordr, NElem);     
 
     % Convert IntD to periodic
%     
    IntD(2:3,:) = zeros(2,NElem+2); 
    IntD(end-2:end-1,:) = zeros(2,NElem+2); 
    
    IntD(3,2:6) = IntD(4,3:7); IntD(3,end-1) = IntD(4,2); 
    IntD(end-2,:) = fliplr(IntD(3,:));    
    
    IntD(2,2:5) = IntD(4,4:7); IntD(2,end-2:end-1) = IntD(4,2:3); 
    IntD(end-1,:) = fliplr(IntD(2,:));
    
    % End conversion IntD to periodic

        
    C(1,1) = 0;
    C(2,1) = 1/6;
    C(3,1) = 4/15;
    C(4,1) = 2/3;
    C(5,1) = 5/6; 
    C(6,1) = 1;    
    C(7,1) = 1/15;
    C(8,1) = 1; 

    A(2,1) = 1/6;
    A(3,1) = 4/75;   A(3,2) = 16/75;
    A(4,1) = 5/6;   A(4,2) = -8/3;     A(4,3) = 5/2;    
    A(5,1) = -165/64; A(5,2) = 55/6; A(5,3) = -425/64;
        A(5,4) = 85/96;
    A(6,1) = 12/5; A(6,2) = -8; A(6,3) = 4015/612; 
        A(6,4) = -11/36; A(6,5) = 88/255;
    A(7,1) = -8263/15000; A(7,2) = 124/75; 
        A(7,3) = -643/680; A(7,4) = -81/250; 
        A(7,5) = 2484/10625; A(7,6) = 0;
    A(8,1) = 3501/1720; A(8,2) = -300/43; 
        A(8,3) = 297275/52632; A(8,4) = -319/2322; 
        A(8,5) = 24068/84065; A(8,6) = 0;
        A(8,7) = 3850/26703; 
                

    B(1,1) = 3/40; 
    B(2,1) = 0;
    B(3,1) = 875/2244;
    B(4,1) = 23/72;
    B(5,1) = 264/1955;
    B(6,1) = 0;
    B(7,1) = 125/11592; 
    B(8,1) = 43/616; 
    

     
    uRRK6(:,1) = u0; % initial value
    y = u0; 
    tRRK6(1,1) = 0; % initial time    
    gamRRK6(1,1) = 1;
    E0 = y'*G'*G*y; 
    eRRK6(1,1) = E0; 
    
    iCount = 2; 
    t = dt; 
    
    while t <= tEnd+dt
        z1 = y;        
        [k1] = - D*IntD*z1 ; %   
        % Convert to periodic boundaries
        k1(1,1) = k1(2:4,1)'*IntD(4,5:7)' + ...
                k1(end-3:end-1)'*IntD(4,2:4)'; 
        k1(end,end) = k1(1,1);     

        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;        
        [k2] = - D*IntD*z2 ; % 
        % Convert to periodic boundaries
        k2(1,1) = k2(2:4,1)'*IntD(4,5:7)' + ...
                k2(end-3:end-1)'*IntD(4,2:4)'; 
        k2(end,end) = k2(1,1);     

%
        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2);             
        [k3] = - D*IntD*z3; %
        % Convert to periodic boundaries
        k3(1,1) = k3(2:4,1)'*IntD(4,5:7)' + ...
                k3(end-3:end-1)'*IntD(4,2:4)'; 
        k3(end,end) = k3(1,1);     

% 
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = - D*IntD*z4 ; %
        % Convert to periodic boundaries
        k4(1,1) = k4(2:4,1)'*IntD(4,5:7)' + ...
                k4(end-3:end-1)'*IntD(4,2:4)'; 
        k4(end,end) = k4(1,1);     

        
        % Step 5    
        z5 = y + dt*(A(5,1)*k1 + A(5,2)*k2 + ...
            A(5,3)*k3 + A(5,4)*k4);   
        [k5] = - D*IntD*z5 ; %
        % Convert to periodic boundaries
        k5(1,1) = k5(2:4,1)'*IntD(4,5:7)' + ...
                k5(end-3:end-1)'*IntD(4,2:4)'; 
        k5(end,end) = k5(1,1);     

        
        % Step 6    
        z6 = y + dt*(A(6,1)*k1 + A(6,2)*k2 + ...
            A(6,3)*k3 + A(6,4)*k4 + A(6,5)*k5);   
        [k6] = - D*IntD*z6 ; %
        % Convert to periodic boundaries
        k6(1,1) = k6(2:4,1)'*IntD(4,5:7)' + ...
                k6(end-3:end-1)'*IntD(4,2:4)'; 
        k6(end,end) = k6(1,1);     
                
        % Step 7    
        z7 = y + dt*(A(7,1)*k1 + A(7,2)*k2 + ...
            A(7,3)*k3 + A(7,4)*k4 + A(7,5)*k5 + A(7,6)*k6);   
        [k7] = - D*IntD*z7 ; %
        % Convert to periodic boundaries
        k7(1,1) = k7(2:4,1)'*IntD(4,5:7)' + ...
                k7(end-3:end-1)'*IntD(4,2:4)'; 
        k7(end,end) = k7(1,1);     

        
        % Step 8    
        z8 = y + dt*(A(8,1)*k1 + A(8,2)*k2 + ...
            A(8,3)*k3 + A(8,4)*k4 + A(8,5)*k5 + A(8,6)*k6 ...
            + A(8,7)*k7);   
        [k8] = - D*IntD*z8 ; %
        % Convert to periodic boundaries
        k8(1,1) = k8(2:4,1)'*IntD(4,5:7)' + ...
                k8(end-3:end-1)'*IntD(4,2:4)'; 
        k8(end,end) = k8(1,1);     

              
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4 + B(5,1)*k5 + B(6,1)*k6 + ...
            B(7,1)*k7 + B(8,1)*k8;   
    
        switch RKFlag
            case 1                
                U = y; 
                dU = BKsum; 
                
                AA = U'*G'*G*U ;
                BB = 2*(U'*G'*G*dU);
                CC = dU'*G'*G*dU;
                
                gam = (E0 - AA - BB)/(dt* CC); 

            case 0
                gam = 1; 
        end
        
        uNew = y  + gam*dt*BKsum;             
                       

        uRRK6(:,iCount) = uNew;   
        y = uRRK6(:,iCount);
        tRRK6(iCount,1) = t;
        gamRRK6(iCount,1) = gam; 
        eRRK6(iCount,1) = uNew'*G'*G*uNew; 
        
        iCount = iCount+1; 
        t = t + gam*dt;      
        
        
    end
    

end

function energyOut = EnergyCalc(uRRK,dh)

    
    for i = 1:size(uRRK,1)
       energyOut(i,1) = uRRK(i,1)-uRRK(1,1); 
    end
    
end



