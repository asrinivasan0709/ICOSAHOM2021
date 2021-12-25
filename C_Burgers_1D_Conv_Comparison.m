clc; clear; close all;
% Convergence comparison for 1D Burgers eqn.
% ASrinivasan, 19Dec21, Use for MOLE library
% Compares the spatial Mimetic 4th with Wicker's scheme
% Time integration is performed using RK4 & RRK. 
% Periodic boundaries are assumed
% Exact soln -- https://arxiv.org/pdf/1503.09079.pdf

addpath('../mole_MATLAB')

Ne = [50,100,200]'; %,400,800,1600]'; 

for i = 1:size(Ne,1)
   NElem = Ne(i,1);
   [nOutWIC, nOutMIM, nOutWIC6, nOutMIM6] = NormCalcs(NElem); 
   normWIC(i,:) = nOutWIC;
   normMIM(i,:) = nOutMIM;
   normWIC6(i,:) = nOutWIC6;
   normMIM6(i,:) = nOutMIM6;
    
end

normMat = [normWIC; normMIM; normWIC6; normMIM6]; 


%%% Plots
f1 = figure; %RRK, inf norm 
loglog(normWIC(:,2), normWIC(:,5), '--sb', 'LineWidth', 3)
hold on
loglog(normMIM(:,2), normMIM(:,5), '-+k', 'LineWidth', 3)
hold on   %line for p=4
loglog([1E-3, 2E-3], [4e-7, 64E-7], '-*k', 'linewidth', 1)
hold on
loglog(normWIC6(:,2), normWIC6(:,5), '--sr', 'LineWidth', 3)
hold on
loglog(normMIM6(:,2), normMIM6(:,5), '-+r', 'LineWidth', 3)
hold on   %line for p=6
loglog([1E-3, 2E-3], [4e-10, 256E-10], '-ok', 'linewidth', 1)

set(gca,'FontSize',10)
legend('WIC4-RK4', 'MIM4-RK4',  ...
        'p=4', ...
        'WIC6-RK6', 'MIM6-RK6', 'p=6', ...
        'location','best');
title({'1D Burgers Eqn u_t + \nabla\cdot(u^2/2) = 0', ...
    'Mimetic vs Finite Difference RK Methods'}); 
set(gca, 'YScale', 'log');
xlabel('\Deltax')
ylabel('Error Convergence, \infty-norm')
movegui(f1,'southwest');
grid on



function [nOutWIC, nOutMIM, nOutWIC6, nOutMIM6] = NormCalcs(NElem)
 
xmin = 0; xmax = 1;

dh = (xmax-xmin)/NElem; 
tEnd = 0.1;  %10
CFL = 0.1; %2.02
dt =  CFL*dh; 

% init val
xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
xNod = [xmin:dh:xmax]';
[u0Grid, u0Nod, uExGrid, uExNod] = initval(xGrid, xNod, tEnd); 


% Wicker's spatial calcs
%%% Wicker 4th order scheme
    Ordr = 4;   
    [uWICRK4, tWICRK4,gamWICRK4, eWICRK4] = WIC(0,NElem,dh,dt,tEnd,u0Nod);         
    energyWICRK4 = EnergyCalc(eWICRK4, dh);
    uWICRK4end = InterpEnd(uWICRK4, tEnd, tWICRK4); 
    normWICRK4inf = norm(uWICRK4end-uExNod, 'inf');
    normWICRK4L2 = norm(uWICRK4end-uExNod, 2);
    nOutWIC = [NElem, dh, dt, normWICRK4inf, ...
                normWICRK4L2]; 


%%% Mimetic 4th order scheme
    Ordr = 4;     
    [uRK4, tRK4,gamRK4, eRK4] = RRK(0,NElem,dh,dt,tEnd,Ordr,u0Grid);    
    energyRK4 = EnergyCalc(eRK4,dh);
    uRK4end = InterpEnd(uRK4, tEnd, tRK4);
    normRK4inf = norm(uRK4end-uExGrid, 'inf');
    normRK4L2 = norm(uRK4end-uExGrid, 2);

    nOutMIM = [NElem, dh, dt, normRK4inf, ...
                normRK4L2]; 

                        
%%% Wicker 6th order scheme
    Ordr = 6;   
    [uWICRK6, tWICRK6,gamWICRK6, eWICRK6] = WIC6(0,NElem,dh,dt,tEnd,u0Nod);    
    energyWICRK6 = EnergyCalc(eWICRK6, dh);
    uWICRK6end = InterpEnd(uWICRK6, tEnd, tWICRK6); 
    normWICRK6inf = norm(uWICRK6end-uExNod, 'inf');
    normWICRK6L2 = norm(uWICRK6end-uExNod, 2);
    
    nOutWIC6 = [NElem, dh, dt, normWICRK6inf, ...
                normWICRK6L2]; 

             
%%% Mimetic 6th order scheme
%%% RRK - 6th order
    Ordr = 6; 
    
    [uRK6, tRK6,gamRK6, eRK6] = RRK6(0,NElem,dh,dt,tEnd,Ordr,u0Grid);    
    energyRK6 = EnergyCalc(eRK6,dh);
    uRK6end = InterpEnd(uRK6, tEnd, tRK4);
    normRK6inf = norm(uRK6end-uExGrid, 'inf');
    normRK6L2 = norm(uRK6end-uExGrid, 2);

    nOutMIM6 = [NElem, dh, dt, normRK6inf, ...
                normRK6L2]; 
                   
            
end         
    

function [u0Grid, u0Nod, uExGrid, uExNod] = initval(xGrid, xNod, tEnd)

            
    u0Grid = zeros(size(xGrid,1),1);        
    u0Grid = sin(2*pi*xGrid); 

    u0Nod = zeros(size(xNod,1),1);        
    u0Nod = sin(2*pi*xNod); 
    % Exact solution at end, Grid
    
    beta = 1; 

    for j = 1:size(xGrid,1)
        xGridIn  = xGrid(j,1);
        myfun = @(x,xGridIn,beta,tEnd) ...
            xGridIn - (x - log(1 - beta*tEnd*sin(2*pi*x)))/beta;   % parameterized function
        fun = @(x) myfun(x,xGridIn,beta,tEnd);    % function of x alone
        x = fzero(fun,0.1);
        xOut(j,1) = x; 
    end

    uExGrid = sin(2*pi*xOut)./(1 - beta*tEnd*sin(2*pi*xOut)); 

    % Exact at end, Nodal Grid       

    for j = 1:size(xNod,1)                
        xNodIn  = xNod(j,1);
        myfunNod = @(x,xNodIn,beta,tEnd) ...
            xNodIn - (x - log(1 - beta*tEnd*sin(2*pi*x)))/beta;   % parameterized function
        funNod = @(x) myfunNod(x,xNodIn,beta,tEnd);    % function of x alone
        xNodVal = fzero(funNod,0.1);
        xOutNod(j,1) = xNodVal; 
    end

    uExNod = sin(2*pi*xOutNod)./(1 - beta*tEnd*sin(2*pi*xOutNod)); 
  
        
end

function [uWIC, tWIC, gamWIC, eWIC] = WIC(RKFlag,NElem,dh,dt,tEnd,u0)
    
% Wicker & Skamarock spatial with RK4
          
    [FP, FM] = FluxWIC(NElem); 
    FP = 1/dh*FP;
    FM = 1/dh*FM; 
    
               
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

     
    uWIC(:,1) = u0; % initial value
    y = u0; 
    tWIC(1,1) = 0; % initial time    
    gamWIC(1,1) = 1;
    E0 = norm(y,2)^2; 
    eWIC(1,1) = E0; 
    
    iCount = 2; 
    t = dt; 

    while t <= tEnd+dt
        z1 = y; 
        [k1] = z1.^2 - z1.*(FP-FM)*z1; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = z2.^2 - z2.*(FP-FM)*z2;

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = z3.^2 - z3.*(FP-FM)*z3 ;
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = z4.^2 - z4.*(FP-FM)*z4 ;     
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4; %   
    
       
        switch RKFlag
            case 1
                U = y; 
                dU = BKsum; 
                
                AA = U'*U ;
                BB = 2*(U'*dU);
                CC = dU'*dU;               
                
                gam = (E0 - AA - BB)/(dt* CC); 

            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;
                
        uWIC(:,iCount) = uNew;   
        y = uWIC(:,iCount);
        tWIC(iCount,1) = t;
        gamWIC(iCount,1) = gam; 
        eWIC(iCount,1) = norm(uNew,2)^2; 
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end

function [uWIC6, tWIC6, gamWIC6, eWIC6] = WIC6(RKFlag,NElem,dh,dt,tEnd,u0)
    
% Wicker & Skamarock 6th order spatial with RK4
          
    [FP, FM] = FluxWIC6(NElem); 
    FP = 1/dh*FP;
    FM = 1/dh*FM; 
    
               
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
    
     
    uWIC6(:,1) = u0; % initial value
    y = u0; 
    tWIC6(1,1) = 0; % initial time    
    gamWIC6(1,1) = 1;
    E0 = norm(y,2)^2; 
    eWIC6(1,1) = E0; 
    
    iCount = 2; 
    t = dt; 
    
    
    while t <= tEnd+dt
        z1 = y; 
        [k1] = z1.^2 - z1.*(FP-FM)*z1; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = z2.^2 - z2.*(FP-FM)*z2;

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = z3.^2 - z3.*(FP-FM)*z3 ;
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = z4.^2 - z4.*(FP-FM)*z4 ;
        

        % Step 5    
        z5 = y + dt*(A(5,1)*k1 + A(5,2)*k2 + A(5,3)*k3 + ...
            A(5,4)*k4);         
        [k5] = z5.^2 - z5.*(FP-FM)*z5 ;


        % Step 6    
        z6 = y + dt*(A(6,1)*k1 + A(6,2)*k2 + A(6,3)*k3 + ...
            A(6,4)*k4 + A(6,5)*k5);         
        [k6] = z6.^2 - z6.*(FP-FM)*z6 ;

        % Step 7    
        z7 = y + dt*(A(7,1)*k1 + A(7,2)*k2 + A(7,3)*k3 + ...
            A(7,4)*k4 + A(7,5)*k5 + A(7,6)*k6);         
        [k7] = z7.^2 - z7.*(FP-FM)*z7 ;
 
        
        % Step 8    
        z8 = y + dt*(A(8,1)*k1 + A(8,2)*k2 + A(8,3)*k3 + ...
            A(8,4)*k4 + A(8,5)*k5 + A(8,6)*k6 + A(8,7)*k7) ;         
        [k8] = z8.^2 - z8.*(FP-FM)*z8 ;
        
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4 + B(5,1)*k5 + B(6,1)*k6 + B(7,1)*k7 + ...
            B(8,1)*k8;   
    
       
        switch RKFlag
            case 1
                U = y; 
                dU = BKsum; 
                
                AA = U'*U ;
                BB = 2*(U'*dU);
                CC = dU'*dU;               
                
                gam = (E0 - AA - BB)/(dt* CC);                 

            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;
                
        uWIC6(:,iCount) = uNew;   
        y = uWIC6(:,iCount);
        tWIC6(iCount,1) = t;
        gamWIC6(iCount,1) = gam; 
        eWIC6(iCount,1) = norm(uNew,2)^2; 
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end


function [FP, FM] = FluxWIC(NElem)

% Flux matrices for Wicker & Skamarock, 4th order spatial
% Uses periodic boundary

    FM = zeros(NElem+1, NElem+1); 
    for i = 3:NElem
       FM(i,i-2:i+1) = [-1,7,7,-1]; 
    end

    FM(1,1:2) = [7,-1]; FM(1,end-2:end-1) = [-1,7]; 
    FM(2,1:3) = [7,7,-1]; FM(2,end-1) = -1;
    FM(end,:) = fliplr(FM(2,:)); 
    FM = 1/12*FM; 

    FP = zeros(NElem+1, NElem+1); 
    for i = 2:NElem-1
       FP(i,i-1:i+2) = [-1,7,7,-1]; 
    end
    FP(1,1:3) = [7,7,-1]; FP(1,end-1) = -1;
    FP(end-1,:) = fliplr(FP(1,:));
    FP(end,2:3) = [7,-1]; FP(end,end-1:end) = [-1,7];
    FP = 1/12*FP; 
end

function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tEnd,Ordr,u0)
    
% Relaxation RK4
   
    D = div(Ordr,NElem,dh);
    IntD = interpDMat(Ordr, NElem);  
    
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
    E0 = y'*y;  
    eRRK(1,1) = E0;
     
    iCount = 2; 
    t = dt; 
    DIntD = D*IntD; 
    
    while t <= tEnd+dt
        z1 = y; 
        [k1] = z1.^2 - z1.*DIntD*z1; %        
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = z2.^2 - z2.*DIntD*z2;

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = z3.^2 - z3.*DIntD*z3 ;
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = z4.^2 - z4.*DIntD*z4 ;
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4;    
    
      
        switch RKFlag
            case 1
                U = y; 
                dU = BKsum; 
                
                AA = U'*U; 
                BB = 2*(U'*dU); 
                CC = dU'*dU; 
                                
                gam = (E0 - AA - BB)/(dt*CC); % 
                
            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;              
               
        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        eRRK(iCount,1) = (uNew'*uNew); 
        
        iCount = iCount+1; 
        t = t + gam*dt;        %   
        
    end   

end

function [uRRK6, tRRK6, gamRRK6, eRRK6] = RRK6(RKFlag,NElem,dh,dt,tEnd,Ordr,u0)
    
% Relaxation RK6
       
    D = div(Ordr,NElem,dh);
    IntD = interpDMat(Ordr, NElem);  
    
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
    E0 = y'*y;  
    eRRK6(1,1) = E0;
     
    iCount = 2; 
    t = dt; 
    DIntD = D*IntD; 
    
    while t <= tEnd+dt
        z1 = y; 
        [k1] = z1.^2 - z1.*DIntD*z1; %        
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = z2.^2 - z2.*DIntD*z2;


        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = z3.^2 - z3.*DIntD*z3 ;

        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = z4.^2 - z4.*DIntD*z4 ;

        % Step 5    
        z5 = y + dt*(A(5,1)*k1 + A(5,2)*k2 + A(5,3)*k3 + A(5,4)*k4);         
        [k5] = z5.^2 - z5.*DIntD*z5 ;

        % Step 6    
        z6 = y + dt*(A(6,1)*k1 + A(6,2)*k2 + A(6,3)*k3 + A(6,4)*k4 + ...
            A(6,5)*k5);         
        [k6] = z6.^2 - z6.*DIntD*z6 ;

        % Step 7    
        z7 = y + dt*(A(7,1)*k1 + A(7,2)*k2 + A(7,3)*k3 + A(7,4)*k4 + ...
            A(7,5)*k5 + A(7,6)*k6);         
        [k7] = z7.^2 - z7.*DIntD*z7 ;

        % Step 8    
        z8 = y + dt*(A(8,1)*k1 + A(8,2)*k2 + A(8,3)*k3 + A(8,4)*k4 + ...
            A(8,5)*k5 + A(8,6)*k6 + A(8,7)*k7);         
        [k8] = z8.^2 - z8.*DIntD*z8 ;
       
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4 + B(5,1)*k5 + B(6,1)*k6 + B(7,1)*k7 + ...
            B(8,1)*k8 ;    
    
      
        switch RKFlag
            case 1
                U = y; 
                dU = BKsum; 
                
                AA = U'*U; 
                BB = 2*(U'*dU); 
                CC = dU'*dU; 
                                
                gam = (E0 - AA - BB)/(dt*CC); % 
                
            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;              
               
        uRRK6(:,iCount) = uNew;   
        y = uRRK6(:,iCount);
        tRRK6(iCount,1) = t;
        gamRRK6(iCount,1) = gam; 
        eRRK6(iCount,1) = (uNew'*uNew); 
        
        iCount = iCount+1; 
        t = t + gam*dt;        %   
        
    end   

end


function [FP, FM] = FluxWIC6(NElem)

% Flux matrices for Wicker & Skamarock, 6th order spatial
% Uses periodic boundary

    FM = zeros(NElem+1, NElem+1); 
    for i = 4:NElem-1
       FM(i,i-3:i+2) = [1,-8,37,37,-8,1]; 
    end

    FM(1,1:3) = [37,-8,1]; FM(1,end-3:end-1) = [1,-8,37]; 
    FM(2,1:4) = [37,37,-8,1]; FM(2,end-2:end-1) = [1,-8];
    FM(3,1:5) = [-8,37,37,-8,1]; FM(3,end-1) = 1; 
     
    FM(end-1,end-4:end) = [1,-8,37,37,-8]; 
    FM(end-1,2) = 1; 
    
    FM(end, end-3:end) = [1,-8,37,37];
    FM(end, 2:3) = [-8,1]; 
%     FM(end-2,:) = fliplr(FM(3,:)); 
    
    FM = 1/60*FM; 

        
    FP = zeros(NElem+1, NElem+1); 
    for i = 3:NElem-2
       FP(i,i-2:i+3) = [1,-8,37,37,-8,1]; 
    end
    FP(1,1:4) = [37,37,-8,1]; FP(1,end-2:end-1) = [1,-8];
    FP(2,1:5) = [-8,37,37,-8,1]; FP(2,end-1) = 1; 
    
    FP(end-2,:) = fliplr(FP(2,:));
    FP(end-1,:) = fliplr(FP(1,:)); 
    FP(end,2:4) = [37,-8,1]; FP(end,end-2:end) = [1,-8,37]; 
    FP = 1/60*FP; 
    
    
    
end


function energyOut = EnergyCalc(uRRK,dh)
    
    for i = 1:size(uRRK,1)
       energyOut(i,1) = uRRK(i,1)-uRRK(1,1); 
    end
    
end

function uWICRRKend = InterpEnd(uWIC, tEnd, tWIC)
    % Ending value for uRRK, if tEnd > tWICRRK(end)
    if tWIC(end,1) ~= tEnd
       u2 = uWIC(:,end);  
       u1 = uWIC(:,end-1);
       t2 = tWIC(end,1);
       t1 = tWIC(end-1,1);
       uWICRRKend = (u2-u1)*(tEnd-t1)/(t2-t1) + u1;    
    else 
       uWICRRKend = uWIC(:,end);
    end

end



