%% ELEC 4700- Assignment 3
% Saifuddin Mohammed, #101092039

%% P1- Modification of the Monte Carlo Model from Assignment 1


set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth',2); 


clear all;
close all;

%Given Paramters
T = 300;  %Semicondctor Temperature 
C_m0 = 9.10938356e-31; %Rest mass of the Electron
C_m = 0.26*C_m0; %Given Effective Mass of Electrons
C_q = -1.60217662e-19; % Charge on electron


%Given Nominal Dimensions of Semiconductor 200 nm x 100 nm
X_R = 200e-9;
Y_R = 100e-9;

% Updated Conditions include the voltage being applied in x direction to be
% 0.1v
V_X = 0.1; 

% No changes to the voltage being applied in Y direction
V_Y = 0;


% The electron concentration is given as
elec_conc = 1e15*100^2;

% Calculation of Thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt
V_T = sqrt(2*K*T/C_m);


electron_population = 30000;
electron_num = 50;

% Setting the Step Size
TStep = Y_R/V_T/100;
iter = 200;

% The scattering probabily is given by
P_scat = 1 - exp(-TStep/0.2e-12);

%Using Maxewell-Boltzmann Distribution to Generate random velocities
Distr_MB = makedist('Normal', 0, sqrt(K*T/C_m));



animation_plot = 0;


% The scattering probabily is given by
P_scat = 1 - exp(-TStep/0.2e-12);


% Calculating the Electric Field 
Electric_Field = V_X/X_R; 
fprintf(' The Electric field experienced by the electrons is %i\n',Electric_Field);

% Calculating the Force
Force = Electric_Field*C_q;
fprintf(' The Force experienced by the electrons is %i\n',Force);

% Calculating the Accelaration 
Accelaration = Force/C_m; 
fprintf(' The accelaration of the electrons is %i\n',Accelaration);

%The Electric Field and the Force in the Y direction calculated as 
Electric_FieldY = V_Y/Y_R
ForceY = C_q*Electric_FieldY


% Calculating velocity at each time step of the electrons
Del_X = Force*TStep/C_m;
Del_Y = ForceY*TStep/C_m;
Del_X = Del_X.*ones(electron_population,1);
Del_Y = Del_Y.*ones(electron_population,1);


% Defining the Specular conditions at the Top and the Bottom boundary of
% the region
Spec_T = 0;
Spec_B = 0;



%Initialization of the Position of Random Particles
%Setting up x and positions
post = zeros(electron_population, 4);



%Keeping track of trajectories
Trj = zeros(iter, electron_num*2);


%Recording the temeperatures
Temp = zeros(iter,1);


% The current density can be stored as
Curr_Dens = zeros(iter,2);



%Initialization of the Position of Random Particles
%Setting up x and positions
for ni = 1:electron_population
    
    theta = rand*2*pi;
    
    post(ni,:) = [X_R*rand Y_R*rand random(Distr_MB) random(Distr_MB)];
    
end


% Plotting the trajectories of the electrons
figure(1);
plot([],[]);
title('Electron Trajectories, SM101092039'); 
xlabel('X position');
ylabel('Y position');



% Generating a plot of current v time
figure(2);
current_plot =  animatedline('Color','r','LineWidth',3);
set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5; 
a_x.GridColor = [1, 1, 1]; 
title('Plot of Current v Time, SM101092039');
xlabel('Time (s)');
ylabel('Current (A/m)');
grid on;


% Simulating the Monte Carlo Model

for ni = 1:iter                        % Performs the simulation with Random Velocities with updating positions. 
    
    
    %Utilzing the velocities of the electrons to determine positions
    post(:,3) = post(:,3) + Del_X;
    post(:,4) = post(:,4) + Del_Y;
    
    
     %Boundary Conditions 
    post(:,1:2) = post(:,1:2) + TStep.*post(:,3:4);

    nk = post(:,1) > X_R;
    post(nk,1) = post(nk,1) - X_R;
    
    nk = post(:,1) < 0;
    post(nk,1) = post(nk,1) + X_R;
    
    nk = post(:,2) > Y_R;

    
    % At the Top of Region Boundary 
    if(Spec_T)                                  % for specular condition 
        
        post(nk,2) = 2*Y_R - post(nk,2);
        post(nk,4) = -post(nk,4);
        
        
    else 
        
        post(nk,2) = 100e-9;                                        %Diffusive condition has been met
        Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        post(nk,3) = Z.*cos(theta);
        post(nk,4) = -abs(Z.*sin(theta));
        
        
    end
    
    
    nk = post(:,2) < 0;
    
    
    % At the bottom of the region boundary 
      if(Spec_B)
                                                          % for specular condition 
        post(nk,2) = -post(nk,2);
        post(nk,4) = -post(nk,4);
        
        
      else  
                                                           
        post(nk,2) = 0;                                       %Diffusive condition has been met
        Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
        
        
        theta = rand([sum(nk),1])*2*pi;
        post(nk,3) = Z.*cos(theta);
        post(nk,4) = abs(Z.*sin(theta));
        
        
    end
    
    
    
    
    nk = rand(electron_population, 1) < P_scat;           % Random dsitribution of particles within the scattering probability limit
    post(nk,3:4) = random(Distr_MB, [sum(nk),2]);          %Scattering particles using the exponential scattering probability 
    
    
    
    %Calculating and Recording the Temperature of the electrons
    Temp(ni) = (sum(post(:,3).^2) + sum(post(:,4).^2))*C_m/K/2/electron_population;           
    
   
    
    for nk=1:electron_num
        
        %Storing the positions of the electrons on the trajector 
        Trj(ni, (2*nk):(2*nk+1)) = post(nk, 1:2);
        
    end
    
    
    
    
    Curr_Dens(ni, 1) = C_q.*elec_conc.*mean(post(:,3));             % Calculate the current density 
    
    Curr_Dens(ni, 2) = C_q.*elec_conc.*mean(post(:,4));
    
    

  
    addpoints(current_plot, TStep.*ni, Curr_Dens(ni,1));          % Combining all the locations at points for the total current 

    
    %Plotting the electron trajectories and temeperature variance
    % The animation updates for every 10 iterations
    if(animation_plot && mod(ni,10) == 0)
        figure(1);
        hold off;
        plot(post(1:electron_num,1), post(1:electron_num,2), 'o');       
        hold on;
        title(' Electron Trajectories, SM101092039'); 
        xlabel('X position)');
        ylabel('Y position');
        pause(0.05);
        
    end
    
end



%Final plotting of the trajectories after all the iterations have been
%completed

figure(1);
title('Electron Trajectories, SM101092039');
xlabel('X position');
ylabel('Y position');

grid on;
hold on;

%Storing the trajectory after completion 
for ni=1:electron_num
    
    plot(Trj(:,ni*2), Trj(:,ni*2+1), '.');
    
end



% Plot the Density Map 

%Utilizing the Hist3 command to generate the density map 
elec_conc = hist3(post(:,1:2),[200 100])';

 %Utilzing the bins from the historgram to plot the the Temperature and Density map 
hist_bins = 10;    

[x_elec y_elec] = meshgrid(round(-hist_bins/2):round(hist_bins/2), round(-hist_bins/2):round(hist_bins/2));
M=exp(-x_elec.^2/(2*1^2)-y_elec.^2/(2*1^2));
M=M./sum(M(:));

figure(3);
elec_conc = conv2(elec_conc,M,'same');       %2D convolution for generating a 2D map;
elec_conc = elec_conc/(Y_R./size(elec_conc,1)*X_R./size(elec_conc,2));
surf(conv2(elec_conc,M,'same'));
title('Electron Density, SM101092039');
xlabel('X position');
ylabel('Y position');


%% Temperature Map

X_Total_Temp = zeros(ceil(X_R/1e-9),ceil(Y_R/1e-9));
Y_Total_Temp = zeros(ceil(X_R/1e-9),ceil(Y_R/1e-9));
Total_Temp = zeros(ceil(X_R/1e-9),ceil(Y_R/1e-9));




for ni=1:electron_population
    
   % Determination of the bin 
    x_elec = floor(post(ni,1));
    y_elec = floor(post(ni,2));
    
    
    if(x_elec==0)
        
        x_elec = 1;
        
    end
    
    if(y_elec==0)
        
        y_elec= 1;
        
    end
    
   
    Y_Total_Temp(x_elec,y_elec) = Y_Total_Temp(x_elec,y_elec) + post(ni,3)^2;
    
    X_Total_Temp(x_elec,y_elec) = X_Total_Temp(x_elec,y_elec) + post(ni,4)^2;
    
     %Summing all the components of the velocity 
    Total_Temp(x_elec,y_elec) = Total_Temp(x_elec,y_elec) + 1;
    
end


%The temperature can be calculated as:

temp = (X_Total_Temp + Y_Total_Temp).*C_m./K./2./Total_Temp;
temp= temp;  




%Generate The temperature plot
hist_bins = 10;
[x_elec , y_elec] = meshgrid(round(-hist_bins/2):round(hist_bins/2), round(-hist_bins/2):round(hist_bins/2));
M=exp(-x_elec.^2/(2*1^2)-y_elec.^2/(2*1^2));
M=M./sum(M(:));
figure(4);
imagesc(conv2(temp,M,'same'));             % Utilizng surf to plot the 2D convolution map
view(0,90)
title('Temperature Map, SM101092039');
xlabel('X position');
ylabel('Y position');







%% P2 - Utilizing the Finite Difference Method to setup Electric Field for the Monte Carlo Model

clear

% Setting up Length and Width as per the ratio L/W
L = 1;
W = 1;

V_0 = 1;    %Assigning V0 a value of 1;


nx = 100*L;  
ny = 100*W;

W_BN = 0.4;     %Width of the rectangle box
L_BN = 0.4;     %Length of the rectangle box


sig_1 = 1;        %Value of sigma outside the box
sig_2 = 1e-2;     %Value of sigma inside the box





C = zeros(ny,nx);   % Map of Conductivity 




for h = 1:ny
    
    for g = 1:nx
        
        %Assigning the values of sigma as per the location of the element 
        
        
        if(g >= nx*W_BN && g <= nx-nx*W_BN && (h >= ny-ny*L_BN || h <= ny*L_BN)) 
            
            
            % If inside the box, then sigma = 1e^-2
            
            C(h,g) = sig_2;
            
       
        else                                         %This is for outside the box. sigma = 1
            
            C(h,g) = sig_1;
            
        end
        
        
    end
    
    
end



%  Creating the G matrix and B vector for the GV = F solution

 G = sparse(nx*ny);
 F = zeros(nx*ny,1);



for g = 1:nx
    
    for h = 1:ny
        
       %Mapping of the nodes equation  
        n = h + (g - 1)*ny;
        
        
        
       %Local Mapping of the nodes around g and h
        nxm = h + (g - 2)*ny;
        nxp = h + g*ny;
        nym = (h - 1) + (g - 1)*ny;
        nyp = (h + 1) + (g - 1)*ny;
        
       
        
        if(g == 1 || g == nx) 
            
            
            G(n,n) = 1;       %Left Side Set Voltage 
            
       
            
       
        elseif (h == 1)    %Evalutation at the bottom region 
            
            U_Y = (C(h,g)+C(h+1,g))/2;
            U_X = (C(h,g)+C(h,g+1))/2;
            UX_D = (C(h,g)+C(h,g-1))/2;
            
            G(n,n) = -(U_Y + U_X + UX_D);
            G(n,nyp) = U_Y;
            G(n,nxp) = U_X;
            G(n,nxm) = UX_D;
            
            
       % Evaluation at the top region   
       
       elseif (h == ny)
           
              
              YDY_E = (C(h,g)+C(h-1,g))/2;
              UDX_E = (C(h,g)+C(h,g+1))/2;
              DDX_E = (C(h,g)+C(h,g-1))/2;
              
              G(n,n) = -(YDY_E + UDX_E + DDX_E);
              G(n,nym) = YDY_E;
              G(n,nxp) = UDX_E;
              G(n,nxm) = DDX_E;
              
            
            
        else 
            
            
           %The finite difference method is being applied to evaluate the
           %regions potetntial outside the range as well. 
           
            U_Y = (C(h,g)+C(h+1,g))/2;
            D_Y = (C(h,g)+C(h-1,g))/2;
            U_X = (C(h,g)+C(h,g+1))/2;
            UX_D = (C(h,g)+C(h,g-1))/2;
            
            
            G(n,n) = -(U_Y + D_Y + U_X + UX_D);
            G(n,nyp) = U_Y;
            G(n,nym) = D_Y;
            G(n,nxp) = U_X;
            G(n,nxm) = UX_D;
            
            
        end
        
        
    end
    
    
end





for g = 1:nx
    
    for h = 1:ny
     
        
        
        %Node Mapping Equation 
        n = h + (g - 1)*ny;
       
        
        if (g == 1) %Indicating a shift towards left so the value must be set equal to V_0
            
          
            F(n) = V_0;
            
            
        end
        
        
        
    end
    
end

% Utilizing  GV = F to solve the equation 

V = G\F;


for g = 1:nx
    
    for h = 1:ny
        
        % Node mapping to put entries into the correct place
        n = h + (g - 1)*ny;
        
        Vmap(h,g) = V(n);
        
        
    end
end



% Plotting the Voltage V across the region
figure(5)
surf(Vmap)

xlabel('L (um)')
ylabel('W (um)')
title({'Voltage V(x,y) Plot, SM101092039'})



[E_x,E_y] = gradient(-Vmap);   %Electric Field of the regions can be determined from using the gradient




% Plotting the Electric Field over the region 


figure(6)
quiver(E_x,E_y)
xlabel('L (um)')
ylabel('W (um)')
title({'Electric Field in the Region, SM101092039'})



%Given Paramters
T = 300;  %Semicondctor Temperature 
C_m0 = 9.10938356e-31; %Rest mass of the Electron
C_m = 0.26*C_m0; %Given Effective Mass of Electrons
C_q = -1.60217662e-19; % Charge on electron


%Given Nominal Dimensions of Semiconductor 200 nm x 100 nm
X_R = 200e-9;
Y_R = 100e-9;

% Updated Conditions include the voltage being applied in x direction to be
% 0.1v
V_X = 0.1; 

% No changes to the voltage being applied in Y direction
V_Y = 0;


% The electron concentration is given as
elec_conc = 1e15*100^2;

% Calculation of Thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt
V_T = sqrt(2*K*T/C_m);


electron_population = 10000;
electron_num = 100;

% Setting the Step Size
TStep = Y_R/V_T/100;
iter = 200;

% The scattering probabily is given by
P_scat = 1 - exp(-TStep/0.2e-12);

%Using Maxewell-Boltzmann Distribution to Generate random velocities
Distr_MB = makedist('Normal', 0, sqrt(K*T/C_m));



animation_plot = 0;


% The scattering probabily is given by
P_scat = 1 - exp(-TStep/0.2e-12);


% Calculating the Electric Field 
Electric_Field = V_X/X_R; 
fprintf(' The Electric field experienced by the electrons is %i\n',Electric_Field);

% Calculating the Force
Force = Electric_Field*C_q;
fprintf(' The Force experienced by the electrons is %i\n',Force);

% Calculating the Accelaration 
Accelaration = Force/C_m; 
fprintf(' The accelaration of the electrons is %i\n',Accelaration);

%The Electric Field and the Force in the Y direction calculated as 
Electric_FieldY = V_Y/Y_R
ForceY = C_q*Electric_FieldY


% Calculating velocity at each time step of the electrons
Del_X = Force*TStep/C_m;
Del_Y = ForceY*TStep/C_m;
Del_X = Del_X.*ones(electron_population,1);
Del_Y = Del_Y.*ones(electron_population,1);


% Defining the Specular conditions at the Top and the Bottom boundary of
% the region
Spec_T = 0;
Spec_B = 0;



%Initialization of the Position of Random Particles
%Setting up x and positions
post = zeros(electron_population, 4);



%Keeping track of trajectories
Trj = zeros(iter, electron_num*2);


%Recording the temeperatures
Temp = zeros(iter,1);


st_size = 1e-9;
boxes = st_size.*[80 120 0 40; 80 120 60 100];
spec_boxes = [0 1];

%Initialization of the Position of Random Particles
%Setting up x and positions
for ni = 1:electron_population
    
    theta = rand*2*pi;
    
    post(ni,:) = [X_R*rand Y_R*rand random(Distr_MB) random(Distr_MB)];
    
end



% Simulating the Monte Carlo Model

for ni = 1:iter                        % Performs the simulation with Random Velocities with updating positions. 
    
    
    %Utilzing the velocities of the electrons to determine positions
    post(:,3) = post(:,3) + Del_X;
    post(:,4) = post(:,4) + Del_Y;
    
    
     %Boundary Conditions 
    post(:,1:2) = post(:,1:2) + TStep.*post(:,3:4);

    nk = post(:,1) > X_R;
    post(nk,1) = post(nk,1) - X_R;
    
    nk = post(:,1) < 0;
    post(nk,1) = post(nk,1) + X_R;
    
    nk = post(:,2) > Y_R;

    
    % At the Top of Region Boundary 
    if(Spec_T)                                  % for specular condition 
        
        post(nk,2) = 2*Y_R - post(nk,2);
        post(nk,4) = -post(nk,4);
        
        
    else 
        
        post(nk,2) = 100e-9;                                        %Diffusive condition has been met
        Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        post(nk,3) = Z.*cos(theta);
        post(nk,4) = -abs(Z.*sin(theta));
        
        
    end
    
    
    nk = post(:,2) < 0;
    
    
    % At the bottom of the region boundary 
      if(Spec_B)
                                                          % for specular condition 
        post(nk,2) = -post(nk,2);
        post(nk,4) = -post(nk,4);
        
        
      else  
                                                           
        post(nk,2) = 0;                                       %Diffusive condition has been met
        Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
        
        
        theta = rand([sum(nk),1])*2*pi;
        post(nk,3) = Z.*cos(theta);
        post(nk,4) = abs(Z.*sin(theta));
        
        
    end
    
    
    %Figuring out if the particles have moved into to the box. Updates the
  %position and restores the location of the particles
    
    for nk= 1:electron_num
        bottle_neck = box_no(post(nk,1:2), boxes);
        
        
        % Checking for the collision with a box and determining the
        % location of the box of collision. 
        
        while(bottle_neck ~= 0)
            
            dist_X = 0;                  %Finding and updating the X position
            
            X_updated = 0;
            
            
            if(post(nk,3) > 0)
                
                dist_X = post(nk,1) - boxes(bottle_neck,1);
                X_updated = boxes(bottle_neck,1);
                
            else
                
                dist_X = boxes(bottle_neck,2) - post(nk,1);
                X_updated = boxes(bottle_neck,2);
                
                
            end

            dist_Y = 0;                  %Finding and updating the Y position
            Y_updated = 0;
            
            if(post(nk,4) > 0)
                
                dist_Y = post(nk,2) - boxes(bottle_neck, 3);
                Y_updated = boxes(bottle_neck, 3);
                
            else
                
                dist_Y = boxes(bottle_neck, 4) - post(nk,2);
                Y_updated = boxes(bottle_neck, 4);
                
            end

            if(dist_X < dist_Y)
                
                post(nk,1) = X_updated;
                
                if(~spec_boxes(bottle_neck))
                    
                    sgn = -sign(post(nk,3));
                    Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
                    
                    theta = rand()*2*pi;
                    post(nk,3) = sgn.*abs(Z.*cos(theta));
                    post(nk,4) = Z.*sin(theta);
                    
                    
                else 
                    
                    %For specular condition
                    
                    post(nk,3) = -post(nk,3);
                    
                end
                
                
            else
                
                
                post(nk,2) = Y_updated;
                if(~spec_boxes(bottle_neck))
                    
                    sgn = -sign(post(nk,4));
                    Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
                    theta = rand()*2*pi;
                    
                    post(nk,3) = Z.*cos(theta);
                    post(nk,4) = sgn.*abs(Z.*sin(theta));
                    
                else 
                    
                    %For speuclar condition
                    
                    post(nk,4) = -post(nk,4);
                    
                end
            end
            

            bottle_neck = box_no(post(nk,1:2), boxes);
            
            
        end
        
    end
    
     nk = rand(electron_population, 1) < P_scat;           % Random dsitribution of particles within the scattering probability limit
     post(nk,3:4) = random(Distr_MB, [sum(nk),2]);          %Scattering particles using the exponential scattering probability 
    
    
    
    %Calculating and Recording the Temperature of the electrons
    Temp(ni) = (sum(post(:,3).^2) + sum(post(:,4).^2))*C_m/K/2/electron_population;           
    
   
    
    for nk=1:electron_num
        
        %Storing the positions of the electrons on the trajector 
        Trj(ni, (2*nk):(2*nk+1)) = post(nk, 1:2);
        
    end
    
    
    
    
    %Plotting the electron trajectories
    % The animation updates for every 10 iterations
    if(animation_plot && mod(ni,10) == 0)
        figure(7);
        hold off;
        plot(post(1:electron_num,1), post(1:electron_num,2), 'o');       
        hold on;
        
        
          for nk=1:size(boxes,1)          %Plotting the rectangular boxes 
            
           plot([boxes(nk, 1) boxes(nk, 1) boxes(nk, 2) boxes(nk, 2) boxes(nk, 1)],...
               [boxes(nk, 3) boxes(nk, 4) boxes(nk, 4) boxes(nk, 3) boxes(nk, 3)], 'k-');
          
          end
        
        title('Electron Trajectories, SM101092039'); 
        xlabel('X position)');
        ylabel('Y position');
        pause(0.05);
        
    end
    
end


%Final Plotting of the boxes after completion of iterations
for nk=1:size(boxes,1)
    
   plot([boxes(nk, 1) boxes(nk, 1) boxes(nk, 2) boxes(nk, 2) boxes(nk, 1)],...
       [boxes(nk, 3) boxes(nk, 4) boxes(nk, 4) boxes(nk, 3) boxes(nk, 3)], 'k-');
   
end

%Final plotting of the trajectories after all the iterations have been
%completed

figure(7);
title('Electron Trajectories, SM101092039');
xlabel('X position');
ylabel('Y position');

grid on;
hold on;

%Storing the trajectory after completion 
for ni=1:electron_num
    
    plot(Trj(:,ni*2), Trj(:,ni*2+1), '.');
    
end





%% Part 3 Combining the Simulators to extract simple paramters

% a) 
% Plot the Density Map 

%Utilizing the Hist3 command to generate the density map 
elec_conc = hist3(post(:,1:2),[200 100])';

 %Utilzing the bins from the historgram to plot the the Temperature and Density map 
hist_bins = 10;    

[x_elec y_elec] = meshgrid(round(-hist_bins/2):round(hist_bins/2), round(-hist_bins/2):round(hist_bins/2));
M=exp(-x_elec.^2/(2*1^2)-y_elec.^2/(2*1^2));
M=M./sum(M(:));

figure(8);
elec_conc = conv2(elec_conc,M,'same');       %2D convolution for generating a 2D map;
elec_conc = elec_conc/(Y_R./size(elec_conc,1)*X_R./size(elec_conc,2));
surf(conv2(elec_conc,M,'same'));
title('Electron Density @ 0.8V , SM101092039');
xlabel('X position');
ylabel('Y position');

% Setting the Length and Width paramters of the rectangular region
%Utilizng the Length and width to configure the elements ranges
nx = 100*L;
ny = 100*W;




for  num_Iter = 1:5
    
    
      C = zeros(ny,nx); % Map of Conductivity 
      
      
      
     W_BN = 0.4*(1+(num_Iter/20));      % Adjusting the Width of the bottle neck boxes with the iteration factor
     L_BN = 0.4*(1+(num_Iter/20));      % Adjusting the Length of the bottle neck boxes with the iteration factor
     
     

    for h = 1:ny
        
        % Adjusting the sigma value for the elements based on the location 
         for g = 1:nx
           
            
            if(g >= nx*W_BN && g <= nx-nx*W_BN && (h >= ny-ny*L_BN || h <= ny*L_BN))               % If Inside the box then sigma = 10e-2
                C(h,g) = sig_2;
                
            
            else
                
                % If outside the box, sigma = 1 
                C(h,g) = sig_1;
                
                
            end
            
        end
        
    end
    
    
    
    % Initialize and create the G matrix
    G = sparse(nx*ny,nx*ny);
    
    
    %Initialize and create G Matrix 
    F = zeros(nx*ny,1);
    

    
    for g = 1:nx
        
        for h = 1:ny

                                           
            n = h + (g - 1)*ny;              % Node mapping equation
            
            
            

            % Local Mapping of Nodes around g and h 
            
            nxm = h + (g - 2)*ny;
            nxp = h + g*ny;
            nym = (h - 1) + (g - 1)*ny;
            nyp = (h + 1) + (g - 1)*ny;
            
            

           
            if(g == 1 || g == nx) 
                G(n,n) = 1;
                
                
            %Evaluation at the bottom region 
            
            elseif (h == 1)   
                
                U_Y = (C(h,g)+C(h+1,g))/2;
                U_X = (C(h,g)+C(h,g+1))/2;
                UX_D = (C(h,g)+C(h,g-1))/2;

                G(n,n) = -(U_Y + U_X + UX_D);
                G(n,nyp) = U_Y;
                G(n,nxp) = U_X;
                G(n,nxm) = UX_D;

            
                
            %Evaluation at the top region 
            
            elseif (h == ny) 

                D_Y = (C(h,g)+C(h-1,g))/2;
                U_X = (C(h,g)+C(h,g+1))/2;
                UX_D = (C(h,g)+C(h,g-1))/2;

                G(n,n) = -(D_Y + U_X + UX_D);
                G(n,nym) = D_Y;
                G(n,nxp) = U_X;
                G(n,nxm) = UX_D;
                
                
           
            %The finite difference method is being applied to evaluate the
            %regions potetntials outside as well.      
                
            else  
                
                U_Y = (C(h,g)+C(h+1,g))/2;
                D_Y = (C(h,g)+C(h-1,g))/2;
                U_X = (C(h,g)+C(h,g+1))/2;
                UX_D = (C(h,g)+C(h,g-1))/2;

                G(n,n) = -(U_Y + D_Y + U_X + UX_D);
                G(n,nyp) = U_Y;
                G(n,nym) = D_Y;
                G(n,nxp) = U_X;
                G(n,nxm) = UX_D;
                

            end
            
        end
        
    end
    
    
   
    for g = 1:nx
        
        for h = 1:ny
            
          
            n = h + (g - 1)*ny;   %Node Mapping Equation 

          
            
            if (g == 1) 
                
                %Indicating a shift towards left so the value must be set equal to V_0
                F(n) = V_0;
                
                
            end

           
            
        end
        
    end
    
    

    %Utilizing the relationship GV = F to solve for equation
    V = G\F;

    
    %Initliazing Vmap to 0 for setting up the matrix
    
    
    Vmap = 0;
    
    for g = 1:nx
        
        for h = 1:ny
            
            % Node mapping equation
            n = h + (g - 1)*ny;

            
            %Calculating the Gradient of V
            
            
            Vmap(h,g) = V(n);
            
        end
        
        
    end
    
   

    [E_x,E_y] = gradient(-Vmap);     %The electric field can be calculated from the gradient same as the Assignment 2
    
   
    
    % Calculating the current flow in the different locations and
    % directions and examining the variation with respect to the boxes.
   
    
    box_SL(num_Iter) = sum(C(:,1).*E_x(:,1));
    box_SR(num_Iter) = sum(C(:,nx).*E_x(:,nx));
    
    
end


% c) 

% Plotting the Variation of Current with respect to the changes in the
% dimension of the boxes
figure(9)
plot(linspace(1,5,5),box_SR,'b:')
hold on
plot(linspace(1,5,5),box_SR,'b:')
set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;  
a_x.GridColor = [1, 1, 1]; 


hold off
title({'Plot of Current v Changing Box Width, SM101092039'})
xlabel('Box Size Changes')
ylabel({'I (A)'})
lg.TextColor = 'white'; 




%%
function box_num = box_no(pos, boxes)   


%%function for determining the position of the position of the boxes and plotting

    box_num = 0;  %initializing at 0 
    
    for i=1:size(boxes,1)
        
        if(pos(1) > boxes(i,1) && pos(1) < boxes(i,2) && pos(2) > boxes(i,3) && pos(2) < boxes(i,4))
            
            box_num = i;
            
            return;
            
        end
        
    end
    
end








