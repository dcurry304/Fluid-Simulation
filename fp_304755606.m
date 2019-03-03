%Final Project 
%David Curry
%ID: 304755606
clear all;
clc;
close all;
rng('default');

dt = 0.003;  %timestep
tf = 0.7;  %final time  / 1.5 for a
B = 0.40; %damping wall coefficeint
rho0 = 1500;  %initial density
n = 100;  %number of particles  /  120 for part a
m = rho0/n;   %mass of a single particle
k = 80;  %stifness constant  / 80 for part a
mu = 0.8;  %viscosity  / 0.8 for a
a = 5; %makes gravity force stronger

xmax = 4; %define size of domain x /4 for a
ymax = 4;  %define size of domain y
h = 0.1;  %define smoothing radius / 0.10 for a
Nx = floor(xmax/h);  %num x bins
Ny = floor(ymax/h);  %num y bins
dx = xmax/Nx;  %spacing of bins in x direction
dy = ymax/Ny;  %spacing of bins in y direction

%create particles structure and bins structure
particles(1:n) = struct('pos',[0,0],'vel',[0,0],'force',[0,0],'rho',[],'neigh',[]);
bins(1:Nx*Ny) = struct('partIDs',[],'adjBins',[]);

%part 1 values
xa = xmax/2;    xa2 = xmax;
xb = 0;         xb2 = xmax/2;
ya = ymax;
yb = 0;

%fill postion part of particle structure with random x and y values for
%part a
for k2 = 1:n
    particles(k2).pos = [(xa-xb)*rand + xb,(ya-yb)*rand+yb];
end


% %part 2: obstacles
% xf = 3.9;
% yf = 3;
% length = sqrt(yf^2 + xf^2);
% spacing = 0.03*h;
% num = floor(length/spacing);
% xo = zeros(1,num);
% yo = zeros(1,num);
% count = 0;
% slope = -1;
% for e = 1:num
%     xo(e) = 0 + count;
%     yo(e) = slope*xo(e) + yf;
%     count = count + spacing;
% end
% normSlope = (1/sqrt(1 + slope^2))*[1,-slope];
% 
% 
% xa = xmax/2;
% xb = 0;
% ya = ymax;
% yb = 3;
% for k2 = 1:n
%       particles(k2).pos = [(xa-xb)*rand + xb,(ya-yb)*rand+yb];
% end

%%part c water balloon
% xa1 = 0.5;   xa2 = 0.5;   xa3 = 0.5;
% xb1 = 0;     xb2 = 0;     xb3 = 0;
% ya1 = 3;     ya2 = 2;     ya3 = 4;
% yb1 = 3.3;   yb2 = 2.3;   yb3 = 3.7;
% for k2 = 1:n
%     if k2 <= n/3
%         particles(k2).pos = [(xa1-xb1)*rand + xb1,(ya1-yb1)*rand+yb1];
%     elseif k2 > 2*n/3
%         particles(k2).vel = [30,5];
%         particles(k2).pos = [(xa2-xb2)*rand + xb2,(ya2-yb2)*rand+yb2];
%
%     else
%         particles(k2).vel = [10,-10];
%         particles(k2).pos = [(xa3-xb3)*rand + xb3,(ya3-yb3)*rand+yb3];
%
%     end
% end

figure;
%save the video if true, just graph the video if false
savevid = true;
if savevid == true
    vid = VideoWriter('trial3', 'MPEG-4');
    vid.FrameRate = 30;
    vid.Quality = 80;
    open(vid)
end

%x and y coordinates of particles to be used later
x = zeros(1, n);
y = zeros(1, n);

x1 = [];
y1 = [];
A = ones(Ny,Nx);

%fill the adjacentbins part of bins with the bins adjacent to each k3 bin
%this doesnt change throughout
for k3 = 1: Nx*Ny
    if k3 == 1 %top left corner
        bins(k3).adjBins = [(k3 + 1), (k3 + Ny), (k3 +Ny +1)];
    elseif k3 < Ny  %left side
        bins(k3).adjBins = [(k3 - 1), (k3 + 1), (k3 + Ny - 1), (k3 + Ny), (k3 + Ny +1)];
    elseif k3 == Ny  %bottom left corner
        bins(k3).adjBins = [(k3 - 1), (k3 + Ny), (k3 + Ny - 1)];
    elseif k3 == ((Nx*Ny) - Ny + 1)  %top right corner
        bins(k3).adjBins = [(k3 - Ny), (k3 - Ny + 1), (k3 + 1)];
    elseif mod(k3,Ny) == 1 %top side
        bins(k3).adjBins = [(k3 - Ny), (k3 - Ny + 1), (k3 + 1), (k3 + Ny), (k3 + Ny + 1)];
    elseif k3 == Nx*Ny  %bottom right corner
        bins(k3).adjBins = [(k3 - Ny), (k3 - Ny - 1), (k3 - 1)];
    elseif mod(k3,Ny) == 0  %bottom side
        bins(k3).adjBins = [(k3 - Ny), (k3 - Ny - 1), (k3 - 1), (k3 + Ny - 1), (k3 + Ny)];
    elseif k3 > ((Nx*Ny) - Ny + 1)  %left side
        bins(k3).adjBins = [(k3 - Ny - 1), (k3 - Ny), (k3 - Ny + 1), (k3 - 1), (k3 + 1)];
    else                      %middle
        bins(k3).adjBins = [(k3 - Ny - 1), (k3 - Ny), (k3 - Ny + 1), (k3 - 1), (k3 + 1), (k3 + Ny - 1), (k3 + Ny), (k3 + Ny + 1)];
    end
end

binNum = zeros(1,n);

%loop thru all time loops
%stop scrolling it starts here
for t = 0:dt:tf
    %fill the A matrix with zeros
    A = zeros(Ny,Nx);
    %empty the particle Ids field of bins
    for z = 1:Ny*Nx
        bins(z).partIDs = [];
    end
    
    %put the positions of the particles in the correct bins
    %in the partIDs part of bins 
    for k5 = 1:Ny*Nx  %loop thru bins
        for k4 = 1:n  %loop thru particles
            binNum(k4) = (ceil(particles(k4).pos(1)/dx) - 1)*Ny + ceil((ymax - particles(k4).pos(2))/dy);
            if binNum(k4) == k5
                bins(k5).partIDs = [bins(k5).partIDs,k4];
            end
        end
    end
    %empty the neighbor particles field of particles
    for p2 = 1:n
        particles(p2).neigh = [];
    end
    
    %find neighbors of every initial particle by checking if any particles are
    %within h of that particle. Update the neigh bin of particles.
    for z = 1:Nx*Ny  %loop thru bins
        partBinz = bins(z).partIDs; %find particles in that bin
        if partBinz ~= 0 %check if bin is empty
            adjacent = (bins(z).adjBins);  %find adj bins
            for w = [z,adjacent]  %loop thru all adjacent bins
                partBinw = bins(w).partIDs; %find particles in w adjacent bin
                for k6 = partBinz %loop thru particles in z bin
                    for j = partBinw  %loop thru particles in adj bin w
                        dist = sqrt(((particles(k6).pos(1)) - (particles(j).pos(1)))^2 + ((particles(k6).pos(2)) - (particles(j).pos(2)))^2);
                        if (dist < h) && (k6 ~= j)  %if close enough they are neighbors
                            particles(k6).neigh = [(particles(k6).neigh), j];
                        end
                    end
                end
            end
        end
    end
    
    %find density for all particles and put it in the rho field of particles
    for k7 = 1:n  %loop thru particles
        neighPart = particles(k7).neigh;  %find neighbor particles
        sum = 0;
        for w2 = neighPart  %loop thru neighbors
            norm = sqrt(((particles(k7).pos(1)) - (particles(w2).pos(1)))^2 + ((particles(k7).pos(2)) - (particles(w2).pos(2)))^2);
            sum = sum + ((4*m)/(pi*h^8))*((h^2 - norm^2)^3);  %sum of contributions from neighbors
        end
        particles(k7).rho = ((4*m)/(pi*h^2)) + sum; %density formula
    end
    
    %find force on each particle and put in the force field of particles
    for k8 = 1:n  %loop thru all particles
        f_ext = [0,-9.8]*a*(particles(k8).rho);  %external force
        P_k = k*(particles(k8).rho - rho0);  %pressure of particle k
        f = [0,0];
        neighpart2 = particles(k8).neigh;
        for w3 = neighpart2  %loop thru neighbor particles
            P_j = k*(particles(w3).rho - rho0);  %pressure of neighbor particle
            q = (sqrt((particles(k8).pos(1)-particles(w3).pos(1))^2 + (particles(k8).pos(2)-particles(w3).pos(2))^2)/h);
            f_t = ((m/(pi*h^4*particles(w3).rho))*(1 - q))*((15*k*(particles(k8).rho + particles(w3).rho - 2*rho0)*((1-q)/q)*(particles(k8).pos - particles(w3).pos)) - (40*mu*(particles(k8).vel - particles(w3).vel)));
            f = f + f_t;   %total force from pressure and viscoscity
        end
        particles(k8).force = f + f_ext; %put force in force field of particles
    end
    
    %find new velocity and new positions based on semi explicit euler
    for k9 = 1:n
        particles(k9).vel = particles(k9).vel + (dt*(particles(k9).force))/particles(k9).rho;
        particles(k9).pos = particles(k9).pos + dt*(particles(k9).vel);
        %if the particles go outside the bounds, change the position and
        %velocity to make them in bounds
        if particles(k9).pos(1) >= xmax
            particles(k9).pos(1) = (2*xmax) - particles(k9).pos(1);
            particles(k9).vel(1) = (-B)*(particles(k9).vel(1));
        elseif particles(k9).pos(1) <= 0
            particles(k9).pos(1) =  0 - particles(k9).pos(1);
            particles(k9).vel(1) = (-B)*(particles(k9).vel(1));
        end
        if particles(k9).pos(2) >= ymax
            particles(k9).pos(2) = (2*ymax) - particles(k9).pos(2);
            particles(k9).vel(2) = (-B)*(particles(k9).vel(2));
        elseif particles(k9).pos(2) <= 0
            particles(k9).pos(2) =  0 - particles(k9).pos(2);
            particles(k9).vel(2) = (-B)*(particles(k9).vel(2));
        end
    end
%     %%part 2
%     for e = 1:num
%         for k =  1:n
%             if particles(k).pos(1) == xo(e) && particles(k).pos(2) == yo(e)
%             particles(k).vel = -particles(k).vel(1).*normSlope;
%             
%             end
%         end
%     end
    
    %put number of particles in bins into a matrix 
%     r = 1;
%     rr = 0;
%     while r < Nx*Ny
%         target = bins(r).partIDs;
%         c = floor(r/Ny)+1;
%         if r < Ny
%             rr = r;
%         elseif mod(r,Ny) == 1 && r > Ny
%             rr = 1;
%         end 
%         A(rr,c) = length(target);
%         r = r + 1 ;
%         rr = rr + 1;
%     end
%     y1 = zeros(1,Ny);
%     s1 = 1;
%     sx = 0;
%     while sx <= xmax
%         y1(s1) = sx;
%         sx = sx + dx;
%         s1 = s1 + 1;
%     end
%     x1 = zeros(1,Nx);
%     s2 = 1;
%     sy = 0;
%     while sy <= ymax
%         x1(s2) = sy;
%         sy = sy + dy;
%         s2 = s2 + 1;
%     end
%     contour(x1,y1,A,1);

   
    
    %plot the points every timestep
    for s = 1:n
        x(s) = particles(s).pos(1);
        y(s) = particles(s).pos(2);
        plot(x,y,'bo');
    end
    %put grid lines and max size of graph
    xlim([0 xmax]);
    ylim([0 ymax]);
    grid on
    

    %write the video 
    if savevid == true
        writeVideo(vid,getframe(gcf))
    else
        drawnow
    end
    
end
%close the video
if savevid == true
    close(vid)
end

