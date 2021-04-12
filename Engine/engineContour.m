


function [divE, V, Np, Nn, Vrp, Vrn, Vzp, Vzn, Er, Ez] = engineContour( nr, nz,lr,lz )
%%engineContour reads in results from ion-ion engine simulation
%Collects data into matrix form and creates contour and quiver plots
%Reads from data file named "results.txt"
%Inputs:
%   nr: number of segments in r direction
%   nz: number of segments in z direction
%   lr: max radius length
%   lz: max axial length
%Outputs:
%   V: Voltage in space with filled contour plot
%   Np: Positive Ion denisity with filled contour plot
%   Nn: Negative Ion denisity with filled contour plot
%   Vzp, Vrp: Positive ions' axial and radial velocities and quiver plot
%   Vzn, Vrn: Negative ions' axial and radial velocities and quiver plot
%   Er, Ez: Electric field in r and z direction, quiver plot added to
%                  voltage contour plot
%   conservedMassP: conservation of mass check for positive ions
%                   -This should be 0 at equilib.
%   conservedMassN: conservation of mass check for negative ions
%                   -This should be 0 at equilib.
%
%Format for all plots:
%   x-axis: axial length
%   y-axis: radial length
%
%Format for all tables:
%   (1,1) corresponds to the left, most inner corner of thruster domain

%% Start data collecting and plotting

NR=nr+1;
NZ=nz+1;

dr = lr/nr;
dz = lz/nz;

Z = 0:dz:lz-dz;
R = 0:dr:lr-dr;

Zvolt = 0:dz:lz;
Rvolt = 0:dr:lr;
for i=1:NR
    V(i, :) = dlmread('results.txt', ' ', [0, 1+NZ*(i-1), 0, NZ*i]);
end


for i=1:nr
    Np(i, :) = dlmread('results.txt', ' ', [3, 1+nz*(i-1), 3, nz*i]);
    Nn(i, :) = dlmread('results.txt', ' ', [4, 1+nz*(i-1), 4, nz*i]);
    Vrp(i, :) = dlmread('results.txt', ' ', [5, 1+nz*(i-1), 5, nz*i]);
    Vrn(i, :) = dlmread('results.txt', ' ', [6, 1+nz*(i-1), 6, nz*i]);
    Vzp(i, :) = dlmread('results.txt', ' ', [7, 1+nz*(i-1), 7, nz*i]);
    Vzn(i, :) = dlmread('results.txt', ' ', [8, 1+nz*(i-1), 8, nz*i]);
    Er(i, :) = dlmread('results.txt', ' ', [9, 1+nz*(i-1), 9, nz*i]);
    Ez(i, :) = dlmread('results.txt', ' ', [10, 1+nz*(i-1), 10, nz*i]);
end

Zdiff = dz:dz:lz-dz;
Rdiff = dr:dr:lr-dr;
for i=1:nr-1
    for j = 1:nz-1
        divE(i,j) = ((Er(i+1,j)-Er(i,j))/dr)+(Er(i,j)/(i*dr))+((Ez(i,j+1)-Ez(i,j))/dz);
    end
end
%divE = divergence(Z, R, Ez, Er);

figure
% Voltage and E-Field plot
contourf(Zvolt, Rvolt,V, 50, 'LineCOlor', 'none');
colormap(jet)
%caxis([-7500 7500])
shading flat
hold on
quiver(Z,R,Ez,Er,'w')
title('Voltage and E-Fields')
figure

%E-Field Divergence w/ constant r
contourf(Zdiff,Rdiff,divE, 100, 'LineColor', 'none');
title('Divergence of E-Field')
figure

% Positive Ion Density Plot
contourf(Z,R,Np, 50, 'LineColor', 'none');
title('Positive Density')
figure

% Negative Ion Density Plot
contourf(Z,R,Nn, 50, 'LineColor', 'none');
title('Negative Density')
figure

%Density Total
contourf(Z,R,Np+Nn, 50, 'LineColor', 'none');
title('Density Total (Np+Nn)')
figure

contourf(Z,R,Np-Nn, 50, 'LineColor', 'none');
title('Density Difference (Np-Nn)')
figure

quiver(Z,R,Vzp, Vrp, 'r')
title('Velocity')
hold on
quiver(Z,R,Vzn, Vrn, 'k')
legend('Positive Ions', 'Negative Ions')
axis([0 lz 0 lr])

%% Start Conservation Checks

%Positive Ion Mass Check
leftNp=0;
rightNp=0;
topNp=0;

for i=1:nr
    leftNp = leftNp + Np(i, 1)*Vzp(i,1)*i*dr*dr;
    rightNp = rightNp + Np(i, nz)*Vzp(i, nz)*i*dr*dr;
end
for j=1:nz
    topNp = topNp + Np(nr, j)*Vrp(nr,j)*lr*dz;
end
NetFluxP = 2*pi*(leftNp+topNp+rightNp)

%Negative Ion Mass Check
leftNn=0;
rightNn=0;
topNn=0;

for i=1:nr
    leftNn = leftNn + Nn(i, 1)*Vzn(i,1)*i*dr*dr;
    rightNn = rightNn + Nn(i, nz)*Vzn(i, nz)*i*dr*dr;
end
for j=1:nz
    topNn = topNn + Nn(nr, j)*Vrn(nr,j)*lr*dz;
end
NetFluxN = 2*pi*(leftNn+topNn+rightNn)


end