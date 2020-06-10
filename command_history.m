clear; %clear work space and all figures.
clc;
close all;

load shock.mat;

plot(RVS); % Plot the acceleration history of shock 'RVS'. Here only 'RVS' is
           % demonstrated as an example. Other three shocks are also provided.

M=RVS.srs(100,10); % Calculate the resposne matrix 'M' from 100Hz, with damping Q=10.

M.plot; % Plot the SRS of 'RVS'.

M.contourf; % Plot shock response contour, i.e. SRC.

N1=M.svd(1);  % SVD decomposition of SRC. This command will return the matrix N1.
M.svd(1);     % The dual spectrum and margin l, along with the weighted
              % singular evctor in time domain plot will show up here.
           
[Mx, Nx]=M.mdof(MI); % This shows how to calculate structural resposne with
[~,N1x]=N1.mdof(MI);     % response matrix. Here 'MI' is modal information, as 
                     % shown in Table 1 in the paper.
                     
Mx.plot; % These commands show the time history plot of Mx, Nx and N1x.
Nx.plot;
N1x.plot;


