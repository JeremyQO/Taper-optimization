%% Optical Fibre Toolbox two-layer demo
% This script demonstrates the use of Optical Fibre Toolbox functions for
% calculation of two-layer fibres modes.
% 
% (cc-by) K. Karapetyan et al., AG Meschede, Uni Bonn, 2008--2011
% 
% kotya.karapetyan@gmail.com, http://agmeschede.iap.uni-bonn.de

%%
clc
clear
close all

tStart = tic;
addpath('..')
%% Specify the fibre parameters

% Note: diameter is specified in micrometers, wavelength in nanometers

% Fibre materials (core, cladding)
materials = {1.46, 'air'};

% Fibre structure description
fibre = struct(...
    'materials', {materials});

%% Field distribution
poiType = 'wvl';
poi = 635; %nm
par = 0.3; %diamètre fibre µm
task = struct('modetype', 'HYBRID', 'modeindex', [1 1]);
infomode = false;
getWholeMode = false;
n = neff(poiType, poi, par, fibre, task, infomode, getWholeMode);

% 
d = par;
window = 3 * d; 
lambda = poi;
Nr = [500, 1000]; % inn, out
Nphi = [64, 64]; % inn, out
dr_inn = d/2 / Nr(1);
dr_out = (window/2 - d/2) / Nr(2);

dphi_inn = 2 * pi / Nphi(1);
dphi_out = 2 * pi / Nphi(2);

%%
F_circ = struct('dr', [dr_inn dr_out], 'dphi', [dphi_inn dphi_out], ...
    'diam', [d window], 'E1', [], 'H1', [], 'E2', [], 'H2', []);
FG = fieldGrid(F_circ);
% [F.E1, F.H1] = modeField(lambda, d, n, fibre, task, FG.R1, FG.PHI1);
% [F.E2, F.H2] = modeField(lambda, d, n, fibre, task, FG.R2, FG.PHI2);
[F_circ.E1, F_circ.H1] = modeFieldCirc(lambda, d, n, fibre, task, FG.R1, FG.PHI1);
[F_circ.E2, F_circ.H2] = modeFieldCirc(lambda, d, n, fibre, task, FG.R2, FG.PHI2);
F_circ.FG = FG;

displayField2(F_circ, d)
P_circ = fieldPower(F_circ);
%%
F = struct('dr', [dr_inn dr_out], 'dphi', [dphi_inn dphi_out], ...
    'diam', [d window], 'E1', [], 'H1', [], 'E2', [], 'H2', []);
FG = fieldGrid(F);
[F.E1, F.H1] = modeField(lambda, d, n, fibre, task, FG.R1, FG.PHI1);
[F.E2, F.H2] = modeField(lambda, d, n, fibre, task, FG.R2, FG.PHI2);
F.FG = FG;

displayField2(F, d)
P = fieldPower(F);
%%
% Convert to Cartesian coordinates
[X1,Y1] = pol2cart(F.FG.PHI1,F.FG.R1);
[X2,Y2] = pol2cart(F.FG.PHI2,F.FG.R2);

EX1 = F.E1(:,:,1) .* cos(F.FG.PHI1) - F.E1(:,:,2) .* sin(F.FG.PHI1);
EY1 = F.E1(:,:,1) .* sin(F.FG.PHI1) + F.E1(:,:,2) .* cos(F.FG.PHI1);
EZ1 = F.E1(:,:,3);
EX2 = F.E2(:,:,1) .* cos(F.FG.PHI2) - F.E2(:,:,2) .* sin(F.FG.PHI2); 
EY2 = F.E2(:,:,1) .* sin(F.FG.PHI2) + F.E2(:,:,2) .* cos(F.FG.PHI2); 
EZ2 = F.E2(:,:,3);

%%
% Convert to Cartesian coordinates
[X1_circ,Y1_circ] = pol2cart(F_circ.FG.PHI1,F_circ.FG.R1);
[X2_circ,Y2_circ] = pol2cart(F_circ.FG.PHI2,F_circ.FG.R2);

EX1_circ = F_circ.E1(:,:,1) .* cos(F_circ.FG.PHI1) - F_circ.E1(:,:,2) .* sin(F_circ.FG.PHI1);
EY1_circ = F_circ.E1(:,:,1) .* sin(F_circ.FG.PHI1) + F_circ.E1(:,:,2) .* cos(F_circ.FG.PHI1);
EZ1_circ = F_circ.E1(:,:,3);
EX2_circ = F_circ.E2(:,:,1) .* cos(F_circ.FG.PHI2) - F_circ.E2(:,:,2) .* sin(F_circ.FG.PHI2); 
EY2_circ = F_circ.E2(:,:,1) .* sin(F_circ.FG.PHI2) + F_circ.E2(:,:,2) .* cos(F_circ.FG.PHI2); 
EZ2_circ = F_circ.E2(:,:,3);
%%
max(max(sqrt(EX2.^2+EY2.^2)))
max(max(sqrt(EX2_circ.^2+EY2_circ.^2)))
max(max(sqrt(EX2.^2+EY2.^2)))
min(min(sqrt(EX2_circ.^2+EY2_circ.^2)))
max(max(abs(F.E2))).^2
max(max(abs(F_circ.E2))).^2

figure;
hold on
surf( F.E2,'EdgeAlpha',0,'FaceAlpha',1);
surf( X2, Y2, EX2,'EdgeAlpha',0,'FaceAlpha',1);
theta = 0:.1:2*pi;
polar(theta,ones(size(theta)) * d/2)
title('E_x (a.u.)');xlabel ('X (\mum)'); ylabel('Y(\mum)'); axis on; view([0 90]);
colorbar;
%% Mode dispersion in a nanofibre
% Calculate mode dispersion diagram for a fixed wavelength and varying
% fibre diameter
%
% Specify the fibre parameters

materials = {'silica'; 'air'};

nanofibre = struct(...
    'materials', {materials});

argument = struct(...
    'type', 'dia',... % calculate vs. wavelength
    'min', 0.1,... % calculate from
    'max', 1); % calculate to

modeTask = struct(...
    'nu', [0 1 2],... % first modal index
    'type', {{'hybrid', 'te', 'tm'}},... % mode types
    'maxmode', 3,... % how many modes of each type and NU to calculate
    'lambda', 900);%,... % parameter, structure diameter, if argument is wavelength
    %'region', 'cladding');
    
modes = buildModes(argument, nanofibre, modeTask, false);
%%
% Display calculated dispersion curves
figure;
showModes(modes, 'Modal dispersion in a nanofibre');
%%
% Show which modes have been calculated
fprintf('Modes found:\n');
for i = 1:numel(modes)
   fprintf('%s\n', modeDescription(modes(i), false));
end
%%
% Unlike in the weakly guidance example above, this nanofibre has a high
% refractive index step (about 0.5) and is therefore strongly guiding.
% All four modes HE11, HE21, TE01 and TM01 are therefore clearly distinct.
% The LP approximation is not valid in this case any more. The Optical
% Fibre Toolbox correctly simulates this situation using the full vector
% solution of the Maxwell equations.

toc(tStart)