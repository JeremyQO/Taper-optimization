%%
% clear all; close all; clc
calculatemodes = true; % decide to enable or disable (for speed) calculation
if ~calculatemodes % if calculation is disabled the data should be previously stored in a mat file
    load tutorial3ls.mat modesMonerieCore modesTsaoTClad modesTwoClad modesTwoCore modesErdoganClad
end

%materials = {1.4615, 1.458, 'air'}; % 1.4574, 1.4525 for SM 800
materials = {1.45625, 1.45282, 'air'}; % 1.4574, 1.4525 for SM 800 at 852nm
%materials = {1.4574, 1.4525, 'air'}; % Takao
fibre.materials = materials;
fibre.coreCladdingRatio = 4.8/125; % 3.3/125 for SM 600 , 4.8/125 for SM 800

Lambda = 780.24;
argument = struct(...
    'type', 'dia',... % calculate mode curve n_eff vs. fibre diameter
    'harmonic', 1,... % required 
    'min', 0.01,... % minimum diameter
    'max', 600); % maximum diameter
figure();
%% Monerie modes
% Let's now calculate the Monerie modes in this fibre. In a separate demo
% on Monerie modes, I show that they should be traced from the core. We are
% interested in the LP01 mode corresponding to the HE11 mode as well as in
% LP11 corresponding to TE01 and TM01 modes. Our fibre core is single-mode,
% the LP11 mode is not supported. In order to trace it, we artificially
% increase the maximum considered diameter so that this mode is found and
% traced. Note: if argument.max is set to 180, there is a mode jump as
% explained in the Monerie tutorial
task = struct(...
    'nu', [0 1 2],... % first modal index
    'type', {{'monerie'}},... % mode types
    'maxmode', Inf,... % how many modes of each type and NU to calculate
    'lambda', Lambda,...
    'region', 'core');
argument.max = 125;
if calculatemodes
    modesMonerieCore = buildModes(argument, fibre, task, false);
end
hMonerieCore = showModes(modesMonerieCore);
set(hMonerieCore, 'Color', [1 0.5 0.5], 'LineWidth', 1)
% set(hMonerieCore(2), 'LineStyle', '--')
% uistack(hTwoCore, 'top')
% uistack(hTwoClad, 'top')
xlim([0 125])
% Which modes have been found?
t = modeDescription(modesMonerieCore, false, false);
fprintf('Monerie (nu=[0 1 2]) modes found in the cladding: ');
if iscell(t) 
    for i = 1:numel(t)-1, fprintf(' %s,', t{i}), end, fprintf(' %s\n', t{end});
else
    fprintf(' %s\n', t);
end
%% Erdogan mode
% The three-layer fundamental mode (HE11) can be calculated using the
% solultion of ERDOGAN (1997, 2000).
task = struct(...
    'nu', 1,... % first modal index
    'type', {{'erdogan'}},... % mode types
    'maxmode', 3,... 
    'lambda', Lambda,...
    'region', 'cladding');
if calculatemodes
modesErdoganClad = buildModes(argument, fibre, task, false);
end
% modesErdoganClad = addPointsToMode(modesErdoganClad, [0 3]); 
% modesErdoganClad = addPointsToMode(modesErdoganClad, [0 3]); 
hErdoganClad = showModes(modesErdoganClad);
set(hErdoganClad, 'Color', 'cyan', 'LineWidth', 1)
set(hErdoganClad, 'Visible', 'off') % Make TE curve dashed
set(hErdoganClad(3), 'Visible', 'on', 'Color', 'blue')
%%
% HE11 mode with erdogan modes
neff_HE11 = modesMonerieCore.NEFF;
D_HE11 = modesMonerieCore.ARG;

neff_HE12 = modesErdoganClad(3).NEFF;
% neff_HE11 = modesTwoClad(1).NEFF;

D_HE12 = modesErdoganClad(3).ARG;
% D_HE11 = modesTwoClad(1).ARG;

B_HE11 = 2*pi*neff_HE11/(Lambda*1e-9);
B_HE12 = 2*pi*neff_HE12/(Lambda*1e-9);

B_HE11_interp = interp1(D_HE11,B_HE11, D_HE12);

Omega = D_HE12*1e-6/2/2/pi.*(B_HE11_interp-B_HE12); % 4.4 is a non-anderstand paramater i have to put in order to match my results with those of Ryutaro Nagai and Takao Aoki.
% I corrected this line to be (2pi)^-1 as expected. A second 1/2 to go from
% diameter to radius


figure();
plot(D_HE12,neff_HE12); hold on
plot(D_HE11,neff_HE11);

% figure();
% plot(D_HE12,B_HE12); hold on
% plot(D_HE12,B_HE11_interp);

figure();
plot(D_HE12/2, Omega); grid on

% %% HE11 mode using 2-layer theory
% neff_HE11_ = modesMonerieCore.NEFF;
% D_HE11_ = modesMonerieCore.ARG;
% 
% neff_HE12_ = modesTwoClad(3).NEFF;
% % neff_HE11 = modesTwoClad(1).NEFF;
% 
% D_HE12_ = modesTwoClad(3).ARG;
% % D_HE11 = modesTwoClad(1).ARG;
% 
% B_HE11_ = 2*pi*neff_HE11_/600e-9;
% B_HE12_ = 2*pi*neff_HE12_/600e-9;
% 
% B_HE11_interp_ = interp1(D_HE11_,B_HE11_, D_HE12_);
% 
% Omega_ = D_HE12_*1e-6/pi.*(B_HE11_interp_-B_HE12_);
% 
% figure();
% plot(D_HE12_,neff_HE12_); hold on
% plot(D_HE11_,neff_HE11_);
% 
% % figure();
% % plot(D_HE12,B_HE12); hold on
% % plot(D_HE12,B_HE11_interp);
% 
% figure();
% plot(D_HE12_/2, Omega_);