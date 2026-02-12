%% SSRL 10-1 gratings deesign 

%%
% awojdyla@lbl.gov
% December 2025


%% Beamline parameters

% % storage ring energy
% Esr_GeV = 2.0;
% % Storage ring current
% Isr_A = 500e-3;
% % energy spread
% dEpE = 0.98e-3;
% 
% 
% % undulator
% undulator_N = 55;
undulator_L_m = 2.071;
% Lambda0_m = 38e-3;

% operating range
Es_eV = linspace(50,1600,31);
lambdas_m = 1.2398e-06./Es_eV;

% vertical source-slit distance
pv_m = 18;
qv_m = 7;

% horizontal source-slit distance
ph_m = 27;
qh_m = 1;

%% Beam size and divergence
% 
% % ALS-U electron beam properties
% [Sx_m, Sxp_rad, Sy_m, Syp_rad] = Gen4.beam_size(Es_eV, undulator_L_m);
% 
% 
% plot(Es_eV, Sx_m*1e6, Es_eV, Sy_m*1e6)
% hold on
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(Es_eV, 12.33+Es_eV*0,'--')
% ax = gca;
% ax.ColorOrderIndex = 2;
% plot(Es_eV, 15.10+Es_eV*0,'--')
% hold off
% title('COSMIC-U undulator beam size')
% xlabel('photon energy [eV]')
% ylabel('beam size [\mum-rms]')
% legend('horizontal', 'vertical','electron H','electron V')
% 
% xlim([250 2500])
% set(gca,'xTick',0:250:2500)
% 
% ylim([0 30])
% set(gca,'yTick',0:5:50)
% 
% grid on
% 
% %%
% plot(Es_eV, Sxp_rad*1e6, Es_eV, Syp_rad*1e6)
% hold on
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(Es_eV, 6.16+Es_eV*0,'--')
% ax = gca;
% ax.ColorOrderIndex = 2;
% plot(Es_eV, 5.03+Es_eV*0,'--')
% hold off
% title('COSMIC-U undulator beam divergence')
% xlabel('photon energy [eV]')
% ylabel('beam divergence [\murad-rms]')
% legend('horizontal', 'vertical','electron H','electron V')
% 
% 
% xlim([250 2500])
% set(gca,'xTick',0:250:2500)
% 
% ylim([0 40])
% set(gca,'yTick',0:5:50)
% 
% grid on


%% VIA-VLS trajectories


g1_lpm = 300e3;
c1 = 3;
E1min_eV = 50;
E1opt_eV = 800;
E1max_eV = 1600;

g1_lpm = 400e3;
c1 = 2.5;
E1min_eV = 50;
E1opt_eV = 800;
E1max_eV = 1600;

g1_lpm = 400e3;
c1 = 1.5;
E1min_eV = 50;
E1opt_eV = 800;
E1max_eV = 1600;

mask_g1 = Es_eV>=E1min_eV & Es_eV<=E1max_eV;

[alpha1s_rad, beta1s_rad, theta1s_rad] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);

plot(Es_eV(mask_g1), alpha1s_rad(mask_g1)*180/pi,'r', Es_eV(mask_g1), -beta1s_rad(mask_g1)*180/pi,'b', Es_eV(mask_g1), theta1s_rad(mask_g1)*180/pi, 'k')
legend('entrance angle \alpha', 'exit angle |\beta|', 'half included angle \theta','location','Southeast')
xlabel('photon energy [eV]')
ylabel('incendence angle [deg]')
title('SSRL 10-1 300l/mm, cff=3 VIA-VLS trajectory')

xlim([0 1700])
ylim([85 90])
set(gca,'xTick',0:200:2000)
set(gca,'yTick',85:0.5:90)
grid on

%% Footprints

Syp_rad = 0.2e-3;
pv_m = 18;
fpM201_1a_m = sqrt((Syp_rad*pv_m).^2)./sin(pi/2-theta1s_rad);
fpG201a_m = sqrt((Syp_rad*pv_m).^2)./sin(pi/2-alpha1s_rad);



plot(Es_eV, fpM201_1a_m*1e3, Es_eV, fpG201a_m*1e3)
xlabel('photon energy [eV]')
ylabel('beam footprint [mm]')
title('vert ap 0.2 mrad SSRL 10-1 300l/mm, cff=3')
legend('grating','pre-mirror','location','SouthEast')
xlim([0 1600])
set(gca,'xTick',0:200:2000)
set(gca,'yTick',0:100:2000)
ylim([0 500])
grid



%%

%load('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/PDR/COSMICU/data/R2D_Rh_20keV.mat')
load('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/PDR/MAESTROU/data/R2D_Au_20keV_10deg.mat')
[E,T] = meshgrid(E_eV, ga_rad);
Reflec_Rh = @(En_eV, grazing_rad) interp2(E, T, R', En_eV, grazing_rad)';

theta1p_deg = 90- theta1s_rad*180/pi;

ga_rad = (0:0.01:5)*pi/180;
[E,G] = meshgrid(Es_eV, ga_rad);


clf
imagesc(Es_eV, ga_rad*180/pi, Reflec_Rh(E,G)')
colorbar
title('SSRL 10-1 Au pre-mirror reflectivity for 300l/mm, cff=3 VIA-VLS trajectory')
xlabel('photon energy [eV]')
ylabel('grazing angle [deg]')

hold on

h3 = plot(Es_eV(mask_g1),theta1p_deg(mask_g1));
set(h3, 'linewidth',2.25)

% h4 = plot(Es_eV(mask_g2), theta2p_deg(mask_g2));
% set(h4, 'linewidth',2.25)
% 
% h5 = plot(Es_eV(mask_g3), theta3p_deg(mask_g3));
% set(h5, 'linewidth',2.25)

h2 = plot(1./Es_eV*3.9*1e3, Es_eV,'--','color',[1 1 1]/2);

set(h2, 'linewidth',2.25)

set(gca,'YDir','normal')
xlabel('photon energy [eV]')
ylabel('grazing angle [deg]')

xlim([50 1600])
set(gca,'xTick',50:200:2000)

legend('pre-mirror angle')
hold off
colormap gray

%% SiC premirror
material = 'SiC';
data_folder = '/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/PDR/harmonic_contamination/data/';
load([data_folder,'SiC_2D.dat']);

E_eV = 200:10:8200;
ga_rad = (8:0.1:88)*1e-3;
Rr = SiC_2D(:,3);
R = reshape(Rr,[801,801])';
[E,T] = meshgrid(E_eV, ga_rad);
Reflec_SiC = @(En_eV, grazing_rad) interp2(E, T, R, En_eV, grazing_rad);
[E,G] = meshgrid(Es_eV, ga_rad);
imagesc(Es_eV, ga_rad*180/pi, Reflec_SiC(E,G))

clf
imagesc(Es_eV, ga_rad*180/pi, Reflec_SiC(E,G))
colorbar
title('SSRL 10-1 SiC pre-mirror reflectivity for 300l/mm, cff=3 VIA-VLS trajectory')
xlabel('photon energy [eV]')
ylabel('grazing angle [deg]')


hold on

h3 = plot(Es_eV(mask_g1),theta1p_deg(mask_g1));
set(h3, 'linewidth',2.25)

h2 = plot(1./Es_eV*3.9*1e3, Es_eV,'--','color',[1 1 1]/2);

set(h2, 'linewidth',2.25)

set(gca,'YDir','normal')
xlabel('photon energy [eV]')
ylabel('grazing angle [deg]')

xlim([50 1600])
set(gca,'xTick',50:200:2000)

legend('pre-mirror angle')
hold off
colormap gray


%% Lamellar gratings efficiency

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))

% lamellar
thickness1_m = 17e-9;
dc1 = 0.3;
material = 'Au';

Egs_eV = linspace(50,1600,31);
lambdags_m = 1.2398e-06./Egs_eV;
eta1gs = zeros(1,length(lambdags_m));

for i_l=1:length(lambdags_m)
    [alphas_rad,~,~] = Gen4.vls_trajectory(lambdags_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
     eta1gs(i_l) = Blazr.efficiency_lamellar_duty(1/g1_lpm, thickness1_m, ...
            lambdas_m(i_l), pi/2-alphas_rad(i_l), material,dc1);
end

%
eta_g1s = interp1(Egs_eV,eta1gs,Es_eV);

mask_g1 = Es_eV<=(E1max_eV*1.1);

clf
plot(Es_eV(mask_g1), eta_g1s(mask_g1))
title({'lamellar gratings diffraction efficiency',sprintf('g=%1.0fl/mm, cff=%1.2f, mat=%s, h=%1.0fnm, lgr=%1.0fp.c.', ...
    g1_lpm*1e-3,c1, material, thickness1_m*1e9,dc1*100)})
xlabel('photon energy [eV]')
ylabel('diffraction efficiency (+1 order)')

xlim([50 1600])
set(gca,'xTick',0:200:2000)

ylim([0 0.5])
set(gca,'yTick',0:0.1:1)
grid on

%% Lamellar harmonics

%% Lamellar gratings efficiency

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))

% lamellar
thickness1_m = 15e-9;
dc1 = 0.5;
material = 'SiC';

Egs_eV = linspace(50,1600,31);
lambdags_m = 1.2398e-06./Egs_eV;
eta1gs = zeros(1,length(lambdags_m));
eta2gs = zeros(1,length(lambdags_m));
eta3gs = zeros(1,length(lambdags_m));

for i_l=1:length(lambdags_m)
    [alphas_rad,~,~] = Gen4.vls_trajectory(lambdags_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
     eta1gs(i_l) = Blazr.efficiency_lamellar_higher(1/g1_lpm, thickness1_m, ...
            lambdas_m(i_l), pi/2-alphas_rad(i_l), material,dc1,1);
    eta2gs(i_l) = Blazr.efficiency_lamellar_higher(1/g1_lpm, thickness1_m, ...
    lambdas_m(i_l)/2, pi/2-alphas_rad(i_l), material,dc1,2);
       eta3gs(i_l) = Blazr.efficiency_lamellar_higher(1/g1_lpm, thickness1_m, ...
    lambdas_m(i_l)/3, pi/2-alphas_rad(i_l), material,dc1,3);
end

%
eta_g1s = interp1(Egs_eV,eta1gs,Es_eV);
eta_g2s = interp1(Egs_eV,eta2gs,Es_eV);
eta_g3s = interp1(Egs_eV,eta3gs,Es_eV);

mask_g1 = Es_eV<=(E1max_eV*1.1);

clf
plot(Es_eV(mask_g1), eta_g1s(mask_g1), ...
     Es_eV(mask_g1), eta_g2s(mask_g1), ...
     Es_eV(mask_g1), eta_g3s(mask_g1))
title({'lamellar gratings diffraction efficiency',sprintf('g=%1.0fl/mm, cff=%1.2f, mat=%s, h=%1.0fnm, lgr=%1.0fp.c.', ...
    g1_lpm*1e-3,c1, material, thickness1_m*1e9,dc1*100)})
xlabel('photon energy [eV]')
ylabel('diffraction efficiency (+1 order)')
legend('fundamental','second harmonic','third harmonic')

xlim([50 1600])
set(gca,'xTick',0:200:2000)

ylim([0 0.5])
set(gca,'yTick',0:0.1:1)
grid on


%% Blazed gratings efficiency

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))


% blazed angle
theta2_blaze_rad = 0.60*pi/180;
thickness2_m = 1/g1_lpm*tan(theta2_blaze_rad);
material = 'Au';


eta_blaze = zeros(1,length(lambdags_m));
for i_l=1:length(lambdags_m)
    [alphas_rad,~,~] = Gen4.vls_trajectory(lambdags_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
    eta_blaze(i_l) = Blazr.efficiency_blazed(1/g1_lpm, thickness2_m, ...
        lambdags_m(i_l), pi/2-alphas_rad(i_l), material);
end

eta_blaze_gs = interp1(Egs_eV,eta_blaze,Es_eV);

plot(Es_eV(mask_g1), eta_blaze_gs(mask_g1))
title({'blazed gratings diffraction efficiency',sprintf('g=%1.0fl/mm, cff=%1.2f, mat=%s, blaze=%1.2fdeg', ...
    g1_lpm*1e-3,c1,material, theta2_blaze_rad*180/pi)})
xlabel('photon energy [eV]')
ylabel('diffraction efficiency (+1 order)')

xlim([50 1600])
set(gca,'xTick',0:200:2000)

ylim([0 0.5])
set(gca,'yTick',0:0.1:1)
grid on


%% harmonics

% blazed angle
theta2_blaze_rad = 0.4*pi/180;
thickness2_m = 1/g1_lpm*tan(theta2_blaze_rad);
material = 'SiC';

eta_blaze = zeros(1,length(lambdags_m));
eta_harm2 = zeros(1,length(lambdags_m));
eta_harm3 = zeros(1,length(lambdags_m));
for i_l=1:length(lambdags_m)
    [alphas_rad,~,~] = Gen4.vls_trajectory(lambdags_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
    eta_blaze(i_l) = Blazr.efficiency_blazed(1/g1_lpm, thickness2_m, ...
        lambdags_m(i_l), pi/2-alphas_rad(i_l), material);
    eta_harm2(i_l) = Blazr.efficiency_blazed_higher(1/g1_lpm, thickness2_m, ...
        lambdags_m(i_l)/2, pi/2-alphas_rad(i_l), material,2);
    eta_harm3(i_l) = Blazr.efficiency_blazed_higher(1/g1_lpm, thickness2_m, ...
        lambdags_m(i_l)/3, pi/2-alphas_rad(i_l), material,3);
end

eta_blaze_gs = interp1(Egs_eV,eta_blaze,Es_eV);
eta_harm2_gs = interp1(Egs_eV,eta_harm2,Es_eV);
eta_harm3_gs = interp1(Egs_eV,eta_harm3,Es_eV);

%
plot(Es_eV(mask_g1), eta_blaze_gs(mask_g1), ...
     Es_eV(mask_g1), eta_harm2_gs(mask_g1), ...
     Es_eV(mask_g1), eta_harm3_gs(mask_g1))
title({'blazed gratings diffraction efficiency',sprintf('g=%1.0fl/mm, cff=%1.2f, mat=%s, blaze=%1.2fdeg', ...
    g1_lpm*1e-3,c1,material, theta2_blaze_rad*180/pi)})
xlabel('photon energy [eV]')
ylabel('diffraction efficiency (+1 order)')
legend('fundamental','second harmonic','third harmonic')

xlim([50 1600])
set(gca,'xTick',0:200:2000)

ylim([0 0.5])
set(gca,'yTick',0:0.1:1)
grid on





%% c=0.1 VIA-VLS trajectories


g1_lpm = 3000e3;
c1 = 0.1;
E1min_eV = 50;
E1opt_eV = 800;
E1max_eV = 1600;

mask_g1 = Es_eV>=E1min_eV & Es_eV<=E1max_eV;

[alpha1s_rad, beta1s_rad, theta1s_rad] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
alphap1s_rad = sqrt(2*lambdas_m*g1_lpm/(1-c1^2));
alpha1s_rad = pi/2-alphap1s_rad;

plot(Es_eV(mask_g1), alpha1s_rad(mask_g1)*180/pi,'r', Es_eV(mask_g1), -beta1s_rad(mask_g1)*180/pi,'b', Es_eV(mask_g1), theta1s_rad(mask_g1)*180/pi, 'k')
legend('entrance angle \alpha', 'exit angle |\beta|', 'half included angle \theta','location','Southeast')
xlabel('photon energy [eV]')
ylabel('incendence angle [deg]')
title('SSRL 10-1 300l/mm, cff=3 VIA-VLS trajectory')

xlim([0 1700])
ylim([80 90])
set(gca,'xTick',0:200:2000)
set(gca,'yTick',85:0.5:90)
grid on

%%

%load('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/PDR/COSMICU/data/R2D_Rh_20keV.mat')
load('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/PDR/MAESTROU/data/R2D_Au_20keV_10deg.mat')
[E,T] = meshgrid(E_eV, ga_rad);
Reflec_Rh = @(En_eV, grazing_rad) interp2(E, T, R', En_eV, grazing_rad)';

theta1p_deg = 90- theta1s_rad*180/pi;

ga_rad = (0:0.01:5)*pi/180;
[E,G] = meshgrid(Es_eV, ga_rad);

clf
imagesc(Es_eV, ga_rad*180/pi, Reflec_Rh(E,G)')
colorbar
title('SSRL 10-1 Au pre-mirror reflectivity for 300l/mm, cff=3 VIA-VLS trajectory')
xlabel('photon energy [eV]')
ylabel('grazing angle [deg]')

hold on

h3 = plot(Es_eV(mask_g1),theta1p_deg(mask_g1));
set(h3, 'linewidth',2.25)

% h4 = plot(Es_eV(mask_g2), theta2p_deg(mask_g2));
% set(h4, 'linewidth',2.25)
% 
% h5 = plot(Es_eV(mask_g3), theta3p_deg(mask_g3));
% set(h5, 'linewidth',2.25)

h2 = plot(1./Es_eV*3.9*1e3, Es_eV,'--','color',[1 1 1]/2);

set(h2, 'linewidth',2.25)

set(gca,'YDir','normal')
xlabel('photon energy [eV]')
ylabel('grazing angle [deg]')

xlim([50 1600])
set(gca,'xTick',50:200:2000)

legend('pre-mirror angle')
hold off
colormap gray

%% Lamellar gratings efficiency

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))

% lamellar
thickness1_m = 7e-9;
dc1 = 0.4;
material = 'Au';

Egs_eV = linspace(50,1600,31);
lambdags_m = 1.2398e-06./Egs_eV;
eta1gs = zeros(1,length(lambdags_m));

for i_l=1:length(lambdags_m)
     eta1gs(i_l) = Blazr.efficiency_lamellar_higher(1/g1_lpm, thickness1_m, ...
            lambdas_m(i_l), pi/2-alpha1s_rad(i_l), material,dc1,-1);
end

%
eta_g1s = interp1(Egs_eV,eta1gs,Es_eV);

mask_g1 = Es_eV<=(E1max_eV*1.1);

clf
plot(Es_eV(mask_g1), eta_g1s(mask_g1))
title({'lamellar gratings diffraction efficiency',sprintf('g=%1.0fl/mm, cff=%1.2f, mat=%s, h=%1.0fnm, lgr=%1.0fp.c.', ...
    g1_lpm*1e-3,c1, material, thickness1_m*1e9,dc1*100)})
xlabel('photon energy [eV]')
ylabel('diffraction efficiency (-1 order)')

xlim([50 1600])
set(gca,'xTick',0:200:2000)

ylim([0 0.1])
set(gca,'yTick',0:0.01:0.1)
grid on



%% Blazed gratings efficiency

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))


% blazed angle
theta2_blaze_rad = 2.0*pi/180;
thickness2_m = 1/g1_lpm*tan(theta2_blaze_rad);
material = 'Au';


eta_blaze = zeros(1,length(lambdags_m));
for i_l=1:length(lambdags_m)
    eta_blaze(i_l) = Blazr.efficiency_blazed_higher(1/g1_lpm, thickness2_m, ...
        lambdags_m(i_l), pi/2-alphas_rad(i_l), material,-1);
end

eta_blaze_gs = interp1(Egs_eV,eta_blaze,Es_eV);

plot(Es_eV(mask_g1), eta_blaze_gs(mask_g1))
title({'blazed gratings diffraction efficiency',sprintf('g=%1.0fl/mm, cff=%1.2f, mat=%s, blaze=%1.2fdeg', ...
    g1_lpm*1e-3,c1,material, theta2_blaze_rad*180/pi)})
xlabel('photon energy [eV]')
ylabel('diffraction efficiency (-1 order)')

xlim([50 1600])
set(gca,'xTick',0:200:2000)

ylim([0 0.5])
set(gca,'yTick',0:0.1:1)
grid on

%%
eta_g1s = interp1(Egs_eV,eta1gs,Es_eV);
eta_blaze_gs = interp1(Egs_eV,eta_blaze,Es_eV);
eta_g3s = interp1(Egs_eV,eta3gs,Es_eV);

mask_g1 = Es_eV<=(E1max_eV*1.1);
mask_g2 = Es_eV>=(E2min_eV*0.9) & Es_eV<=(E2max_eV*1.1);
mask_g3 = Es_eV>=(E3min_eV*0.9);

clf
plot(Es_eV(mask_g1), eta_g1s(mask_g1), ...
     Es_eV(mask_g2), eta_blaze_gs(mask_g2), ...
     Es_eV(mask_g3), eta_g3s(mask_g3))

xlabel('photon energy [eV]')
ylabel('diffraction efficiency (+1 order) \eta_1')

xlim([250 2500])
set(gca,'xTick',0:250:2500)

ylim([0 0.5])
set(gca,'yTick',0:0.05:0.5)
grid on

hold on
ax = gca;

ax.ColorOrderIndex = 1;
plot(Es_eV, eta_g1s,'--')

ax.ColorOrderIndex = 2;
plot(Es_eV, eta_blaze_gs,'--')

ax.ColorOrderIndex = 3;
plot(Es_eV, eta_g3s,'--')

title('COSMIC-U gratings diffraction efficiency')
legend('G101 (lam., gold)','G102 (blazed, gold)','G103 (blazed, rhodium)',...
    'location','South')



%% Blaze scanning


RP =10000;
L_und = 2;
p_m = 18;
cff = linspace(1.2, 3, 11);
g_lpm = RP^2/p_m^2*L_und./(cff.^2-1)*(2*2.35/pi)^2;

Es_eV = linspace(50,1600,16);

E0_eV = 500;

lambda0_m = 1239e-9/E0_eV; 
% entrance angle (G1)
alpha0p_rad = sqrt(2*g_lpm*lambda0_m./(cff.^2-1));
% exit angle
beta0p_rad = cff.*alpha0p_rad;
% half included angle (M2 angle)
theta0_rad = (alpha0p_rad+beta0p_rad)/2;
beta_rad = -(alpha0p_rad-beta0p_rad)/2;

plot(cff, theta0_rad*180/pi)

plot(cff, beta_rad*180/pi)
%ylim([0 10])

%%

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))


E0_eV = 200;
RP =15000;

pv_m = 18;
qv_m = 7;

g_lpm = RP^2/p_m^2*L_und./(cff.^2-1)*(2*2.35/pi)^2;

Egs_eV = linspace(50,1600,31);
Egs_eV = linspace(50,300,31);
lambdas_m = 1.2398e-06./Egs_eV;

efficiency_map = zeros(length(lambdags_m), length(cff));
for i_c=1:length(cff)

    lambda0_m = 1239e-9/E0_eV; 
    % entrance angle (G1)
    alpha0p_rad = sqrt(2*g_lpm(i_c)*lambda0_m./(cff(i_c).^2-1));
    % exit angle
    beta0p_rad = cff(i_c).*alpha0p_rad;
    % half included angle (M2 angle)
    theta0_rad = (alpha0p_rad+beta0p_rad)/2;
    beta_rad = -(alpha0p_rad-beta0p_rad)/2;

    % blazed angle
    thickness_m = 1/g_lpm(i_c)*tan(beta_rad);
    material = 'Au';

    [alphas_rad, ~, ~] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E0_eV, g_lpm(i_c), pv_m, qv_m, cff(i_c));


    etas_blaze = zeros(1,length(lambdas_m));
    for i_l=1:length(lambdas_m)
        etas_blaze(i_l) = Blazr.efficiency_blazed(1/g_lpm(i_c), thickness_m, ...
            lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
    end
    
    efficiency_map(:,i_c) = etas_blaze;
    fprintf("%i",i_c)
end

%%
imagesc(Egs_eV,cff,efficiency_map')
xlabel("photon energy [eV]")
ylabel("cff")
colorbar
title(sprintf("efficiency of blazed (gold) for RP=%1.0f optimized at %1.0f eV", RP, E0_eV))

%%


RP =10000;
L_und = 2;
p_m = 18;
cff = linspace(1.2, 3, 11);
g_lpm = RP^2/p_m^2*L_und./(cff.^2-1)*(2*2.35/pi)^2;

Es_eV = linspace(50,1600,16);

E0_eV = 500;

lambda0_m = 1239e-9/E0_eV; 
% entrance angle (G1)
alpha0p_rad = sqrt(2*g_lpm*lambda0_m./(cff.^2-1));
% exit angle
beta0p_rad = cff.*alpha0p_rad;
% half included angle (M2 angle)
theta0_rad = (alpha0p_rad+beta0p_rad)/2;
beta_rad = -(alpha0p_rad-beta0p_rad)/2;

plot(cff, theta0_rad*180/pi)

plot(cff, beta_rad*180/pi)
%ylim([0 10])

%% focus on blazed

addpath(genpath('C:\Users\antoi\Documents\MATLAB\alsu-matlab\CD2_ESD\grating_efficiency'))
addpath(genpath('/Users/awojdyla/Documents/MATLAB/alsu-matlab/CD2_ESD/grating_efficiency'))


E0_eV = 150;
RP = 15000;

pv_m = 18;
qv_m = 7;
cff = 2.0;

g_lpm = RP^2/p_m^2*L_und./(cff.^2-1)*(2*2.35/pi)^2;

Egs_eV = linspace(50,1600,156);
%Egs_eV = linspace(50,1600,31);

lambdas_m = 1.2398e-06./Egs_eV;

lambda0_m = 1239e-9/E0_eV; 
% entrance angle (G1)
alpha0p_rad = sqrt(2*g_lpm*lambda0_m./(cff.^2-1));
% exit angle
beta0p_rad = cff.*alpha0p_rad;
% half included angle (M2 angle)
theta0_rad = (alpha0p_rad+beta0p_rad)/2;
beta_rad = -(alpha0p_rad-beta0p_rad)/2;

% blazed angle
thickness_m = 1/g_lpm*tan(beta_rad);
material = 'Au';

[alphas_rad, ~, thetas_rad] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E0_eV, g_lpm, pv_m, qv_m, cff);


etas_blaze = zeros(1,length(lambdas_m));
for i_l=1:length(lambdas_m)
    etas_blaze(i_l) = Blazr.efficiency_blazed(1/g_lpm, thickness_m, ...
        lambdas_m(i_l), pi/2-alphas_rad(i_l), material);
end

plot(Egs_eV, etas_blaze)


%save('HEG_a.mat','Es_eV', 'thetas_rad', 'etas_blaze', 'g_lpm', 'cff')
save('LEG_a.mat','Egs_eV', 'thetas_rad', 'etas_blaze', 'g_lpm', 'cff')

%% Pre-mirror reflectivity

[~,~,thetas_rad] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
R1_M102  = Reflec_Rh(Es_eV, pi/2-thetas_rad).';
[~,~,thetas_rad] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E2opt_eV, g2_lpm, pv_m, qv_m, c2);
R2_M102  = Reflec_Rh(Es_eV, pi/2-thetas_rad).';
[~,~,thetas_rad] = Gen4.vls_trajectory(lambdas_m, 1239e-9./E3opt_eV, g3_lpm, pv_m, qv_m, c3);
R3_M102  = Reflec_Rh(Es_eV, pi/2-thetas_rad).';

plot(Es_eV, R1_M102, Es_eV, R2_M102, Es_eV, R3_M102)
xlabel('photon energy [eV]')
ylabel('reflectivity R')
title('COSMIC-U pre-mirror reflectivity')
legend('G101','G102','G103','location','SouthEast')
xlim([250 2500])
ylim([0 1])
set(gca,'xTick',0:250:2500)
set(gca,'yTick',0:0.1:1)
grid on

%% Mirror reflectivity
R_mirror  = Reflec_Rh(Es_eV, 1.25*pi/180).';

plot(Es_eV, R_mirror,'k')
xlabel('photon energy [eV]')
ylabel('reflectivity R')
%title('COSMIC-U Si stripe reflectivity')
title('COSMIC-U mirror reflectivity (Rh, 1.25deg)')

xlim([250 2500])
set(gca,'xTick',0:250:2500)

ylim([0 1])
set(gca,'yTick',0:0.1:1)

grid on
hold off


%% Total light efficiency
eta1s = eta_g1s.*R1_M102.*R_mirror.^2;
eta2s = eta_blaze_gs.*R2_M102.*R_mirror.^2;
eta3s = eta_g3s.*R3_M102.*R_mirror.^2;

mask_g1 = Es_eV>=E1min_eV & Es_eV<=E1max_eV*1.1;
mask_g2 = Es_eV>=E2min_eV*0.9 & Es_eV<=E2max_eV*1.1;
mask_g3 = Es_eV>=E3min_eV*0.9 & Es_eV<=E3max_eV;

clf
plot(Es_eV, eta1s, Es_eV, eta2s, Es_eV, eta3s)
plot(Es_eV(mask_g1), eta1s(mask_g1), Es_eV(mask_g2), eta2s(mask_g2), Es_eV(mask_g3), eta3s(mask_g3))

xlim([0 2600])
ylim([0 0.25])
set(gca,'xTick',0:250:2500)
set(gca,'yTick',0:0.025:0.25)
grid on
xlabel('photon energy [eV]')
ylabel('light efficiency')
title('COSMIC-U overall light efficiency')
legend('G101','G102','G103','orientation','horizontal')

    
%% Resolution

S1y_m = 15e-6;
S1y_m = 80e-6./sqrt(Es_eV/100);
pv_m = 18;
qv_m = 8; % kind of arbitrary
E1opt_eV = 500;

g1_lpm = 400e3;
c1 = 2.8;

Es_eV = linspace(50,1600,31);
lambdags_m = 1.2398e-06./Egs_eV;
[alpha1s_rad,beta1s_rad,theta1s_rad] = Gen4.vls_trajectory(lambdags_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
dlambda_source1_m = S1y_m./pv_m.*cos(alpha1s_rad)./g1_lpm;
dlambda_slits1_m  = S1y_m.*(qv_m/(c1*pv_m))./qv_m.*cos(beta1s_rad)./g1_lpm;


RP1 = 1239e-9./Es_eV./sqrt(dlambda_source1_m.^2+dlambda_slits1_m.^2)/2.35;

plot(Es_eV, RP1)

title('SSRL 10-1 HEG')
xlabel('photon energy [eV]')
ylabel('resolving power E/\DeltaE (FHWM)')
%legend('10-1HEG','orientation','vertical','location','SouthWest')


set(gca,'xTick',000:200:2000)
grid on
xlim([0 1600])
ylim([0 20000])



%% Footprints
% footprint on grating

S1yp_rad = 200e-6;
S1yp_rad = 800e-6./sqrt(Es_eV/100);

fpM102_m = sqrt(S1y_m.^2+(S1yp_rad*pv_m).^2)./sin(pi/2-theta1s_rad);
fpG10x_m = sqrt(S1y_m.^2+(S1yp_rad*pv_m).^2)./sin(pi/2-alpha1s_rad);
plot(Es_eV, fpM102_m*1e3, Es_eV, fpG10x_m*1e3)

title('SSRL 10-1 HEG Monochromator footprints')
xlabel('photon energy [eV]')
ylabel('beam footprint [mm-6\sigma]')

grid on
xlim([0 1600])
ylim([0 300])
set(gca,'yTick',0:20:300)
legend('M2','HEG','location','SouthEast')
grid on


%% Resolution

S1y_m = 80e-6;
S1y_m = 80e-6./sqrt(Es_eV/100);
pv_m = 18;
qv_m = 8; % kind of arbitrary
E1opt_eV = 500;

g1_lpm = 1000e3;
c1 = 2.0;

Es_eV = linspace(50,1600,31);
lambdags_m = 1.2398e-06./Egs_eV;
[alpha1s_rad,beta1s_rad,theta1s_rad] = Gen4.vls_trajectory(lambdags_m, 1239e-9./E1opt_eV, g1_lpm, pv_m, qv_m, c1);
dlambda_source1_m = S1y_m./pv_m.*cos(alpha1s_rad)./g1_lpm;
dlambda_slits1_m  = S1y_m.*(qv_m/(c1*pv_m))./qv_m.*cos(beta1s_rad)./g1_lpm;


RP1 = 1239e-9./Es_eV./sqrt(dlambda_source1_m.^2+dlambda_slits1_m.^2)/2.35;

plot(Es_eV, RP1)

title('SSRL 10-1 LEG')
xlabel('photon energy [eV]')
ylabel('resolving power E/\DeltaE (FHWM)')
%legend('10-1LEG','orientation','vertical','location','SouthWest')


set(gca,'xTick',000:200:2000)
grid on
xlim([0 1600])
ylim([0 20000])



%% Footprints
% footprint on grating

S1yp_rad = 800e-6;
S1yp_rad = 800e-6./sqrt(Es_eV/100);

fpM102_m = sqrt(S1y_m.^2+(S1yp_rad*pv_m).^2)./sin(pi/2-theta1s_rad);
fpG10x_m = sqrt(S1y_m.^2+(S1yp_rad*pv_m).^2)./sin(pi/2-alpha1s_rad);
plot(Es_eV, fpM102_m*1e3, Es_eV, fpG10x_m*1e3)

title('SSRL 10-1 LEG Monochromator footprints')
xlabel('photon energy [eV]')
ylabel('beam footprint [mm-6\sigma]')

grid on
xlim([0 1600])
ylim([0 300])
set(gca,'yTick',0:20:300)
legend('M2','LEG','location','SouthEast')
grid on

