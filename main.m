clc
clear all;
close all;

photon_number=10000;
f_v=0.04; %volume fraction 0.04 = 4%
lamda=(300:100:25000)*10^-9; %wavelength in meter
polar_angle=linspace(0,89.99999,30); %theta from 0 to 90 deg
sigma=20; %sigma theta_f in deg
radius=250*10^-9; % radius in meter
pre_process % calculate cyliner scatterg parameters
h=10*10^-3; %thickness of coating in meter
cdf_gamma = transpose(icdf(makedist('Normal','mu',90,'sigma',sigma),linspace(0.000000000001,0.999999999999999,1000)));
reflectance=zeros(length(lamda),length(polar_angle));
cos_gelen_tara=cosd(polar_angle);

rng shuffle
disp('Calculating Monte Carlo method...');
tic
parfor idx = 1 : length(lamda)
    reflectance(idx,:)=par_mc(mu_tot(:,idx),sca_prob(:,idx),inv_cdf(:,:,idx),cos_gelen_tara,photon_number,nang_gel,h,n_cdf_random,lamda(idx),cdf_gamma);
end
toc
absorptance=1-reflectance;
emittance=absorptance;
disp('Monte Carlo ended');
post_process

