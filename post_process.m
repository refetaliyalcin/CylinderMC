close all
clc
lamda_ub=max(lamda)*10^6;
T_sur=300;
T_amb=300;

polar_angle_rad=polar_angle*pi/180;
Solar_l=I_solar(lamda)';

trans_atm_pure=emis_atm_new(lamda)';
emit_atm=1-trans_atm_pure.^(1./cos(polar_angle_rad));

BB_Tamb_l=I_bb(lamda,T_amb)';
BB_Tsur_l=I_bb(lamda,T_sur)';
emittance_rad=trapz(lamda,BB_Tsur_l.*emittance,1);
emittance_atm=trapz(lamda,BB_Tamb_l.*emittance.*emit_atm,1);

P_rad=trapz(polar_angle_rad,2*cos(polar_angle_rad).*sin(polar_angle_rad).*emittance_rad)
P_atm=trapz(polar_angle_rad,2*cos(polar_angle_rad).*sin(polar_angle_rad).*emittance_atm)
P_sol=trapz(lamda,Solar_l.*absorptance(:,1))
P_net=P_rad-P_atm- P_sol





% 
% V = 2900*10^-9;
% A = repmat(lamda,[1 length(V)]);
% [minValue,closestIndex] = min(abs(A-V'));
% 
% figure(4)
% plot(lamda*10^6,Solar_l/max(Solar_l),lamda*10^6,BB_Tsur_l/max(BB_Tsur_l),lamda(closestIndex:end)*10^6,trans_atm_pure(closestIndex:end),lamda*10^6,emittance(:,1),'LineWidth',2)
% xlim([0 lamda_ub])
% xlabel('Wavelength [\mum]')
% legend('Normalized G_S_o_l_a_r','Normalized E_b_b','\tau_a_t_m','\epsilon_c_o_a_t_i_n_g')
% 
% 
% figure(5)
% semilogx(lamda*10^6,Solar_l/max(Solar_l),lamda*10^6,BB_Tsur_l/max(BB_Tsur_l),lamda(closestIndex:end)*10^6,trans_atm_pure(closestIndex:end),lamda*10^6,emittance(:,1),'LineWidth',2)
% xlim([0 lamda_ub])
% xlabel('Wavelength [\mum]')
% legend('Normalized G_S_o_l_a_r','Normalized E_b_b','\tau_a_t_m','\epsilon_c_o_a_t_i_n_g')
% 
% 
% [minValue,closestIndex_min] = min(abs(lamda-8*10^-6)); %find lamda array number for 8 micrometer
% [minValue,closestIndex_max] = min(abs(lamda-13*10^-6));%find lamda array number for 13 micrometer
% 
% 
% emit_2=squeeze(emittance(closestIndex_min:closestIndex_max,:));
% emit_2=mean(emit_2,1);
% figure(8)
% ax1 = polaraxes;
% plot_polar_angle_rad=[-flip(polar_angle_rad) polar_angle_rad];
% plot_emit_2=[flip(emit_2) emit_2];
% polarplot(ax1,plot_polar_angle_rad,plot_emit_2,'LineWidth',2)
% ax1.ThetaDir = 'clockwise';
% ax1.ThetaZeroLocation = 'top';
% rlim([0 1.05])
% thetalim([-90 90])
% title('Average Directional Emissivity at \ between 8-13\mum ')
% 
% e_lamda_nobinder_cylinder_r_200_h_50=emittance(:,1);
% e_teta_nobinder_cylinder_r_200_h_50=emit_2;