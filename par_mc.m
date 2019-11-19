function [reflectance] = par_mc(mu_tot_wl,sca_prob_wl,inv_cdf_wl,cos_gelen_tara,photon_number,nang_gel,h,n_cdf_random,wl,cdf_gamma)
    reflectance=zeros(1,length(cos_gelen_tara));
    for i=1:length(cos_gelen_tara) %loop through incoming angle
        ref=0;
%         abso=0; %only consider reflectance for simplicity as if it is not reflected it is absorbed.
%         tra=0;
        cos_gelen_ev=cos_gelen_tara(i);
        for s=1:photon_number %loop throug photons
            r_no=0;
%             a_no=0;
%             t_no=0;

            x=0;
            y=0;
            z=0;

            s_z = cos_gelen_ev;
            sintheta = sqrt(1.0 - s_z*s_z);
            psi = 2.0*pi*rand();
            cospsi = cos(psi);
            if (psi < pi)
                sinpsi = sqrt(1.0 - cospsi*cospsi); 
            else
                sinpsi = -sqrt(1.0 - cospsi*cospsi);
            end
            s_x = sintheta*cospsi;
            s_y = sintheta*sinpsi;

            alive=1;

            while alive   
                
                costheta_f= cosd(cdf_gamma(ceil(rand()*1000))); % start to determine fiber orientation
                sintheta_f = sqrt(1.0 - costheta_f*costheta_f);
                psi_f = 2.0*pi*rand();
                cospsi_f = cos(psi_f);
                if (psi_f < pi)
                    sinpsi_f = sqrt(1.0 - cospsi_f*cospsi_f); 
                else
                    sinpsi_f = -sqrt(1.0 - cospsi_f*cospsi_f);
                end
                v_x = sintheta_f*cospsi_f;
                v_y = sintheta_f*sinpsi_f;
                v_z = costheta_f; % end of  fiber orientation determination


                zeta_ind=ceil(nang_gel*acosd(abs(s_x*v_x+s_y*v_y+s_z*v_z))/90); %index of angle between fiber and incident photon bundle
                if zeta_ind<1 || imag(zeta_ind)>0
                    zeta_ind=1;
                end

                l_beta=-log(rand())/mu_tot_wl(zeta_ind); %ext length

                if (s_z>0)
                    l_w = (h - z)/s_z; %distance to lower boundary
                else
                    l_w = -z/s_z; %distance to upper boundary
                end

                if l_w<l_beta
                    min_index=1;
                    min_l=l_w;
                else
                    min_index=2;
                    min_l=l_beta;
                end

                if (min_index==1)
                    alive=0;
                    if s_z>0
%                         a_no=1; %black background
                    else
                        r_no=1;
                    end
                else
                    if rand()<sca_prob_wl(zeta_ind)
            %               disp('scattering');

                        x=x+min_l*s_x;
                        y=y+min_l*s_y;
                        z=z+min_l*s_z;
                        
                        zeta_ind=ceil(nang_gel*acosd(abs(s_x*v_x+s_y*v_y+s_z*v_z))/90);
                        if zeta_ind<1 || imag(zeta_ind)>0
                            zeta_ind=1;
                        end
                        if rand()>0.5 
                            a=-1; 
                        else
                            a=1;
                        end
                        theta_s=a*inv_cdf_wl(zeta_ind,ceil(n_cdf_random*rand())); %find angle to scatter 
                        
                        %start of rodrigues' rotation
                        cos_theta_s=cos(theta_s);
                        if theta_s>0
                            sin_theta_s=sqrt(1-cos_theta_s*cos_theta_s);
                        else
                            sin_theta_s=-sqrt(1-cos_theta_s*cos_theta_s);
                        end
                        crosskv_1 = v_y*s_z - v_z*s_y;
                        crosskv_2 = v_z*s_x - v_x*s_z;
                        crosskv_3 = v_x*s_y - v_y*s_x;

                        dot_k_v = v_x*s_x+v_y*s_y+v_z*s_z;
                        s_x = cos_theta_s*s_x + crosskv_1*sin_theta_s + v_x*dot_k_v*(1 - cos_theta_s);
                        s_y = cos_theta_s*s_y + crosskv_2*sin_theta_s + v_y*dot_k_v*(1 - cos_theta_s);
                        s_z = cos_theta_s*s_z + crosskv_3*sin_theta_s + v_z*dot_k_v*(1 - cos_theta_s);
                        %end of rodrigues' rotation and new direction is determined

                    else
            %               disp('absorption');
                        alive=0;
%                         a_no=1;
                    end
                end
            end
            if r_no==1
                ref=ref+1;
            end
        end
        reflectance(i)=ref/photon_number;
    end
end