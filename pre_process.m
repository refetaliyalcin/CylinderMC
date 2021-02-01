disp('Calculating Mie Scattering...');
n_cdf_random=1000;
nang_gid = 500; % number of outgoing far field angles to evaluate
nang_gel = nang_gid/2; % number of incident field angles to evaluate, it is within the coating and is not same with theta
ang_d=transpose(linspace(0,180,nang_gid));
ang_r=pi*ang_d/180;
cos_ang_r=cos(ang_r);


gelen_acilar=transpose(linspace(0,90,nang_gel));
gelen_acilar(1)=10^-5;

cdf2=zeros(nang_gid,1);
Area_fiber=radius*radius*pi;
rng shuffle
inv_cdf=zeros(length(gelen_acilar),n_cdf_random,length(lamda));
% pdf=zeros(length(gelen_acilar),nang_gid);
alfa=zeros(length(gelen_acilar),length(lamda));
beta=zeros(length(gelen_acilar),length(lamda));

nm = 1;         % outer medium refractive index (real) which is air
x2=cos(transpose(linspace(0,pi,nang_gid)));
for i=1:length(lamda)    
%     i*100/length(lamda)
    ns = sio2_n(lamda(i)) + sio2_k(lamda(i))*1j;     % cylinder refractive index (complex)
    alfa_i=zeros(length(gelen_acilar),1);
    beta_i=zeros(length(gelen_acilar),1);
    inv_cdf_i=zeros(length(gelen_acilar),n_cdf_random);
    lamda_i=lamda(i);
    parfor k=1:length(gelen_acilar)
        cdf2=zeros(nang_gid,1);
        [S_NP, C] = calccyl(radius, ns, nm,lamda_i, nang_gid, gelen_acilar(k));
        y2=S_NP./trapz(x2,S_NP);
        for i2=2:nang_gid
            cdf2(i2)=trapz(x2(1:i2),y2(1:i2));
        end
        alfa_i(k)=single(mean(C.sca)*f_v/Area_fiber);
        beta_i(k)=single(mean(C.abs)*f_v/Area_fiber);
        inv_cdf_i(k,:)=single(interp1(cdf2,ang_r,linspace(0,1,n_cdf_random),'spline','extrap'));
    end
    alfa(:,i)=alfa_i; %scattering coef
    beta(:,i)=beta_i; %absorption coef
    inv_cdf(:,:,i)=inv_cdf_i; %incerse of cdf
end

% zero_is_good=sum(isnan(inv_cdf(:)))
mu_tot=beta+alfa; %extinction coeff
sca_prob=alfa./mu_tot; %scattering albedo

