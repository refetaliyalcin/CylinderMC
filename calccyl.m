function [T1, C] = calccyl( r, ns, nm, lambda, nang, zeta)

	%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
	%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
	%   Organization: Institut für Lasertechnologien in der Medizin und
	%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)
	k = 2*pi/lambda*nm;     % wavenumber in outer medium nm
	x = k*r;                % size parameter
	m = ns/nm;              % relative refractive index

	zeta=zeta*pi/180;
	xsin = x*sin(zeta);
	M = ceil((xsin + 4*(xsin^(1/3)) + 2));
	n = 0:M;

	%% Calculate auxiliary variables
	xi = x*sin(zeta);
	eta = x*sqrt(m^2 - cos(zeta)^2);
	jneta = besselj(n, eta);
	djneta = 0.5*(besselj(n-1,eta) - besselj(n+1,eta));
	jnxi = besselj(n, xi);
	djnxi = 0.5*(besselj(n-1,xi) - besselj(n+1,xi));
	hnxi = besselh(n, 1, xi);
	dhnxi = 0.5*(besselh(n-1,1,xi) - besselh(n+1,1,xi));

	An = 1.j*xi*(xi*djneta.*jnxi - eta*jneta.*djnxi);
	Bn = xi*(m^2*xi*djneta.*jnxi - eta*jneta.*djnxi);
	Cn = n*cos(zeta)*eta.*jneta.*jnxi*(xi^2/eta^2 - 1);
	Dn = n*cos(zeta)*eta.*jneta.*hnxi*(xi^2/eta^2 - 1);
	Vn = xi*(m^2*xi*djneta.*hnxi - eta*jneta.*dhnxi);
	Wn = 1.j*xi*(eta*jneta.*dhnxi - xi*djneta.*hnxi);

	wvd = (Wn.*Vn + 1.j*Dn.^2);
	cd = 1.j*Cn.*Dn;

	idx = isnan(wvd);

	%% Calculate expansion coefficients
	anp = (Cn.*Vn - Bn.*Dn)./wvd;
	ann = -(An.*Vn - cd)./wvd;
	bnp = (Wn.*Bn + cd)./wvd;
	bnn = -1.j*(Cn.*Wn + An.*Dn)./wvd;

	anp(idx) = 0.;
	ann(idx) = 0.;
	bnp(idx) = 0.;
	bnn(idx) = 0.;

	n = 1:(numel(anp)-1);
	ang = (0:nang-1)/(nang-1)*pi;

	% Calculate amplitude scattering matrix
	T = zeros(nang,2,2);
	for iang=1:nang
		T(iang,1,1) = bnp(1) + 2*sum(bnp(2:end).*cos(n*ang(iang)));
		T(iang,2,2) = ann(1) + 2*sum(ann(2:end).*cos(n*ang(iang)));
		T(iang,2,1) = -2*1.j*sum(anp(2:end).*sin(n*ang(iang)));
		T(iang,1,2) = -2*1.j*sum(bnn(2:end).*sin(n*ang(iang)));
	end %for iang=1:nang

	% Calculate cross sections
	C.ext(1) = 4/k*real(T(1,1,1));
	C.ext(2) = 4/k*real(T(1,2,2));
	C.sca(1) = 4/k*(abs(bnp(1)).^2 + 2*sum(abs(bnp(2:end)).^2+abs(anp(2:end)).^2));
	C.sca(2) = 4/k*(abs(ann(1)).^2 + 2*sum(abs(ann(2:end)).^2+abs(bnn(2:end)).^2));

	C.abs(1) = C.ext(1) - C.sca(1);
	C.abs(2) = C.ext(2) - C.sca(2);

	T1 = 0.5*(abs(T(:,1,1)).^2 + abs(T(:,1,2)).^2 + abs(T(:,2,1)).^2 + abs(T(:,2,2)).^2);
end