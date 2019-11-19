# CylinderMC
Solution of Radiative Transfer Equation in a fibrous medium with Monte Carlo method

Works with MATLAB

Thickness of coating, volume fraction and radius of fibers can be modified from main.m

Execute main.m to run the code, if the parallel toolbox is available code will be run on parallel

calccyl.m is missing for now, I am still waiting to take the proper permission from Jan Schäfer, however it can be downloaded from the link below:
https://www.mathworks.com/matlabcentral/fileexchange/36831-matscat

However the codes in asmcyl.m should change from present form;

C.sca(1) = 4/k*(abs(bnp(1)).^2 + 2*sum(abs(bnp(2:end)).^2));

C.sca(2) = 4/k*(abs(ann(1)).^2 + 2*sum(abs(ann(2:end)).^2));

to this:

C.sca(1) = 4/k*(abs(bnp(1)).^2 + 2*sum(abs(bnp(2:end)).^2+abs(anp(2:end)).^2));

C.sca(2) = 4/k*(abs(ann(1)).^2 + 2*sum(abs(ann(2:end)).^2+abs(bnn(2:end)).^2));

You can contact me from refetali@gmail.com
