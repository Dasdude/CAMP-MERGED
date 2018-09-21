clc
close all
clear
mu= 2;
omega = 3;
nak_samp = random('Nakagami',mu,omega,[100000,1]);
gam_samp = random('Gamma',mu,omega/mu,[100000,1]);
gam_samp = sqrt(gam_samp);
figure;
histogram(nak_samp)
title('Nakagami')

figure 
histogram(gam_samp);
title('Gamma');

