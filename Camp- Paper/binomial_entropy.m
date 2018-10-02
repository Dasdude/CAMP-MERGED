clc
close all
clear
n=10;
p = 0:.01:1
ent = .5.*(log2(2.*pi.*exp(1).*n.*p.*(1-p)));
plot(p,ent);
var = n.*p.*(1-p);
figure;plot(p,var)