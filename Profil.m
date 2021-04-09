function [sectiondata] =Profil

sectiondata.sectionType = 'Idoublysymmetric';
sectiondata.b = 180e-3; % [m]
sectiondata.houtside=(400)*10^(-3); % [m]
sectiondata.tf = 13.5e-3; % [m]
sectiondata.h = sectiondata.houtside-sectiondata.tf; % [m]
sectiondata.hinside=sectiondata.houtside-2*sectiondata.tf; % [m]
sectiondata.tw = 8.6e-3; % [m]
sectiondata.material.type = 'isotropic';
sectiondata.material.E = 210e9; %[N/m^2]
sectiondata.material.nu = 0.3; %[-]

Iyy=(2*(sectiondata.b)*sectiondata.tf^3 + sectiondata.tw*(sectiondata.h-sectiondata.tf)^3)/12 + 2* (sectiondata.b*sectiondata.tf) *(sectiondata.h/2)^2;
IYY=((sectiondata.b*sectiondata.houtside^3) - (sectiondata.b-sectiondata.tw) *sectiondata.hinside^3)/12;
assert(abs(Iyy/IYY-1)<2e-15,'Iyy caluclations differ')
Area=2*sectiondata.b*sectiondata.tf+sectiondata.hinside*sectiondata.tw;
assert(Area>0,'Area negativ')
Izz=(2*sectiondata.tf*sectiondata.b^3+sectiondata.hinside*sectiondata.tf^3)/12;
assert(Izz>0,'Izz negativ')

end