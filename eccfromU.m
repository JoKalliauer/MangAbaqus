function [sectiondata,yecc,zecc]=eccfromU(UBpUges)


 % A=80.678; %cm^2
 % Iyy=21876.5; %cm^4
 % e=0:1:300; %cm
 % ratio=1./(1+I./(A*e.^2)); %-
 % plot(e,ratio)
 %
 %
 % myratio=.5; %-
 % mye=sqrt(I/(A*(1/myratio-1))); %cm
 
assert(UBpUges>=0,'Energyratio must be larger than 0')
assert(UBpUges<=1,'Energyratio must be smaler then 1')
%1-eps(10000000)=0.999999998137355=1-1-eps(10000000)Izz
%UMpUges=min(UMpUges,1-eps(10000000));
UBpUges=min(UBpUges,1-eps(100000000));

 sectiondata.sectionType = 'Idoublysymmetric';
 sectiondata.b = 180e-3; % [m]
 sectiondata.houtside=(400)*10^(-3); % [m]
 sectiondata.tf = 13.5e-3; % [m]
 sectiondata.heff = sectiondata.houtside-sectiondata.tf; % [m]
 sectiondata.hinside=sectiondata.houtside-2*sectiondata.tf; % [m]
 sectiondata.tw = 8.6e-3; % [m]
 sectiondata.material.type = 'isotropic';
 sectiondata.material.E = 210e9; %[N/m^2]
 sectiondata.material.nu = 0.3; %[-]
 
 Iyy=(2*(sectiondata.b)*sectiondata.tf^3 + sectiondata.tw*(sectiondata.heff-sectiondata.tf)^3)/12 + 2* (sectiondata.b*sectiondata.tf) *(sectiondata.heff/2)^2;
 IYY=((sectiondata.b*sectiondata.houtside^3) - (sectiondata.b-sectiondata.tw) *sectiondata.hinside^3)/12;
 assert(abs(Iyy/IYY-1)<2e-15,'Iyy caluclations differ')
 Area=2*sectiondata.b*sectiondata.tf+sectiondata.hinside*sectiondata.tw;
 assert(Area>0,'Area negativ')
 Izz=(2*sectiondata.tf*sectiondata.b^3+sectiondata.hinside*sectiondata.tf^3)/12;
 assert(Izz>0,'Izz negativ')
 
 BpM=UBpUges/(1-UBpUges);
 yecc=sqrt(BpM*Izz/Area);
 zecc=sqrt(BpM*Iyy/Area);
 
 
end
% Izz