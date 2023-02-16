function [sectiondata,modelpropsprofil] =Profil(name,modelpropsprofil)
%#!/usr/bin/env octave -q
%university:TU Wien
%author:Johannes Kalliauer
%created: 2022

%% Defines predefined Profil-sets

%% Input
% name ... name of the predifend Set
% modelpropsprofil ... profil-data of modelprops before adding more data (optional)

%% Output
% sectiondata ... data of the cross-section
% modelpropsprofil ... profil-data of the predifend profil-set


%% Recent Changes
%2023-02-16 JK: added comments for explanation; added unit-comments

 MV=1; %[m] in and output is in meters, Abaqus might use differnet units

if ~exist('name','var')
 name='default';
end

if strcmp(name,'default')
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
 
elseif strcmp(name,'PavlicekPage93')
 Emodul=2e+11/MV; % [N/m^2]
 EA=5*10^7*MV; %[N]
 EI=10^7*MV^3; %[N * m^2]
 itragheit=sqrt(EI/EA); % [m]
 modelpropsprofil.h=itragheit*sqrt(12);
 A=EA/Emodul;
 modelpropsprofil.b=A/modelpropsprofil.h;
 sectiondata=modelpropsprofil;
elseif strcmp(name,'MalendowskiTLArch')
 modelpropsprofil.h=0.2;%[m]
 modelpropsprofil.b=0.1;%[m]
 sectiondata=modelpropsprofil;
else
 error('MyPrgm:notImplemented','not implemented')
end



end