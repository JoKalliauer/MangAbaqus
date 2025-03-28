function runAbaqus(filename,AbaqusRunsFolder,modelprops)
%#!/usr/bin/env octave -q
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2023)

%% start Abaqus

%% Input
% filename .. name of the inp-file
% AbaqusRunsFolder .. location of the inp-file
% modelprops ... parameters which are used for the Abaqus-run

%% no Output (Abaqus-files are saved on disk)

%% Recent Changes
%2023-02-16 JK: cd(oldpwd)
%2023-03-31 JK: removed duplicate reply0~=127

%% Code

 if nargin<1
  %modelname = 'pureBendingBeamShellElm';
 end
 if sum(strcmp(fieldnames(modelprops), 'ask_delete')) == 0
  modelprops.ask_delete=true;
 end
 
 if ~exist(AbaqusRunsFolder, 'dir')
  if isunix
   warning('MyProgram:OS','AbaqusRunsFolder does not exist')
   mkdir(AbaqusRunsFolder);
  end
  if ispc
   warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder does not exist, therfore skipping')
   return
  end
 end
 oldpwd=pwd;
 cd(AbaqusRunsFolder) %cd /home/jkalliau/Abaqus/ownCloud/Post/MangAbaqus/AbaqusRuns% cd AbaqusRuns
 assert(isfile(strcat(filename,'.inp')),'input-file does not exist in %s',AbaqusRunsFolder)
 if ~(exist([filename,'.sta'], 'file') == 2) || modelprops.forceAbaqus==true
  if isunix
   if modelprops.forceAbaqus==true && modelprops.ask_delete==false;   system(['rm -f ',[filename,'.sta ',filename,'_STIF*.mtx']]); end
   disp(AbaqusRunsFolder)
   if modelprops.ask_delete==false || ~usejava('desktop')
    reply0=system(['exec abq cpus=1 interactive ask_delete=OFF job=',filename]) %#ok<NOPRT>
   else
    reply0=system(['exec abq cpus=1 interactive job=',filename]) %#ok<NOPRT>
   end

    %reply is 127 if abq: not found, otherwise use the reply
    if reply0~=127 
     reply=reply0;  
    else
     %/home/jkalliau/ownCloud/Linu aufgrund von Fig33x/bin/imwsrun
     if modelprops.ask_delete==false || ~usejava('desktop')
      reply=system(['exec abaqus cpus=1 interactive ask_delete=OFF job=',filename]);
     else
      reply=system(['exec abaqus cpus=1 interactive job=',filename]);
     end
    end
    if reply==1
     [~,cmdout]=system('if ifconfig cscotun0 &>/dev/null; then echo "y"; else echo "n"; fi');
     if strcmp('n',cmdout(1))
      warning('MyProgram:NoVPN','Trying to connect to VPN')
      system('bash ~/bin/private/VPN.sh');
      %reply=system(['exec abqJK.sh ',filename]);
      reply=system(['exec abq2020hf5 cpus=1 interactive job=',filename]);
      if reply==1
       warning('MyProgram:NoVPN','Abaqus/Analysis exited with error(s).')
       disp('check if directory has *.lck file')
      end
     else
      warning('MyProgram:Abaqus','Abaqus/Analysis exited with error(s)')
      disp('check if directory has write permission for user')
      disp('check if license is connected')
      disp('check if inputfile exists in the correct folder')
      disp('try to delete all filename')
      disp('check *.dat file, maybe wrong element')
     end% if strcmp('n',cmdout(1))
    elseif reply==127
     %disp(strcat('cd ',AbaqusRunsFolder))
     disp(strcat('exec imwsrun abaqus:2022 cpus=1 interactive job=',filename))
     warning('MyProgram:Abaqus','Abaqus command not found)')
    elseif reply==0
     if ~(exist([filename,'.dat'], 'file') == 2)
      warning('MyProgram:Abaqus','no *.dat-file, check VPN-connection')
     else
      %everything worked
     end
    elseif reply==2
     warning('MyProgram:Abaqus','reply2 syntax error in imwsrun')
    elseif reply==126
     warning('MyProgram:Abaqus','no permission to execute file')
    else
     warning('MyProgram:Abaqus','reply unknown)')
    end%if reply==1
   %else
%     warning('MyProgram:OSDetection','hostname not recogniced')
%     system(['exec abq2020hf5 job=',filename,' cpus=1 interactive'])
%    %end %if strcmp(name,)



  elseif ispc
   system(['call abaqus job=',filename,' cpus=1 interactive'])
  elseif ismac
   warning('MyProgram:Untested','Mac not tested, maybe try linuxcommand')
  else
   warning('MyProgram:Untested','Platform not supported')
  end %if isunix/ispc/ismac
  if (exist([filename,'.sta'], 'file') == 2)
   %cleanup useless staff
   % ls([filename,'_X*.sim'])
   system(['rm ',filename,'_X*.sim'])
   %system(['rm ',filename,'.mdl'])
   %system(['rm ',filename,'.odb'])
  end
 else
  %   if~(exist([filename,'.odb_f'], 'file') == 2)
  %    disp('no .odb_f')
  %   else
  %    disp('.odb_f available')
  %   end
  %   if~(exist([filename,'.log'], 'file') == 2)
  %    disp('no .log')
  %   else
  %    disp('.log available')
  %   end
  %   if~(exist([filename,'.lck'], 'file') == 2)
  %    disp('no .lck')
  %   else
  %    disp('.lck available')
  %   end
  %   if ~exist([model.AbaqusRunsFolder,filename,'_STIF9.mtx'],'file')
  %    %warning('MyProgramm:Missing','_STIF*.mtx missing')
  %    error('MyProgramm:Missing','_STIF*.mtx missing, try rerunning forceAbaqus=true')
  %    %return
  %   end
  %return
 end %if ~(exist([filename,'.sta'], 'file') == 2) || modelprops.forceAbaqus==true

 if isunix %linux-PC
   %cd ~/ownCloud/Post/MangAbaqus/ %cd ..
   cd(oldpwd)
 elseif ispc %windows
   cd(oldpwd)
   %cd 'C:\Users\jokal\OneDrive\Dokumente\GitHub\MangAbaqus\' %cd ..
   %cd '%USERPROFILE%\Documents\MangAbaqus\' %cd ..
   %warning('MyPrgm:OS','Please specify the location, this must be changed to the correct folder')
   %error('MyPrgm:OS','Please specify the location, this must be changed to the correct folder')
 elseif ismac %Mac-PC
     warning('MyProgram:OS','not tested on Mac-pc')
 else
     warning('MyProgram:OS','OS unknown')
 end
 
 if ~exist([AbaqusRunsFolder,filename,'_STIF9.mtx'],'file')
  %warning('MyProgramm:Missing','_STIF*.mtx missing')
  if ~exist([AbaqusRunsFolder,filename,'_STIF7.mtx'],'file')
   AbaqusRunsFolder %#ok<NOPRT>
   %warning('MyProgramm:Missing','_STIF*.mtx missing in %s , try rerunning forceAbaqus=true or switch on VPN',AbaqusRunsFolder)
   error('MyProgramm:Missing','_STIF*.mtx missing in %s , try rerunning forceAbaqus=true or switch on VPN or comand "abaqus" not found',AbaqusRunsFolder)
   %return
  else
   warning('MyProgram:Abaqus','only few stif existing, maybe abaqus failed?')
  end
 end

end
