function runAbaqus(filename,AbaqusRunsFolder,modelprops)
 if nargin<1
  %modelname = 'pureBendingBeamShellElm';
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
 cd(AbaqusRunsFolder) %cd /home/jkalliau/Abaqus/ownCloud/Post/MangAbaqus/AbaqusRuns% cd AbaqusRuns
 if ~(exist([filename,'.sta'], 'file') == 2) || modelprops.forceAbaqus==true
  system(['rm -f ',[filename,'_STIF*.mtx']]);
  if isunix
   [~, name] = system('hostname');
   if strcmp(name(1:17),'jkalliau-Z87M-D3H') || strcmp(name(1:21),'localhost.localdomain') ||strcmp(name(1:22),'dhcp46.imws.tuwien.ac.') ||strcmp(name(1:29),'e250-200.eduroam.tuwien.ac.at')
    reply=system(['exec abq2020hf5 cpus=1 interactive job=',filename]);
    if reply==1
     [~,cmdout]=system('if ifconfig cscotun0 &>/dev/null; then echo "y"; else echo "n"; fi');
     if strcmp('n',cmdout(1))
      warning('MyProgram:NoVPN','Trying to connect to VPN')
      system('bash ~/bin/private/VPN.sh');
      %reply=system(['exec abqJK.sh ',filename]);
      reply=system(['exec abq2020hf5 cpus=1 interactive job=',filename]);
      if reply==1
       warning('MyProgram:NoVPN','Abaqus/Analysis exited with error(s).')
      end
     else
      warning('MyProgram:Abaqus','Abaqus/Analysis exited with error(s)')
      disp('check if directory has write permission for user')
      disp('check if license is connected')
      disp('check if inputfile exists in the correct folder')
     end
    elseif reply==127
     warning('MyProgram:Abaqus','Abaqus command not found)')
    elseif reply==0
     %everything worked
    else
     warning('MyProgram:Abaqus','reply unknown)')
    end
   else
    warning('MyProgram:OSDetection','hostname not recogniced')
    system(['exec abq2020hf5 job=',filename,' cpus=1 interactive'])
    %system(['exec abaqus job=',filename,' cpus=2 interactive'])
    %system(['exec abq2016 job=',filename,' scratch="/scratch/tmp" cpus=1 mp_mode=threads interactive'])
   end
  elseif ispc
   system(['call abaqus job=',filename,' cpus=1 interactive'])
  elseif ismac
   warning('MyProgram:Untested','Mac not tested, maybe try linuxcommand')
  else
   warning('MyProgram:Untested','Platform not supported')
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
 end
 cd ~/ownCloud/Post/MangAbaqus/ %cd ..

end
