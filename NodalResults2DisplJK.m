function [Displ,Rot] = NodalResults2DisplJK(Nres)
%% returns the displacement and the rotations from the nodal-results of Abaqus
%university:TU Wien
%author: Johannes Kalliauer(©2020-2023)

%% Input
% Nres ... Nodal results of Abaqus

%% Output
% Displ  ... Displacements
% Rot ... Rotations

%% recent-change
%2023-02-16 JK: adding error-identifier

%% Code

NreskeysIn = Nres.keys;
if isempty(NreskeysIn)
 %Displ=NaN;
 %Rot=NaN;
 error('MyPrgm:MSG','No Displ-Results, check *.msg-file')
 %return
end
NresDisp=containers.Map(NreskeysIn,Nres.values);
NresRot=containers.Map(NreskeysIn,Nres.values);

if numel(NreskeysIn)>0
 for i = 1:length(NreskeysIn)
  key = NreskeysIn{i};
  if strcmp(key(1),'U')
   if strcmp(key(1:2),'UR')
    remove(NresDisp, key);
   elseif strcmp(key(1:2),'U1') ||  strcmp(key(1:2),'U2') ||  strcmp(key(1:2),'U3')
    remove(NresRot, key);
   else
    error('MyProg:Strange','key not recogniced')
   end
  elseif  strcmp(key(1:2),'WA')
   remove(NresDisp, key);
   remove(NresRot, key);
  elseif  strcmp(key,'NODE')
   remove(Nres, key);
   remove(NresDisp, key);
   remove(NresRot, key);
  else
   error('MyProg:Strange','key not recogniced')
  end
 end
 NresDispkeys = NresDisp.keys;
 RES = NaN(length(NresDispkeys),size(Nres(NresDispkeys{1}),1),size(Nres(NresDispkeys{1}),2)); % Richtungen(max6) x NrKnoten x NrLastschritte
 for i = 1:length(NresDispkeys)
  key = NresDispkeys{i};


  val = Nres(key);
  %val = val(:,2:end);
  %val(abs(val)<=1e-31)=0;%remove numeric issues close to zero
  RES(i,:,:) = val;
 end %for
 RES(abs(RES)<=8.78e-30)=0;%remove numeric issues close to zero
 smallnr= (abs(RES)<=2.25e-12);
 if any(smallnr(:))
  if any(RES(smallnr))
   if any(abs(RES(smallnr))>=eps(5000))
    warning('MyPrgm:Numeric','remove numeric issues close to zero')
   end
   RES(smallnr)=0;%remove numeric issues close to zero
  end
 end
 lam=size(RES,3);
 Displ=cell(lam,1);
 for i=1:lam
  RESi=RES(:,:,i);
  %RESi(abs(RESi)<=1e-31)=0;%remove numeric issues close to zero
  Displ{i} = RESi;
 end

 NresRotkeys = NresRot.keys;
 RES = NaN(length(NresRotkeys),size(Nres(NresRotkeys{1}),1),size(Nres(NresRotkeys{1}),2)); % Richtungen(max6) x NrKnoten x NrLastschritte
 for i = 1:length(NresRotkeys)
  key = NresRotkeys{i};
  val = Nres(key);
  %val = val(:,2:end);
  %val(abs(val)<=1e-31)=0;%remove numeric issues close to zero
  RES(i,:,:) = val;
 end %for
 RES(abs(RES)<=8.78e-30)=0;%remove numeric issues close to zero
 Rot=cell(lam,1);
 for i=1:lam
  RESi=RES(:,:,i);
  %RESi(abs(RESi)<=1e-31)=0;%remove numeric issues close to zero
  Rot{i} = RESi;
 end
else
 Displ=NaN(0);
end %if

end %function