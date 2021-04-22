function Displ = NodalResults2DisplJK(Nres)
NreskeysIn = Nres.keys;

if numel(NreskeysIn)>0
 for i = 1:length(NreskeysIn)
  key = NreskeysIn{i};
  if strcmp(key(1:2),'UR')
   remove(Nres, key);
  end
 end
 Nreskeys = Nres.keys;
 RES = NaN(length(Nreskeys),size(Nres(Nreskeys{1}),1),size(Nres(Nreskeys{1}),2)-1); % Richtungen(max6) x NrKnoten x NrLastschritte
 for i = 1:length(Nreskeys)
  key = Nreskeys{i};
  if strcmp(key(1:2),'UR')
   break
  end

  val = Nres(key); val = val(:,2:end);
  %val(abs(val)<=1e-31)=0;%remove numeric issues close to zero
  RES(i,:,:) = val;
 end %for
 RES(abs(RES)<=8.78e-30)=0;%remove numeric issues close to zero
 lam=size(RES,3);
 Displ=cell(lam,1);
 for i=1:lam
  RESi=RES(:,:,i);
  %RESi(abs(RESi)<=1e-31)=0;%remove numeric issues close to zero
  Displ{i} = RESi;
 end
else
 Displ=NaN(0);
end %if

end %function