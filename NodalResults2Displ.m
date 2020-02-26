function Displ = NodalResults2Displ(Nres)
Nreskeys = Nres.keys;

RES = zeros(length(Nreskeys)*size(Nres(Nreskeys{1}),1),size(Nres(Nreskeys{1}),2)-1);
for i = 1:length(Nreskeys)
    key = Nreskeys{i};
    val = Nres(key); val = val(:,2:end);
    RES(i:length(Nreskeys):end,:) = val;    
end
Displ = mat2cell(RES,size(RES,1),ones(1,size(RES,2)))';
end