function [membrane, nonmembrane] = selectMembraneNonMembrane(val,key)

if strcmpi(key,'epsx')||strcmpi(key,'epsy')||strcmpi(key,'epsz')
    membrane = sum(val,1);
    nonmembrane = 0*membrane;
else
    nonmembrane = sum(val,1);
    membrane = 0*nonmembrane;
end

end