function [membrane, nonmembrane] = GetEnergies(ELres,Nodes,Elements)

if ELres.isKey('STH')
    ShellThick = ELres('STH');
    ShellThick = ShellThick(1:4:end,3);
    RESshell = AbaqusModelsGeneration.EnergiesShell(ELres,Nodes,Elements);
    
    RESshellKeys = RESshell.keys;

    membrane = 0;
    nonmembrane = membrane;
    for i = 1:length(RESshellKeys)
       key =  RESshellKeys{i};
       val1 = RESshell(key);

       [membrane1, nonmembrane1] = AbaqusModelsGeneration.selectMembraneNonMembrane(val1,key);

       membrane = membrane + membrane1;
       nonmembrane = nonmembrane + nonmembrane1;
    end    
else
    RESbeam = AbaqusModelsGeneration.EnergiesBeam(ELres,Nodes,Elements);
    
    if size(Nodes,2)<4
        Nodes(:,4) = zeros(size(Nodes,1),1);
    end
    beamNodes(:,2) = Nodes(:,2) - min(Nodes(:,2));
    beamNodes(:,3) = Nodes(:,3) - min(Nodes(:,3));
    beamNodes(:,4) = Nodes(:,4) - min(Nodes(:,4));
    beamcoords = sqrt(beamNodes(Elements(:,3),2).^2 + beamNodes(Elements(:,3),3).^2 + beamNodes(Elements(:,3),4).^2);
    keys_beam = RESbeam.keys;

    membrane = 0; %zeros(1,size(RESbeam(keys_beam{1}),2));
    nonmembrane = membrane;
    for i = 1:length(keys_beam)
       key =  keys_beam{i};
       val = RESbeam(key);
       [membraneb, nonmembraneb] = AbaqusModelsGeneration.selectMembraneNonMembrane(val,key);
       membrane = membrane + membraneb;
       nonmembrane = nonmembrane + nonmembraneb;
    end
end

end