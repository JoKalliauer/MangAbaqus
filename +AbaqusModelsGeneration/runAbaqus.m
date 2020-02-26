function runAbaqus(filename)
if nargin<1
    modelname = 'pureBendingBeamShellElm';
end

        cd AbaqusRuns
        if ~(exist([filename,'.sta'], 'file') == 2)
           system(['exec abaqus job=',filename,' cpus=2 interactive'])
%            system(['exec abq2016 job=',filename,' scratch="/scratch/tmp" cpus=1 mp_mode=threads interactive'])
        end
        cd ..

end
