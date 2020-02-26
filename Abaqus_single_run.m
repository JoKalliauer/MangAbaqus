function [res,model] = Abaqus_single_run(modelprops,sortType,noplot,forcedeig)
if nargin<1
    % there are following predefined test cases:
    % testcase = 'TL_arch';
    % testcase = 'TL_arch3D';
    % testcase = 'TL_arch_Hinge';
    % testcase = 'TL_arch3D_Hinge';
    % testcase = 'pureBendingBeam';
%     testcase = 'cantilever';
    testcase = 'eccenCompressionBeam';
%     testcase = 'eccenCompressionBeam2D';
    
%     modelprops.length = [];
    modelprops.length = 5;
%     modelprops.ecc = [];
    modelprops.ecc = 0.1;
    
    % possible element types (be aware of 2D and 3D):
%     eltype = 'B22';
    % eltype = 'B22H';
    % eltype = 'B23';
    % eltype = 'B31';
%     eltype = 'B32';
    % eltype = 'B33';
    % eltype = 'B31OS';
%     eltype = 'B32OS';
    % eltype = 'B31OSH';
    eltype = 'B32OSH';
    
    % possible types of analysis
    % typeofanal = 'I'
    % typeofanal = 'CLE'
    typeofanal = 'KNL2';
    % typeofanal = 'KNL3'
    % typeofanal = 'KNL4'
    
    numberofelements = 10;
    
    epsil = 0.005;  % finite difference step
    sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
    noplot = 0; % 1 for no plot, 0 for plot
    forcedeig = [1]; % forced eigenvector number 'none' sorting
    
       modelprops.testcase = testcase;
       modelprops.numofelm = numberofelements;
       modelprops.elementtype = eltype;
       
       modelprops.lambda = [5*epsil:10*epsil:(0.78-4*epsil)]; % do not go over snap-through point
       
       modelprops.epsilon = epsil;
       modelprops.typeofanalysis = typeofanal;
       modelprops.loadfactor = 1.0;
       %
       
end

    model = runEigenProblem(modelprops);
    if (strcmpi(sortType,'none'))&&isempty(forcedeig)
       for k3 = 1:1:10
           res = sortEigenValuesAndGetQuantities(model,sortType,noplot,k3);
       end
    else
       res = sortEigenValuesAndGetQuantities(model,sortType,noplot,forcedeig);
    end
end