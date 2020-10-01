function model = selectModel(modelprops,AbaqusRunsFolder)
if nargin<1
    testcase = 'TL_arch';
    numofelm = 20;
    eltype = 'B21';
    lambda = [.025,.075,.125,.175,.225,.2750,.325,.375,.425,.475,.525,.575,.625,.675,.725];
    %epsil = 0.005;
    %typeofanal = 'K0';
    loadFactor = 1.0;
    len = [];
else
    testcase = modelprops.testcase;
    numofelm = modelprops.numofelm;
    eltype =   modelprops.elementtype;
    lambda =   modelprops.lambda;
    %epsil =    modelprops.epsilon;
    %typeofanal = modelprops.typeofanalysis;
    loadFactor = modelprops.loadfactor;
    len = modelprops.length;
    if isfield(modelprops,'ecc')
        ecc = modelprops.ecc;
    else
        ecc = [];
    end
end

   switch testcase
       case 'TL_arch'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch([],numofelm,lambda,loadFactor,eltype);
       case 'TL_arch3D'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch3D([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder);
       case 'TL_arch_Hinge'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch_Hinge([],numofelm,lambda,loadFactor,eltype);
       case 'TL_arch3D_Hinge'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch3D_Hinge([],numofelm,lambda,loadFactor,eltype);
       case 'pureBendingBeam'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.pureBendingBeam(len,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
       case 'cantilever'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.cantilever(len,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
       case 'eccenCompressionBeam'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.eccenCompressionBeam(len,numofelm,lambda,loadFactor,eltype,ecc);
       case 'eccenCompressionBeam2D'
           [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.eccenCompressionBeam2D(len,numofelm,lambda,loadFactor,eltype,ecc);

   end
   
   model.filename = filename;
   model.lambda = lambda;
   model.BC = BC;
   model.Nodes = Nodes;
   model.Elements = Elements;
   model.AbaqusRunsFolder = AbaqusRunsFolder;
   
end