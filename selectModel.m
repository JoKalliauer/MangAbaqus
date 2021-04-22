function model = selectModel(modelprops,AbaqusRunsFolder)
if nargin<1
    testcase = 'TL_arch';
    numofelm = 20;
    eltype = 'B21';
    lambda = [.025,.075,.125,.175,.225,.2750,.325,.375,.425,.475,.525,.575,.625,.675,.725];
    %epsil = 0.005;
    %typeofanal = 'K0';
    loadFactor = 1.0;
    L = [];
else
    %testcase = modelprops.testcase;
    numofelm = modelprops.numofelm;
    eltype =   modelprops.elementtype;
    lambda =   modelprops.lambda;
    %epsil =    modelprops.epsilon;
    %typeofanal = modelprops.typeofanalysis;
    loadFactor = modelprops.loadfactor;
    if isfield(modelprops,'ecc')
        ecc = modelprops.ecc;
    else
        ecc = [];
    end
end

if ~exist('testcase','var')
 testcase = modelprops.testcase;
end
if ~exist('L','var')
 L = modelprops.length;
end
 if sum(strcmp(fieldnames(modelprops), 'MeterValue')) == 0
   modelprops.MeterValue=1;
 end

if strcmp(testcase,'detKt') && strcmp(eltype(1:2),'B2')
 testcase=[testcase '2D'];
end

   switch testcase
    case 'TL_arch'
     assert(mod(numofelm,2)==0,'number of Elements must be euqal')
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch_old'
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch_old([],numofelm,lambda,loadFactor,eltype);
     AbaqusRunsFolder='AbaqusRuns/';
    case 'TL_arch3D'
     [filename,lambda,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.TL_arch3D([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch3DKg'
     [filename,~,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.TL_arch3DKg([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch3D_sin'
     [filename,lambda,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.TL_arch3D_sin([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch_Hinge'
     assert(~strcmp(eltype(1:2),'B3'),'TL_arch_Hinge is 2D')
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch_Hinge([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch3D_Hinge'
     assert(~(strcmp(eltype,'B32OS') || strcmp(eltype,'B32OSH')),'TL_arch3D_Hinge uses rect-section not open section')
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch3D_Hinge([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'pureBendingBeamJK'
     [filename,lambda,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.pureBendingBeam(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    case 'pureBendingBeamMalendowski'
     M = loadFactor*0.5e6; %[N*m ?]
     model.load=lambda*M;
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.pureBendingBeam_Malendowski(L,numofelm,lambda,loadFactor,eltype);
    case 'cantilever'
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.cantilever(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam'
     [filename,lambda,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.eccenCompressionBeam(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam2D'
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.eccenCompressionBeam2D(L,numofelm,lambda,loadFactor,eltype,ecc);
    case 'twoBeams'
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.twoBeams(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'detKt2D'
     [filename,lambda,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.detKt2D(L,numofelm,lambda,loadFactor,eltype,[],modelprops,AbaqusRunsFolder);
    case 'd2bock'
     [filename,lambda,BC,Nodes,Elements,model.load,model.dofpNode] = AbaqusModelsGeneration.d2bock(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    otherwise 
     warning('MyProgram:unknown','testcase unknown')
   end
   
   model.filename = filename;
   model.lambda = lambda;
   model.BCMatlab = BC;
   model.Nodes = Nodes;
   model.Elements = Elements;
   model.AbaqusRunsFolder = AbaqusRunsFolder;
   
end