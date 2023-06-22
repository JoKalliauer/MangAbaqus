function model = selectModel(modelprops,AbaqusRunsFolder)
%% calls a function based on modelprops.typeofanalysis which creates the inp-files for Abaqus
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2023)

%% Input
% modelprops ... input-parameters from the script

%% Output
% model ... everything that should be returned to Abaqus_single_run

%% Recent Changes
%2023-02-16 JK: isfield(modelprops,'length')

%% Code
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
 if isfield(modelprops,'length')
  L = modelprops.length;
 end
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
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.TL_arch3D([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
     model.xlabelloadname='line load $p$ [N/m]';
    case 'Kreis_arch3D'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.Kreis_arch3D([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch3DKg'
     [filename,~,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.TL_arch3DKg([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch3D_sin'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.TL_arch3D_sin([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch_Hinge'
     assert(~strcmp(eltype(1:2),'B3'),'TL_arch_Hinge is 2D')
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch_Hinge([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch3D_Hinge'
     assert(~(strcmp(eltype,'B32OS') || strcmp(eltype,'B32OSH')),'TL_arch3D_Hinge uses rect-section not open section')
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch3D_Hinge([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'pureBendingBeamJK'
     if modelprops.numofelm==1
      if strcmp(modelprops.elementtype,'B31OSH')
       warning('MyPrgm:Input','B31OSH need at least 2 elements to calculate a stiffnessmatrix')
      end
     end
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.sectiondata.houtside] = AbaqusModelsGeneration.pureBendingBeam(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
     model.xlabelloadname='bending moment $M$ [N\,m]';
%      model.xValload=model.load;
    case 'pureBendingBeamJK10'
     if modelprops.numofelm==1
      if strcmp(modelprops.elementtype,'B31OSH')
       warning('MyPrgm:Input','B31OSH need at least 2 elements to calculate a stiffnessmatrix')
      end
     end
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.sectiondata.houtside] = AbaqusModelsGeneration.pureBendingBeam10(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
     model.xlabelloadname='bending moment $M$ [N\,m]';
    case 'pureBendingBeamMalendowski'
     M = loadFactor*0.5e6; %[N*m ?]
     model.load=lambda*M;
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.pureBendingBeam_Malendowski(L,numofelm,lambda,loadFactor,eltype);
    case 'pureBendingCantilever'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.pureBendingCantilever(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
     model.xlabelloadname='bending moment $M$ [N\,m]';
    case 'mixedCantilever'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.mixedCantilever(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    case 'RotatedCantilever'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.RotatedCantilever(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    case 'cantilever'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.cantilever(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam'
     numofelm=numofelm(1);
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.sectiondata.houtside] = AbaqusModelsGeneration.eccenCompressionBeam(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam64'
     numofelm=numofelm(1);
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.eccenCompressionBeam64(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam2D'
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.eccenCompressionBeam2D(L,numofelm,lambda,loadFactor,eltype,ecc);
    case 'twoBeams'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.twoBeams(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'detKt2D'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.detKt2D(L,numofelm,lambda,loadFactor,eltype,[],modelprops,AbaqusRunsFolder);
     model.xlabelloadname='axial load $N$ [N]';
    case 'detKt2Dneu'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.detKt2Dneu(L,numofelm,lambda,loadFactor,eltype,[],modelprops,AbaqusRunsFolder);
     model.xlabelloadname='axial load $N$ [N]';
    case 'detKt2Dgen'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.detKt2D(L,numofelm,lambda,loadFactor,eltype,[],modelprops,AbaqusRunsFolder);
     model.xlabelloadname='axial load $N$ [N]';
    case 'd2bock'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.JC] = AbaqusModelsGeneration.d2bock(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    case 'd2bockDisp'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.JC] = AbaqusModelsGeneration.d2bockDisp(L,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder);
    otherwise 
     warning('MyProgram:unknown','testcase "%s" unknown',testcase)
   end
   
   model.filename = filename;% model.filename = replace(filename, ' ', '_');
   model.fulllambda = lambda;
   model.lambda = lambda(1:end-3);
   model.BCMatlab = BC;
   model.Nodes = Nodes;
   model.Elements = Elements;
   model.AbaqusRunsFolder = AbaqusRunsFolder;
   %model.fullload=model.load;
   model.load=model.fullload(1:end-3);
   
end