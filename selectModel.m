function model = selectModel(modelprops,AbaqusRunsFolder)
%% calls a function based on modelprops.typeofanalysis which creates the inp-files for Abaqus
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2024)

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
    %lambda = [.025,.075,.125,.175,.225,.2750,.325,.375,.425,.475,.525,.575,.625,.675,.725];
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
    case 'TL_arch2D'
     assert(mod(numofelm,2)==0,'number of Elements must be even')
     RefLastMalendowski=83.3*10^3*10^2; % [N/m]
     if strcmp(modelprops.RefLast,'Eh')
      h = 20e-2;
      E=2e+11;
      RefLast=E*h;%
     else
      RefLast=loadFactor*RefLastMalendowski;%Malendowski's Load
     end
     loadFactor=RefLast/RefLastMalendowski;
     lambda=lambda/loadFactor;
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch([],numofelm,lambda,RefLast,eltype,AbaqusRunsFolder,modelprops);
    case 'TL_arch_old'
     [filename,lambda,BC,Nodes,Elements] = AbaqusModelsGeneration.TL_arch_old([],numofelm,lambda,loadFactor,eltype);
     AbaqusRunsFolder='AbaqusRuns/';
    case 'TL_arch3D'
     % RefLastMalendowski=83.3*10^3*10^2;
     % if strcmp(modelprops.RefLast,'Eh')
     %  h = modelprops.profil.h; %20e-2;
     %  E=modelprops.profil.E; % 2e+11;
     %  RefLast=E*h;%
     % else
     %  RefLast=loadFactor*RefLastMalendowski;%Malendowski's Load
     % end
     % loadFactor=RefLast/RefLastMalendowski;
     % lambda=lambda/loadFactor;
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.TL_arch3D([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
     model.xlabelloadname='line load $p$ [N/m]';
    case 'Kreis_arch3D'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.Kreis_arch3D([],numofelm,lambda,loadFactor,eltype,AbaqusRunsFolder,modelprops);
    case 'Kreis_2024'
     if strcmp(modelprops.RefLast,'EA')
      RefLast=modelprops.profil.h*modelprops.profil.b*modelprops.profil.E;%4000 MN
     elseif strcmp(modelprops.RefLast,'Wende')
      RefLast=150000;%untested 0.15MN
     else
      RefLast=1000000;%first try 1MN
     end
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.Kreis_2024(numofelm,lambda,RefLast,eltype,AbaqusRunsFolder,modelprops);
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
     % MV=modelprops.MeterValue;
     % RefLastMalendowski=0.5e6*MV*MV; % [Nm]
     % modelprops.profil.E=2.1e+11/MV;
     % if sum(strcmp(fieldnames(modelprops), 'orientate')) == 0
     %  modelprops.orientate=5;
     %  %error('not tested')
     % end
     % if strcmp(modelprops.RefLast,'ES')
     %  if modelprops.orientate==5
     %   Sy=653.6*10^-6*MV*MV*MV;
     %   %Sz=55.07*10^-6*MV*MV*MV;
     %   RefLast=modelprops.profil.E*Sy;
     %  end
     % else
     %   RefLast=loadFactor*0.5e6*MV*MV; %[Nm]
     % end
     %loadFactor=RefLast/RefLastMalendowski;
     %lambda=lambda/loadFactor;
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
    case 'ec2CompressionBeam'
     numofelm=numofelm(1);
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.sectiondata.houtside] = AbaqusModelsGeneration.ec2CompressionBeam(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam64'
     numofelm=numofelm(1);
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode] = AbaqusModelsGeneration.eccenCompressionBeam64(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
    case 'eccenCompressionBeam2D'
     [filename,lambda,BC,Nodes,Elements,model.fullload,model.dofpNode,model.sectiondata.houtside] = AbaqusModelsGeneration.eccenCompressionBeam2D(L,numofelm,lambda,loadFactor,eltype,ecc,modelprops,AbaqusRunsFolder);
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
   %modelprops.lambda=lambda;
   model.BCMatlab = BC;
   model.Nodes = Nodes;
   model.Elements = Elements;
   model.AbaqusRunsFolder = AbaqusRunsFolder;
   %model.fullload=model.load;
   model.load=model.fullload(1:end-3);
   
end