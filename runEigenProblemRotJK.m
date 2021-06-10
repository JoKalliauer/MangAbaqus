function model = runEigenProblemRotJK(~,model,Displ,~,~,matches,wbrEP)
 

 displacements = cell(length(matches),1);

 

 
 %% solve EigvalueProblem

 
 
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem EigvalueProblemRot');end

 for i = 1:length(matches)
  if usejava('jvm'); waitbar(i/length(matches),wbrEP,'runEigenProblem EigvalueProblemrot');end
 %for i
  if i==0

   
  else % if i~= 1
   if matches(i)>=1
    displacements_ = Displ{matches(i)};
   end
  end
   
   displacements{i} = displacements_;



 end%for i = 1:length(matches)
 model.rotJK = displacements;

 

 
end %fucntion

