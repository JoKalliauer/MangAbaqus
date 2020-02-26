function [k,lami,poslam] = findStabilityLimit(eigval,lambda0)
if nargin<1
    load('EigProblem_003.mat','eigval','lambda0');
end
   lami = 10e10;
   k = 0;
   poslam = 0;
   for i = 2:length(lambda0)
       Lami = eigval{i}(5,:);
       [ts,pos] = min(abs(Lami));
       if ts<abs(lami)
          lami = Lami(pos);
          poslam = pos;
          k = i;
       end
   end

end