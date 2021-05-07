function [k,lami] = findStabilityLimit(fullEV,lambda0)

%    lami = 10e10;
%    k = 0;
%    poslam = 0;
%    for i = 2:length(lambda0)
% %        Lami = eigval{i}(5,:);
% %        [ts,pos] = min(abs(Lami));
% %        if ts<abs(lami)
% %           lami = Lami(pos);
% % %           poslam = pos;
% % %           k = i;
% %        end
%    end
   Nr2=find(fullEV(1,:)<0,1);
   if isempty(Nr2)
    [y2,Nr2]=min(fullEV(1,:));
   else
    y2=fullEV(1,Nr2);
   end
   x2=lambda0(Nr2);
   Nr1=Nr2-1;
   x1=lambda0(Nr1);
   y1=fullEV(1,Nr1);
   dk=(y1)/(y1-y2);
   k=Nr1+dk;
   lami=x1+(x2-x1)*dk;
   if isempty(lami)
    lami=lambda0(end-1)/(1-fullEV(1,end-1));
   end
end