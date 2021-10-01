function [k,lami] = InterpolateJK(xData,lambda0)

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
   if size(xData,2)==1
    fullEV=xData;
   else
    fullEV=xData(1,:);
   end
   Nr2=find(fullEV<0,1);
   if isempty(Nr2)
    [y2,Nr2]=min(fullEV);
   else
    y2=fullEV(Nr2);
   end
   x2=lambda0(Nr2);
   Nr1=max(Nr2-1,1);
   x1=lambda0(Nr1);
   y1=fullEV(Nr1);
   dk=(y1)/(y1-y2);
   k=Nr1+dk;
   lami=x1+(x2-x1)*dk;
   if isempty(lami)
    lami=lambda0(end-1)/(1-fullEV(end-1));
   end
end