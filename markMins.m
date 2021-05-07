function [] = markMins(lambdaplot,y4,bbb)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   [miny4,idx]=min(y4);
  idxs=islocalmin(y4,'FlatSelection', 'all');
  plot([0 0],[0 0],'LineWidth',eps(0));
  plot(lambdaplot(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  lambdaplot(idxs)
  y4s=y4(idxs);
  maxy4=max(y4);
  miny4=min([min(y4s(y4s>0)) maxy4]);
  if miny4<.01*maxy4 || any(miny4<0)
   untergrenze=(miny4);
   obergrenze=min(bbb.YLim(2),round(100*abs(untergrenze),2,'significant'));
   bbb.YLim = [0,obergrenze];
  end
end

