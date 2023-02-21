function [val,is,indi] = compareEigenmodes(R,r0)

      cres = abs(r0'*R);
      [val,is] = sort(cres,'d');
      if r0'*R(:,is(1))<0
          indi = -1;
      else
          indi = 1;
      end

end