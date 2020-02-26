function integ = integrate1Dbeam(v,evalF)
if nargin < 1
   v = [-3; 4];
   evalF = [1 1];
end

v1 = v(1);
v2 = v(2);

w = [1, 1];
ptGaussRef = [-1/sqrt(3); 1/sqrt(3)];
   
% Shape functions
Psi1=@(xi) -0.5*xi + 0.5;
Psi2=@(xi)  0.5*xi + 0.5;
% Shape function derivatives
dPsi11=@(xi) -0.5;
dPsi21=@(xi)  0.5;

% Gradient matrix
Jacb =@(xi) [dPsi11(xi), dPsi21(xi)];

xi = ptGaussRef(:,1);
evalPsi1 = Psi1(xi);
evalPsi2 = Psi2(xi);

% from the change of variables function
ptGaussDomain = evalPsi1*v1+evalPsi2*v2;

% evaluate Jacobian contribution for each Gauss point
evalDetJacb = zeros(size(xi,1),1);
for i=1:size(xi,1)
    evalDetJacb(i) = abs(det(Jacb(xi(i))*v));
end

% Finally, apply Gauss integration formula
suma=0;
for i=1:size(ptGaussDomain,1)
    suma=suma + abs(w(i)*evalF(i)*evalDetJacb(i));
end
integ = suma;

end