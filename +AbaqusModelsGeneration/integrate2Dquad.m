function integ = integrate2Dquad(v,evalF)
if nargin<1
    % The domain vertices
    v1=[0,0];
    v2=[5,-1];
    v3=[4,5];
    v4=[1,4];
    v = [v1; v2; v3; v4];
    evalF = [1 1 1 1];
end


%% 2D Quad:
v1 = v(1,:);
v2 = v(2,:);
v3 = v(3,:);
v4 = v(4,:);

w = [1 1 1 1];
ptGaussRef = sqrt(3)/3*[-1 -1; 1 -1; -1 1; 1 1];

% Shape functions
Psi1=@(xi,eta)(1-xi).*(1-eta)/4;
Psi2=@(xi,eta)(1+xi).*(1-eta)/4;
Psi3=@(xi,eta)(1+xi).*(1+eta)/4;
Psi4=@(xi,eta)(1-xi).*(1+eta)/4;
% Shape function derivatives
dPsi11=@(xi,eta) -(1-eta)/4;
dPsi21=@(xi,eta) (1-eta)/4;
dPsi31=@(xi,eta) (1+eta)/4;
dPsi41=@(xi,eta) -(1+eta)/4;
dPsi12=@(xi,eta) -(1-xi)/4;
dPsi22=@(xi,eta) -(1+xi)/4;
dPsi32=@(xi,eta) (1+xi)/4;
dPsi42=@(xi,eta) (1-xi)/4;
% Gradient matrix
Jacb =@(xi,eta) [dPsi11(xi,eta), dPsi21(xi,eta),dPsi31(xi,eta),dPsi41(xi,eta);...
                 dPsi12(xi,eta), dPsi22(xi,eta),dPsi32(xi,eta),dPsi42(xi,eta)];
             
xi = ptGaussRef(:,1);
eta = ptGaussRef(:,2);
evalPsi1 = Psi1(xi,eta);
evalPsi2 = Psi2(xi,eta);
evalPsi3 = Psi3(xi,eta);
evalPsi4 = Psi4(xi,eta);
% from the change of variables function
ptGaussDomain = evalPsi1*v1+evalPsi2*v2+evalPsi3*v3+evalPsi4*v4;

% evaluate Jacobian contribution for each Gauss point
evalDetJacb = zeros(size(xi,1),1);
for i=1:size(xi,1)
    evalDetJacb(i) = abs(det(Jacb(xi(i),eta(i))*v));
end

% Finally, apply Gauss integration formula
suma=0;
for i=1:size(ptGaussDomain,1)
    suma=suma + abs(w(i)*evalF(i)*evalDetJacb(i));
end
integ = suma;

end