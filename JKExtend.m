function  JKExtend(model,modelprops)
 %#ok<*NASGU>
aa=model.stiffnessMatrices{1,1};%Kt0
ab=model.stiffnessMatrices{1,2};%Ktprim0;
ac=model.stiffnessMatrices{1,3};%dKtprim0
ba=model.stiffnessMatrices{2,1};
bb=model.stiffnessMatrices{2,2};
bc=model.stiffnessMatrices{2,3};

astrich=(bc-ac)/(2*modelprops.epsilon);

EVaa=model.eigenvalues{1,1};

EVstrich=(EVaa(7,1)-EVaa(6,1))/(modelprops.epsilon);

r1x=model.eigenvectors{1};
r1=transpose(r1x(8,8:end,1));
x1=ab*r1;
x2=EVstrich*aa*r1;

Fehler=2*norm(x2-x1)/norm(x2+x1) %#ok<NOPRT>






end %function