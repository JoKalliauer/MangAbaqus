function [evec0,eval0,numofeigs] = solveCLEforMinEigNew(Kt,Ktprim,Eigres,Kt0_0,typeofanal,iter)
    
        numofeigs = 10;

    if nargin>0
        if numofeigs>size(Kt,1)
            numofeigs = size(Kt,1);
        end
        switch typeofanal
            case 'I'
                Ktprim = speye(size(Kt,1),size(Kt,1));
            case 'CLE'
                Ktprim = Ktprim;
            case 'I-CLE'
                Ktprim = -speye(size(Kt,1),size(Kt,1)) + Ktprim;
            case 'I+CLE'
                Ktprim = speye(size(Kt,1),size(Kt,1)) + Ktprim;
            case 'test'
                Ktprim = (Kt - Kt0_0);
            case 'test2'
                Ktprim = (Kt - Kt0_0);
                Kt = Kt0_0;     
            case 'I-CLE-test3'
                Ktprim = -speye(size(Kt,1),size(Kt,1))*sum(diag(Ktprim))/size(Ktprim,1) + Ktprim;
            case 'I-CLE-test4'
                Ktprim = speye(size(Kt,1),size(Kt,1))*sum(diag(Ktprim))/size(Ktprim,1) - Ktprim;
            case 'diag(CLE)-CLE'
                Ktprim = -diag(Ktprim)+ Ktprim;
            case 'KNL'
                Kt = Kt - Kt0_0;
                if sum(full(Kt(:)))==0
                    Kt = Kt0_0;
                end
                Ktprim = speye(size(Kt,1),size(Kt,1));
            case 'KNL2'
                if iter>1
                   Ktprim =  -Kt0_0;
                else
                   Ktprim = eye(size(Kt,1),size(Kt,1));
                end
            case 'KNL3'
                Ktprim = Kt - Kt0_0;
                Kt = Kt0_0;
            case 'KNL4'
                Kt = Kt - Kt0_0;
                Ktprim = Kt0_0;
                if iter==1
                    Kt = eye(size(Kt));
                end
            case 'Kg'
                eigenvalues = Eigres.eigenvalues{iter};
                eigenvectors = Eigres.eigenvectors;
                numofeigs = length(eigenvalues);
                
                keys = eigenvectors.keys;
                
                test = eigenvectors(keys{1});
                numofnodes = size(test{1},1);
                ndof = length(keys)*numofnodes;
                evec0 = zeros(ndof,numofeigs);
                for i = 1:length(keys)
                    k1 = keys{i};
                    dis1 = eigenvectors(k1);
                    for j = 1:numofeigs
                        evec0(i:length(keys):ndof,j) = dis1{j}(:,iter+1);
                    end
                end
                eval0 = eigenvalues;
                numbers = [1:numofeigs]';
                posit = numbers(eval0>=0);
                negat = numbers(eval0<0);
                eval0 = [eval0(posit);eval0(negat)];
                evec0 = [evec0(:,posit),evec0(:,negat)];
                for ki = 1:size(evec0,2)
                   evec0(:,ki) = evec0(:,ki)/norm(evec0(:,ki));
                end
                return
        end
        
        if full(sum(Ktprim(:)))==0
            Ktprim = speye(size(Kt,1),size(Kt,1));
        end
        if strcmpi(typeofanal,'KNL4')
           [evec,eval] = eigs(Ktprim,Kt,numofeigs,1);
%         [evec,eval] = eig(full(Kt),-full(Ktprim));
        else
           [evec,eval] = eigs(Kt,-Ktprim,numofeigs,0);
%         [evec,eval] = eig(full(Kt),-full(Ktprim));
        end
        eval = diag(eval);

%         k = 0;
%         reals = [];
%         for i=1:length(eval)
%             if isreal(eval(i)) && (abs(eval(i))>1e-10)
%                k = k+1;
%                reals(k) = i;
%             end
%         end
% 
%         evec = evec(:,reals);
%         eval = eval(reals);
        
%         evec(:,eval<0) = [];
%         eval(eval<0) = 10e16;

        [ev,s] = sort(abs(eval));
        evec0 = evec(:,s(1:numofeigs));
        eval0 = eval(s(1:numofeigs));
    else
        evec0 = []; eval0 = [];
    end
    for i = 1:size(evec0,2)
       evec0(:,i) = evec0(:,i)/norm(evec0(:,i));
    end
end