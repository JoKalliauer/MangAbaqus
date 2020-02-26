function [ELres, Nres, EigRes] = getHistoryOutputFromDatFile(path)
if nargin<1
   path = 'AbaqusRuns/testmodel_beam.dat'; 
end
    ELres = containers.Map;
    Nres = containers.Map;
    EIGVECres = containers.Map;
    EigValList = {};

    u = fopen(path,'r');
    
    NodeMap = [];
    nm = 0;
    eigstepnum  = 0;

    while ~feof(u)
        tline = fgetl(u);
        if length(tline)>=37
            if strcmpi(tline(1:37),'GLOBAL TO LOCAL NODE AND ELEMENT MAPS')
                fgetl(u);
                fgetl(u);
                fgetl(u);
                fgetl(u);
                tline = fgetl(u);
                while length(tline)>1
                    nm = nm + 1;
                    NodeMap(nm,:) = sscanf(tline,'%d %d');
                    tline = fgetl(u);
                end
            end
        end
        if length(tline)>=20
            if strcmpi(tline(1:20),' STEP TIME COMPLETED')
                tline = '123456789012345672534868901234567890';
%                 while ~strcmpi(tline(1:25),' TIME INCREMENT COMPLETED')
                while ~strcmpi(tline(1:31),'                        S T E P')
                    tline = fgetl(u);
                    if length(tline)>=15
                        if strcmpi(tline(1:15),'    ELEMENT  PT')
                            C = textscan(tline,'%s','delimiter',{' ','FOOT-'},'MultipleDelimsAsOne',1);
                            C = C{1};

                            fgetl(u); fgetl(u);

                            tline = fgetl(u);
                            Val = []; ecount = 0;
                            while (~isempty(tline))&&(~feof(u))
                                ecount = ecount + 1;
                                Val(ecount,:) = sscanf(tline,'%f',[1 length(C)]);
                                tline = fgetl(u);
                            end
                            for i = 3:length(C)
                               if ELres.isKey(C{i})
                                   ELres(C{i}) = [ELres(C{i}),Val(:,i)];
                               else
                                   ELres(C{i}) = [Val(:,1:2),Val(:,i)];
                               end
                            end
                            if feof(u)
                                break
                            end
            %                 fseek(u,filepos,'bof');
                        end
                    end
                    if length(tline)>=17
                        if strcmpi(tline(1:17),'       NODE FOOT-')
                            C = textscan(tline,'%s','delimiter',{' ','FOOT-'},'MultipleDelimsAsOne',1);
                            C = C{1};

                            fgetl(u); fgetl(u);

                            tline = fgetl(u);
                            Val = []; ecount = 0;
                            while (~isempty(tline))&&(~feof(u))
                                ecount = ecount + 1;
                                Val(ecount,:) = sscanf(tline,'%f',[1 length(C)]);
                                tline = fgetl(u);
                            end
                            for i = 2:length(C)
                               if Nres.isKey(C{i})
                                   toupdate = [Nres(C{i}),zeros(size(Nres(C{i}),1),1)];
                                   toupdate(Val(:,1),end) = Val(:,i);
                                   Nres(C{i}) = toupdate;
                               else
                                   Nres(C{i}) = [Val(:,1),Val(:,i)];
                               end
                            end
                            if feof(u)
                                break
                            end
            %                 fseek(u,filepos,'bof');
                        end
                    end
                    if feof(u)
                        break
                    end
                    if length(tline)<25
                        tline = '12345678901234567890122456874534567890';
                    end
                end
            end
        end
        if length(tline)>=64
            if strcmpi(tline(1:64),'                              E I G E N V A L U E    O U T P U T')
                fgetl(u); fgetl(u); fgetl(u);
                fgetl(u); fgetl(u); fgetl(u);
                fgetl(u); fgetl(u); fgetl(u);
                tline = fgetl(u);
                eigstepnum = eigstepnum + 1;
                nm = 0;
                EigVal = [];
                while length(tline)>1
                    nm = nm + 1;
                    EigVal(nm,:) = sscanf(tline,'%g %g');
                    tline = fgetl(u);
                end
                EigVal(:,1) = [];
                EigValList{eigstepnum} = EigVal;
                fgetl(u); fgetl(u); fgetl(u);
                
                tline = fgetl(u);
                eignum = 0;
                while ~strcmpi(tline(1:31),'                        S T E P')
                    tline = fgetl(u);
                    if length(tline)>=17
                        if strcmpi(tline(1:17),'       NODE FOOT-')
                            eignum = eignum + 1;
                            C = textscan(tline,'%s','delimiter',{' ','FOOT-'},'MultipleDelimsAsOne',1);
                            C = C{1};

                            fgetl(u); fgetl(u);

                            tline = fgetl(u);
                            Val = []; ecount = 0;
                            while (~isempty(tline))&&(~feof(u))
                                ecount = ecount + 1;
                                Val(ecount,:) = sscanf(tline,'%f',[1 length(C)]);
                                tline = fgetl(u);
                            end
                            for i = 2:length(C)
                               if EIGVECres.isKey(C{i})
                                   toupdatelist = EIGVECres(C{i});
                                   if length(toupdatelist)<eignum
                                       toupdatelist{eignum} = [Val(:,1),Val(:,i)];
                                   else
                                       toupdate = [toupdatelist{eignum},zeros(size(toupdatelist{eignum},1),1)];
                                       toupdate(Val(:,1),end) = Val(:,i);
                                       toupdatelist{eignum} = toupdate;
                                   end
                                   EIGVECres(C{i}) = toupdatelist;
                               else
                                   toupdatelist{eignum} = [Val(:,1),Val(:,i)];
                                   EIGVECres(C{i}) = toupdatelist;
                               end
                            end
                            if feof(u)
                                break
                            end
            %                 fseek(u,filepos,'bof');
                        end
                    end
                    if feof(u)
                        break
                    end
                    if length(tline)<25
                        tline = '12345678901234567890122456874534567890';
                    end
                end
            end
        end
    end

    fclose(u);

    EigRes.eigenvalues = EigValList;
    EigRes.eigenvectors = EIGVECres;
end