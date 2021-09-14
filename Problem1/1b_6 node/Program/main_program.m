input_6node;

fp = fopen('6_node_result.txt','w');

% Global stiffness matrix and global load vector
% ----------------------------------------------
K = zeros(nnode,nnode);
F = zeros(nnode,1);

fprintf(fp,'Elemental stiffness matrices and Elemental load vectors:\n');
fprintf(fp,'=======================================================');
% Calculation of element stiffness matrix and element load vector
% ---------------------------------------------------------------
for el = 1:nele     % loop over elements
    
    nd1 = connect(el,2);
    nd2 = connect(el,3);
    nd3 = connect(el,4);
    nd4 = connect(el,5);
    nd5 = connect(el,6);
    nd6 = connect(el,7);
    
    xy = [coord(nd1,2), coord(nd1,3);
          coord(nd2,2), coord(nd2,3);
          coord(nd3,2), coord(nd3,3);
          coord(nd4,2), coord(nd4,3);
          coord(nd5,2), coord(nd5,3);
          coord(nd6,2), coord(nd6,3)];

    vec = [nd1, nd2, nd3, nd4, nd5, nd6];                     % Global DOF location
    
    kele = zeros(6,6);
    fele = zeros(6,1);

    for gpd = 1:ngaussd     % loop over gauss points for domain terms
        
        xi = xivec(gpd);    eta = etavec(gpd);     w = wvec(gpd);
        if el==1 || el==2
            [ked,fed] = domain(xi, eta, xy, k, Q1);
            kele(1:6,1:6) = kele(1:6,1:6) + ked*w;
            fele(1:6) = fele(1:6) + fed*w;
        elseif el == 3 || el==4
            [ked,fed] = domain(xi, eta, xy, k, Q2);
            kele(1:6,1:6) = kele(1:6,1:6) + ked*w;
            fele(1:6) = fele(1:6) + fed*w;
        end
        
    end

    for gpb = 1:ngaussb     % loop over gauss points for boundary terms
        
        beta = betavec(gpb);       ew = ewvec(gpb);
        
        if el==1
            [keh,feh] = gamah(beta, xy, h, Tinf, 1);
            kele(1:6,1:6) = kele(1:6,1:6) + keh*ew;
            fele(1:6) = fele(1:6) + feh*ew;
        elseif el == 2
            [keh,feh] = gamah(beta, xy, h, Tinf, 2);
            kele(1:6,1:6) = kele(1:6,1:6) + keh*ew;
            fele(1:6) = fele(1:6) + feh*ew;
        elseif el == 3
            feh = gamaq(beta, xy, qn, 2);
            fele(1:6) = fele(1:6) + feh*ew;
        elseif el == 4
            [keh,feh] = gamah(beta, xy, h, Tinf, 1);
            kele(1:6,1:6) = kele(1:6,1:6) + keh*ew;
            fele(1:6) = fele(1:6) + feh*ew;
            
            feh = gamaq(beta, xy, qn, 2);
            fele(1:6) = fele(1:6) + feh*ew;
        end

    end
    
    fprintf(fp,'\n\nElement %d:\n',el);
    fprintf(fp,'---------\n');
    fprintf(fp,'\nStiffness matrix:\n');
    for i = 1:6
        fprintf(fp,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',kele(i,1:6));
    end
    fprintf(fp,'Load vector:\n');
    for i = 1:6
        fprintf(fp,'%14.4e\n\n',fele(i));
    end
    
    % Assembly using connectivity data
    for ii = 1:6
        for jj = 1:6
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + kele(ii,jj);
        end
        F(vec(ii)) = F(vec(ii)) + fele(ii);
    end
    
end

fprintf(fp,'\n\nGlobal stiffness matrix:\n');
fprintf(fp,'=======================\n\n');
for i = 1:nnode
    fprintf(fp,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',K(i,1:nnode));   
end

fprintf(fp,'\n\nGlobal load vector:\n');
fprintf(fp,'==================\n\n');
for i = 1:nnode
    fprintf(fp,'%14.4e\n\n',F(i));
end


% Imposition of B.C.
Kreduce = K([2:6,8:12,14:16],[2:6,8:12,14:16]);
Freduce = F([2:6,8:12,14:16]) - (K([2:6,8:12,14:16],1)*200 + K([2:6,8:12,14:16],7)*200 + K([2:6,8:12,14:16],13)*200);

% Finding Solution
Treduce = Kreduce\Freduce;
Tn = [200;Treduce(1:5);200;Treduce(6:10);200;Treduce(11:13)];

fprintf(fp,'\n\nFinal nodal temperatures:\n');
fprintf(fp,'========================\n\n');
for i = 1:nnode
    fprintf(fp,'%14.4e\n\n',Tn(i));
end

fclose(fp);