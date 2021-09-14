input_4node;

fp = fopen('4_node_result.txt','w');

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

    xy = [coord(nd1,2), coord(nd1,3);
          coord(nd2,2), coord(nd2,3);
          coord(nd3,2), coord(nd3,3);
          coord(nd4,2), coord(nd4,3)];

    vec = [nd1, nd2, nd3, nd4];                     % Global DOF location
    
    kele = zeros(4,4);
    fele = zeros(4,1);

    for gp1 = 1:ngauss                  % outer loop over gauss points
        xi = xivec(gp1);     xiw = xi_w(gp1);
        for gp2 = 1:ngauss              % inner loop over gauss points for domain terms
            eta = etavec(gp2);     etaw = eta_w(gp2);
            if el == 1
                [ked,fed] = domain(xi, eta, xy, k, Q1);
                kele(1:4,1:4) = kele(1:4,1:4) + ked*xiw*etaw;
                fele(1:4) = fele(1:4) + fed*xiw*etaw;    
            elseif el == 2
                [ked,fed] = domain(xi, eta, xy, k, Q2);
                kele(1:4,1:4) = kele(1:4,1:4) + ked*xiw*etaw;
                fele(1:4) = fele(1:4) + fed*xiw*etaw;
            end
            
        end
        
        if el==1
            [keh,feh] = gamah(xi, xy, h, Tinf, 1);
            kele(1:4,1:4) = kele(1:4,1:4) + keh*xiw;
            fele(1:4) = fele(1:4) + feh*xiw;
            
            [keh,feh] = gamah(xi, xy, h, Tinf, 3);
            kele(1:4,1:4) = kele(1:4,1:4) + keh*xiw;
            fele(1:4) = fele(1:4) + feh*xiw;
             
        elseif el == 2
            [keh,feh] = gamah(xi, xy, h, Tinf, 1);
            kele(1:4,1:4) = kele(1:4,1:4) + keh*xiw;
            fele(1:4) = fele(1:4) + feh*xiw;
            
            feh = gamaq(xi, xy, qn, 2);
            fele(1:4) = fele(1:4) + feh*xiw;
            
            feh = gamaq(xi, xy, qn, 3);
            fele(1:4) = fele(1:4) + feh*xiw;
        end
    end
    
    fprintf(fp,'\n\nElement %d:\n',el);
    fprintf(fp,'---------\n');
    fprintf(fp,'\nStiffness matrix:\n');
    for i = 1:4
        fprintf(fp,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',kele(i,1:4));
    end
    fprintf(fp,'Load vector:\n');
    for i = 1:4
        fprintf(fp,'%14.4e\n\n',fele(i));
    end
    
    % Assembly using connectivity data
    for ii = 1:4
        for jj = 1:4
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + kele(ii,jj);
        end
        F(vec(ii)) = F(vec(ii)) + fele(ii);
    end
    
end

fprintf(fp,'\n\nGlobal stiffness matrix:\n');
fprintf(fp,'=======================\n\n');
for i = 1:nnode
    fprintf(fp,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t\n\n',K(i,1:nnode));   
end

fprintf(fp,'\n\nGlobal load vector:\n');
fprintf(fp,'==================\n\n');
for i = 1:nnode
    fprintf(fp,'%14.4e\n\n',F(i));
end


% Imposition of B.C.
Kreduce = K(2:6,2:6);
Freduce = F(2:6) - K(2:6,1)*200 - K(2:6,7)*200;

% Finding Solution
Treduce = Kreduce\Freduce;
Tn = [200;Treduce;200];

fprintf(fp,'\n\nFinal nodal temperatures:\n');
fprintf(fp,'========================\n\n');
for i = 1:nnode
    fprintf(fp,'%14.4e\n\n',Tn(i));
end

fclose(fp);