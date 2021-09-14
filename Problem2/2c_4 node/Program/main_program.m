input_4node;

[C] = constitutive(E, nu);

fp = fopen('4_node_result.txt','w');

% Global stiffness matrix and global load vector
% ----------------------------------------------
K = zeros(2*nnode,2*nnode);
F = zeros(2*nnode,1);

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

    vec = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2, 2*nd3-1, 2*nd3, 2*nd4-1, 2*nd4];       % Global DOF location                   
    
    kele = zeros(8,8);
    fele = zeros(8,1);

    for gp1 = 1:ngauss                  % loop over gauss points
        xi = xivec(gp1);    xiw = xi_w(gp1);    
        for gp2 = 1:ngauss               % inner loop over gauss points
            eta = etavec(gp2);     etaw = eta_w(gp2);
            kele(1:8,1:8) = kele(1:8,1:8) + elestiff(xi, eta, xy, C)*xiw*etaw;
        end

        if el==2
            fele(1:8) = fele(1:8) + gamat(xi, xy, P0, P1, theta, 2)*xiw;
        end
        
    end

    
    fprintf(fp,'\n\nElement %d:\n',el);
    fprintf(fp,'---------\n');
    fprintf(fp,'\nStiffness matrix:\n');

    for i = 1:8
        for j = 1:8
            fprintf(fp,'%14.4e\t',kele(i,j));
        end
        fprintf(fp,'\n\n');
    end

    fprintf(fp,'Load vector:\n');
    for i = 1:8
        fprintf(fp,'%14.4e\n\n',fele(i));
    end
    
    % Assembly using connectivity data
    for ii = 1:8
        for jj = 1:8
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + kele(ii,jj);
        end
        F(vec(ii)) = F(vec(ii)) + fele(ii);
    end
    
end

fprintf(fp,'\n\nGlobal stiffness matrix:\n');
fprintf(fp,'=======================\n\n');
for i = 1:2*nnode
    for j = 1:2*nnode
        fprintf(fp,'%14.4e\t',K(i,j));
    end
    fprintf(fp,'\n\n');
end

fprintf(fp,'\n\nGlobal load vector:\n');
fprintf(fp,'==================\n\n');
for i = 1:2*nnode
    fprintf(fp,'%14.4e\n\n',F(i));
end


% Imposition of B.C.
Kreduce = K(3:10,3:10);
Freduce = F(3:10);

% Finding Solution
ureduce = Kreduce\Freduce;
un = [0;0;ureduce;0;0];

% Stresses at Gauss Points
Stress_Gauss = zeros(2,4,3);     % No. of elements x No. of gauss points x 3 stress components
for ele = 1:nele
    vec = connect(ele,2:5);
    
    nd1 = connect(el,2);
    nd2 = connect(el,3);
    nd3 = connect(el,4);
    nd4 = connect(el,5);

    xy = [coord(nd1,2), coord(nd1,3);
          coord(nd2,2), coord(nd2,3);
          coord(nd3,2), coord(nd3,3);
          coord(nd4,2), coord(nd4,3)];
    
    uvec = zeros(8,1);
    for i = 1:4
       uvec(2*i-1) = un(2*vec(i)-1);
       uvec(2*i) = un(2*vec(i));
    end
    ig = 0;
    for i = 1:2
        for j = 1:2
            ig = ig + 1;
            [strs] = stress(xivec(i), etavec(j), xy, C, uvec);
            Stress_Gauss(ele,ig,1:3) = strs;
        end
    end
end



fprintf(fp,'\n\nFinal displacements:\n');
fprintf(fp,'========================\n\n');
for i = 1:2*nnode
    fprintf(fp,'%14.4e\n\n',un(i));
end

fprintf(fp,'\n\nStresses at Gauss Points\n');
fprintf(fp,'=========================\n\n');
for ele = 1:nele
    fprintf(fp,'Element No.: %d\n',ele);
    ig = 0;
    for i = 1:2
        for j = 1:2
            ig = ig + 1;
            fprintf(fp,'xi - %10.3e, eta - %10.3e, stresses - %10.3e   %10.3e   %10.3e\n',xivec(i),etavec(j), Stress_Gauss(ele,ig,1:3));
        end
    end
end

fclose(fp);
