input_6node;

[C] = constitutive(E, nu);

fp = fopen('6_node_result.txt','w');

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
    nd5 = connect(el,6);
    nd6 = connect(el,7);

    xy = [coord(nd1,2), coord(nd1,3);
          coord(nd2,2), coord(nd2,3);
          coord(nd3,2), coord(nd3,3);
          coord(nd4,2), coord(nd4,3);
          coord(nd5,2), coord(nd5,3);
          coord(nd6,2), coord(nd6,3)];

    vec = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2, 2*nd3-1, 2*nd3, 2*nd4-1, 2*nd4, 2*nd5-1, 2*nd5, 2*nd6-1, 2*nd6];       % Global DOF location                   
    
    kele = zeros(12,12);
    fele = zeros(12,1);

    for gpd = 1:ngaussd     % loop over gauss points
        
        xi = xivec(gpd);    eta = etavec(gpd);     w = wvec(gpd);
            kele(1:12,1:12) = kele(1:12,1:12) + elestiff(xi, eta, xy, C)*w;
    end

    for gpb = 1:ngaussb     % loop over gauss points
        
        beta = betavec(gpb);       ew = ewvec(gpb);
        
        if el==2
            fele(1:12) = fele(1:12) + gamat(beta, xy, P0, P1, theta, 1)*ew;
        end

    end
    
    fprintf(fp,'\n\nElement %d:\n',el);
    fprintf(fp,'---------\n');
    fprintf(fp,'\nStiffness matrix:\n');

    for i = 1:12
        for j = 1:12
            fprintf(fp,'%14.4e\t',kele(i,j));
        end
        fprintf(fp,'\n\n');
    end

    fprintf(fp,'Load vector:\n');
    for i = 1:12
        fprintf(fp,'%14.4e\n\n',fele(i));
    end
    
    % Assembly using connectivity data
    for ii = 1:12
        for jj = 1:12
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
Kreduce = K([3:6,9:14,17:18],[3:6,9:14,17:18]);
Freduce = F([3:6,9:14,17:18]);

% Finding Solution
ureduce = Kreduce\Freduce;
un = [0;0;ureduce(1:4);0;0;ureduce(5:10);0;0;ureduce(11:12)];

% Stresses at Gauss Points
Stress_Gauss = zeros(2,4,3);     % No. of elements x No. of gauss points x 3 stress components
for ele = 1:nele
    vec = connect(ele,2:7);
    
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
    
    uvec = zeros(12,1);
    for i = 1:6
       uvec(2*i-1) = un(2*vec(i)-1);
       uvec(2*i) = un(2*vec(i));
    end
    ig = 0;
    for i = 1:4
        for j = 1:4
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
    for i = 1:4
        for j = 1:4
            ig = ig + 1;
            fprintf(fp,'xi - %10.3e, eta - %10.3e, stresses - %10.3e   %10.3e   %10.3e\n',xivec(i),etavec(j), Stress_Gauss(ele,ig,1:3));
        end
    end
    fprintf(fp,"\n");
end

fclose(fp);

