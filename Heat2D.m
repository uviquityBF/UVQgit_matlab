%%%%%%%%% FEM program for 2D heat problem using T3 element %%%%%%%%%%%%%%%%

% Developed by Gifari Zulkarnaen, July 2020
% gifari.asari@gmail.com

% References:
% Logan, D. L. (2011). A first course in the finite element method.
% Khennane, A. (2013). Introduction to finite element analysis using MATLAB® and abaqus.

clear
clear global
close all
global nnd nel nne nodof eldof n Block
global geom connec convect nf Nodal_loads
global Length Width NXE NYE dhx dhy X_origin Y_origin
format long g

%% Inputs
% DOF
    nne = 3;            % Number of nodes each element
    nodof = 1;          % Number of degrees of freedom
    eldof = nne*nodof;  % Element degrees of freedom
% Material
    Kxx = 1;           % Coefficient of thermal conductivity, Btu/(h-ft-F)
    Kyy = Kxx;
    h = 1;             % Heat convection coefficient, Btu/(h-ft-F)
    thick = 1;          % Beam thickness in mm
    
% Geometries
%     Block = 4;          % Number of rectangular block in the model
%     a = 70/7;           % Constant geometry parameter
%     Length = [a 60 60 60]; % Length of blocks
%     Width = [70 a a a]; % Width of blocks
%     X_origin = [0 10 10 10]; % X origin of the global coordinate system
%     Y_origin = [0 10 30 50]; % Y origin of the global coordinate system
%     
%     NAE = 2;            % Number of elements in each "a" parameter
%     NXE = [1 6 6 6]*NAE;% Number of elements in the x direction of each block
%     NYE = [7 1 1 1]*NAE;% Number of elements in the y direction of each block
%     dhx = Length./NXE;	% Element size in the x direction
%     dhy = Width./NYE;	% Element size in the x direction
%  
    Block = 2;          % Number of rectangular block in the model
    a = 70/7;           % Constant geometry parameter
    Length = [a 60 ]; % Length of blocks
    Width = [70 a ]; % Width of blocks
    X_origin = [0 10 ]; % X origin of the global coordinate system
    Y_origin = [0 10 ]; % Y origin of the global coordinate system
    
    NAE = 2;            % Number of elements in each "a" parameter
    NXE = [1 6 ]*NAE;% Number of elements in the x direction of each block
    NYE = [7 1 ]*NAE;% Number of elements in the y direction of each block
    dhx = Length./NXE;	% Element size in the x direction
    dhy = Width./NYE;	% Element size in the x direction
    



%     Block = 4;          % Number of rectangular block in the model
%     a = 70/7;           % Constant geometry parameter
%     Length = [70 a a a]; %  Length of blocks (X)
%     Width = [a 60 60 60]; % Width of blocks  (Y)
%     X_origin = [0 10 30 50]; % X origin of the global coordinate system
%     Y_origin = [0 10 10 10]; % Y origin of the global coordinate system
%     
%     NAE = 2;            % Number of elements in each "a" parameter
%     NXE = [7 1 1 1]*NAE;% Number of elements in the y direction of each block
%     NYE = [1 6 6 6]*NAE;% Number of elements in the x direction of each block
%     dhx = Length./NXE;	% Element size in the x direction
%     dhy = Width./NYE;	% Element size in the x direction
    
     
    
    
% Generate mesh (choose one)
     T3_mesh_final_term_axisymmetric
%    T3_mesh_final_term_unaxisymmetric

% Boundaries & Loading
    % Nodes subjected to given temperature source
    T_fix = 100; % F
    nf = ones(nnd,nodof);
    Nodal_loads= zeros(nnd,nodof);
    for i=1:nnd
        % Find nodes at x = 0
        if geom(i,1) == 0
            Nodal_loads(i,1) = T_fix;
        end
    end
    % Element edges subjected to convection of free-steam temperature
    T_inf = 20; % F
    convect = zeros(nel,nne);
    k = [2 3 1];
    for i=1:nel
        for j=1:nne
            % Find x = 70
            if all(geom(connec(i,[j k(j)]),1) == 70)
                convect(i,j) = 1;
            end
            
            % Find y = long block top/bottom & x >= 10
            find_y = [0 70];
            for p = 1:length(find_y)
                isConvect1(p) = all(geom(connec(i,[j k(j)]),2) == find_y(p));
            end
            if any(isConvect1) && all(geom(connec(i,[j k(j)]),1) >= 0)
                  convect(i,j) = 1;
            end
            
            % Find y = long block top/bottom & x >= 10
            find_y = [10 20 30 40 50 60];
            for p = 1:length(find_y)
                isConvect2(p) = all(geom(connec(i,[j k(j)]),2) == find_y(p));
            end
            if any(isConvect2) && all(geom(connec(i,[j k(j)]),1) >= 10)
                convect(i,j) = 1;
            end
            
            % Find x = 10 & y is not inside block
            find_y = [0 10; 20 30; 40 50; 60 70];
            for p = 1:4
                isConvect3(p) = all([geom(connec(i,[j k(j)]),2) >= find_y(p,1);
                                    geom(connec(i,[j k(j)]),2) <= find_y(p,2)]);
            end
            if any(isConvect3) && all(geom(connec(i,[j k(j)]),1) == 10)
                convect(i,j) = 1;
            end 
        end
    end


%% Pre-processing
% Form the material property matrix
dee = [Kxx    0 ;
         0  Kyy];
     
% Counting of the free degrees of freedom
n=0;
for i=1:nnd
    for j=1:nodof
        if nf(i,j) ~= 0
            n=n+1;
            nf(i,j)=n;
        end
    end
end

%% Matrix assembly
% Assemble the global stiffness matrix and convection force vector
kk = zeros(n,n); %
fg = zeros(n,nodof);
for i=1:nel
    [bee,cee,fee,g,A] = elem_T3_heat(i); % Form stiffness matrix, force
                                         % vector, and steering vector
    
    % Stiffness matrix                                     
    kc = thick*A*bee'*dee*bee;      % Compute conduction matrix
    kh = thick*h*cee;               % Compute convection matrix
    ke = kc + kh;                   % Compute local stiffness matrix
    kk = form_kk(kk,ke,g);          % Assemble global stiffness matrix
    
    % Force vector
    fe = thick*h*T_inf*fee;
    fg = form_fg(fg,fe,g);
end

% Modify the global force vector due to boundaries
idT = find(Nodal_loads~=0);
for i=1:nnd
    if ~any(idT==i)
        fg(i) = fg(i) - sum(kk(i,idT)*Nodal_loads(idT));
    end
end
for i=1:nnd
    if Nodal_loads(i) ~= 0
        fg(i) = Nodal_loads(i); % Given temperature source
    end
end
[kk(idT,:),kk(:,idT)] = deal(0);
for i=1:length(idT)
    kk(idT(i),idT(i)) = deal(1);
end

%% Solve Problem
tg = kk\fg ; % solve for unknown displacements

figure(1)
    patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',tg, ...
    'Facecolor','interp','Marker','.')
    colorbar;
    set(gcf,'position',[400 500 500 400])
    title('Temperature Profile ^{o}C')
    xlabel('X (mm)')
    ylabel('Y (mm)')    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions List
function T3_mesh_LoganExample
    global geom connec nnd nel
    geom = [0 0; 2 0; 2 2; 0 2; 1 1];
    connec = [1 2 5; 1 5 4; 4 5 3; 2 3 5];
    nnd = 5;
    nel = 4;
end

function T3_mesh_final_term_unaxisymmetric
global nnd NXE NYE dhx dhy X_origin Y_origin geom connec nel Block
    k = 0;
    for b = 1:Block
        for i = 1:NXE(b)
            for j=1:NYE(b)
                k = k + 1;
                if b == 1
                    n1 = j + (i-1)*(NYE(b)+1);
                    n2 = j + i*(NYE(b)+1);
                    n3 = n1 + 1;
                    n4 = n2 + 1;
                    geom(n1,:) = [(i-1)*dhx(b)+X_origin(b)	(j-1)*dhy(b)+Y_origin(b)];
                    geom(n3,:) = [(i-1)*dhx(b)+X_origin(b)	j*dhy(b)+Y_origin(b)    ];
                else
                    if i == 1
                        n1 = find(and(geom(:,1)==X_origin(b), ...
                             round(geom(:,2),2)==round(Y_origin(b)+(j-1)*dhy(b),2)));
                        n2 = nnd + j;
                        n3 = n1 + 1;
                        n4 = n2 + 1;
                    else
                        n1 = nnd + j + (i-2)*(NYE(b)+1);
                        n2 = nnd + j + (i-1)*(NYE(b)+1);
                        n3 = n1 + 1;
                        n4 = n2 + 1;
                    end
                end
                nnd_b(b) = n4;
                geom(n2,:) = [i*dhx(b)+X_origin(b)      (j-1)*dhy(b)+Y_origin(b)];
                geom(n4,:) = [i*dhx(b)+X_origin(b)   	j*dhy(b)+Y_origin(b)    ];
                nel = 2*k;
                m = nel -1;
                connec(m,:) = [n1 n2 n3];
                connec(nel,:) = [n2 n4 n3];
            end
        end
        nnd = n4;
    end
end

function T3_mesh_final_term_axisymmetric
global nnd NXE NYE dhx dhy X_origin Y_origin geom connec nel Block
    k = 0;
    for b = 1:Block
        for i = 1:NXE(b)
            for j=1:NYE(b)
                k = k + 1;
                if b == 1
                    n1 = j + (i-1)*(2*NYE(b)+1);
                    n2 = j + i*(2*NYE(b)+1);
                    n3 = n1 + 1;
                    n4 = n2 + 1;
                    n5 = j + (2*i-1)*NYE(b)+i;
                    geom(n1,:) = [(i-1)*dhx(b)+X_origin(b)	(j-1)*dhy(b)+Y_origin(b)];
                    geom(n3,:) = [(i-1)*dhx(b)+X_origin(b)	j*dhy(b)+Y_origin(b)    ];
                else
                    if i == 1
                        n1 = find(and(geom(:,1)==X_origin(b), ...
                             round(geom(:,2),2)==round(Y_origin(b)+(j-1)*dhy(b),2)));
                        n2 = nnd + j + NYE(b);
                        n3 = n1 + 1;
                        n4 = n2 + 1;
                        n5 = nnd + j;
                    else
                        n1 = nnd + j + (i-2)*(2*NYE(b)+1) + NYE(b);
                        n2 = nnd + j + (i-1)*(2*NYE(b)+1) + NYE(b);
                        n3 = n1 + 1;
                        n4 = n2 + 1;
                        n5 = nnd + j + (2*i-3)*NYE(b) + i-1 + NYE(b);
                    end
                end
                nnd_b(b) = n4;
                geom(n2,:) = [i*dhx(b)+X_origin(b)      (j-1)*dhy(b)+Y_origin(b)];
                geom(n4,:) = [i*dhx(b)+X_origin(b)   	j*dhy(b)+Y_origin(b)    ];
                geom(n5,:) = mean([geom(n1,:); geom(n4,:)]);
                nel = 4*k;
                m = nel -3;
                connec(m,:) = [n1 n2 n5];
                m = nel -2;
                connec(m,:) = [n2 n4 n5];
                m = nel -1;
                connec(m,:) = [n4 n3 n5];
                connec(nel,:) = [n3 n1 n5];
            end
        end
        nnd = n4;
    end
end

% This function returns the coordinates of the nodes of element i
% and its steering vector
function [bee,cee,fee,g,A] = elem_T3_heat(i)
global nne nodof geom connec nf convect

    x1 = geom(connec(i,1),1); y1 = geom(connec(i,1),2);
    x2 = geom(connec(i,2),1); y2 = geom(connec(i,2),2);
    x3 = geom(connec(i,3),1); y3 = geom(connec(i,3),2);
    
    A = (0.5)*det([1 x1 y1; ...
                   1 x2 y2; ...
                   1 x3 y3]);
               
    L = [((x2-x1)^2 + (y2-y1)^2)^.5 ...
         ((x3-x2)^2 + (y3-y2)^2)^.5 ...
         ((x1-x3)^2 + (y1-y3)^2)^.5];
    
    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);
    
    % Conduction stiffness matrix
    bee = [m12 m22 m32;
           m13 m23 m33];
    
    % Convection stiffness matrix
    cee = 1/6*(convect(i,1)*L(1)*[2 1 0; 1 2 0; 0 0 0] + ...
               convect(i,2)*L(2)*[0 0 0; 0 2 1; 0 1 2] + ...
               convect(i,3)*L(3)*[2 0 1; 0 0 0; 1 0 2]);  
  
    % Heat force vector
    fee = 1/2*(convect(i,1)*L(1)*[1 1 0]' + ...
               convect(i,2)*L(2)*[0 1 1]' + ...
               convect(i,3)*L(3)*[1 0 1]');
    
    % Steering vector
    l=0;
    for k=1:nne
        for j=1:nodof
            l=l+1;
            g(l)=nf(connec(i,k),j);
        end
    end
end

% This function assembles the global stiffness matrix
function kk = form_kk(kk,ke,g)
global eldof
    for i=1:eldof
        if g(i) ~= 0
            for j=1:eldof
                if g(j) ~= 0
                    kk(g(i),g(j))= kk(g(i),g(j)) + ke(i,j);
                end
            end
        end
    end
end

% This function assembles the global force vector
function fg = form_fg(fg,fe,g)
global eldof
    for i=1:eldof
        if g(i) ~= 0
            fg(g(i))= fg(g(i)) + fe(i);
        end
    end
end