function [p4, signk3, k123, pb] = trilateration( p, l )
%  trilateration(p,l) Performs trilateration between 3 points.
%     TRILATERATION(P,L) Performs trilateration between 3 points P1, P2 and
%     P3 to find a fourth P4 point lying in the intersection of 3 spheres
%     with centers at the three mentioned points with radius L1, L2 and L3.
%     
%     P4=TRILATERATION(P,L) Returns the fourth point P4.
% 
%     [P4,SIGNK3]=TRILATERATION(P,L) Also returns wheter P1, P2 and P3 are
%     deployed in CCW rotation or CW as seen from P4.
%
%     TRILATERATION(P) is used to get either the equation for P4 or to
%     verify that the trilateration was performed correctly.
%        a) If three points are provided the equation for P4 is returned as
%           function of x, y and z, P4=f(x,y,z).
%        b) If four points are provided P4 is verified.
% 
%     INPUTS
%         P   A matrix containing points as vectors:
%                 P=[P1 P2 P3]
%                 P=[P1 P2 P3 P4]
%         L   A row vector with the distances from the centers of three
%             points to a fourth unknown position
%                 L=[L1 L2 L3].
% 
%     OUTPUTS
%         LAR FIXME The equation for is no longer available because it was
%         very slow.
%         P4       P4, either as an equation function or numeric value.
%         SIGNK3   States if the points P1, P2 and P3 are deployed in CCW
%                  or CC rotation as seen from P4. Not available for
%                  verification or equation form. TODO make it available.
%                      0  Verify algorithm
%                      1  CCW rotation
%                     -1  CC  rotation
%         K        A vector, values of k1 k2 and k3, k3 already has sign.
%         PB       The projection of P4 on the base plane
% 
%     Based on the paper:
% 
%     Revisiting Trilateration for Robot Localization, Federico Thomas and
%     Lluis Ros. 2005

%% Validate instruction
if size(p,1)~=3
    error(['myApp:argChk', 'p must be declared in 3 dimensional space,'...
           ' it has the positions of the vectors as p = [p1 p2 p3],'...
           ' p4 may be included']);
end
if size(p,2)<3 || size(p,2)>4
    error('myApp:argChk', 'Only 3 or 4 points are required');
end

%% Base plane

% Points
p1 = p(:,1);
p2 = p(:,2);
p3 = p(:,3);

Ab4 = caley_menger( [p1 p2 p3] );
% check_err = 1e-3;
if Ab4 == 0 % No sense to keep going
    error( 'myApp:argChk', 'p1, p2 and p3 are aligned, p4 is undefined' );
% elseif abs( Ab4 ) < check_err
%     display( 'There may be a huge error' );
end

% Vectors, section III in paper
v12 = p2 - p1; % From p1 to p2
v13 = p3 - p1; % From p1 to p3
crossv12v13 = cross( v12, v13 );

%% Trilateration algorithm

signk3 = 0;

if nargin==1
    %% Prove trilateration
    %
    % If no distances are provided, trilateration will not be performed, it
    % means that either the prove for a solution has been requested or the
    % equation to get the fourth point is needed
    
    if size(p,2)==4 % Prove algorithm
        p4 = p(:,4);
    elseif size(p,2)==3 % Get equation
        syms x y z
        p4 = [x y z].';
    end
    
    % p4 has two symmetric locations, hence this option gives two solutions

    p4_1 = p1 + 1/caley_menger([p1 p2 p3])*...
         (-caley_menger([p1 p2 p3],[p1 p3 p4])*v12...
          +caley_menger([p1 p2 p3],[p1 p2 p4])*v13...
          +sqrt(caley_menger([p1 p2 p3 p4]))*crossv12v13);

    p4_2 = p1 + 1/caley_menger([p1 p2 p3])*...
         (-caley_menger([p1 p2 p3],[p1 p3 p4])*v12...
          +caley_menger([p1 p2 p3],[p1 p2 p4])*v13...
          -sqrt(caley_menger([p1 p2 p3 p4]))*crossv12v13);

    p4 = [p4_1 p4_2];

else
    %% Execute Trilateration
    %
    % If distances are provided, perform trilateration, p4 MUST NOT be used
    % explicity in any of the calculations, otherwise what are looking for?

    % Validate instruction
    if size(p,2)<size(l,2) || size(l,1)~=1 || size(l,2)~=3
        error('myApp:argChk',...
              ['p and l do not have the same number of elements or'...
               ' the distances are not defined as scalars. Only'...
               ' three points are needed.']);
    end
    
    % Assign distances
    l14 = l(1);
    l24 = l(2);
    l34 = l(3);

    k1 = ...
  -caley_menger([p1 p2 p3 p1 p3],[l14 l24 l34]) / caley_menger([p1 p2 p3]);
    k2 = ...
   caley_menger([p1 p2 p3 p1 p2],[l14 l24 l34]) / caley_menger([p1 p2 p3]);
       
    % Vector normal to the base
    base_normal = crossv12v13;
    base_normal = base_normal/norm(base_normal); % Normalize vector

    n_err = 3e-3; % The operations do not give exact values, and maybe with
                  % noise and real sensor values this has to be bigger.

    % Equation, mm, between (11) and (12) in paper
    k3 = sqrt( abs(...
         caley_menger([p1 p2 p3],[l14 l24 l34])))/caley_menger([p1 p2 p3]);

    % k3 is the height divided by the norm of ( v12 x v13 ), the
    % normal vector of the base plane, that means the term in equation (9)
    % "k3(v12xv13)" needed for p4 is nothing but the height with direction.
    % So, if p1, p2 and p3 are deployed in CCW fashion when they are seen
    % from p4, chose the "+" sign because that means k3 is in the direction
    % of the normal vector to the base, otherwise is negative and the "-"
    % sign must be chosen.

    % TODO For this approach to work it is assumed that p4 is always above
    % the base plane in the real world. If that were not the case then
    % wheter p4 is above or below the base plane MUST be passed as an
    % argument and the if condition changes to p1(3)-signk3(3)>=n_err when
    % p4 is below the base plane. Abstractly speaking, from the point of
    % view of the base normal p4 can be below the base plane when it is not
    % in the real world.

    % Add the normal to any base point to get an offset, up or down.
    signk3 = p1 + base_normal;
    if p1(3)-signk3(3)<=n_err % CCW
        signk3=1;
    else % CC
        signk3=-1;
    end
    k3 = signk3 * k3;

    % Compute the solution for the projection of p4 on the base, part of
    % the derivation of equation (9)
    pb = p1 + k1*v12 + k2*v13;

    % Trilateration solution p4
    % Equation (9), euqation(12) is also the same, but it doesn't
    % calculate k1, k2 and k3
    p4 = pb + k3 * crossv12v13;

    % k vector
    k123 = [k1 k2 k3].';
end

end
