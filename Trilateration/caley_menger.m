function D = caley_menger(p,q)
%  caley_menger(p,q) Calculates the Caley Menger determinant.
%     CALEY_MENGER(P,Q) Calculates the Caley Menger determinant between two
%     sets of points P and Q, where the columns indicate the number of
%     points.
%     
%     D=CALEY_MENGER(P,Q) return the value of the Caley Menger determinant
%     in the variable D
% 
%     CALEY_MENGER(P) makes Q=P
% 
%     Based on the paper:
% 
%     Revisiting Trilateration for Robot Localization, Federico Thomas and
%     Lluis Ros. 2005
% 
%     Some of the equations are from the appendix.

%% Validate Input

if nargin==1 % Same vector is used
    q = p;
end

if size(p,1)==3 && size(p,2)==size(q,2) && size(q,1)==1
    l = q;
elseif size(p,1)==3 && size(p,2)==5 && size(q,2)==3 && size(q,1)==1
    l = q;
elseif size(p)~=size(q)
    error('myApp:argChk', 'p and q do not have the same number of points');
end

n = size(p,2); % Number of vectors
if n==1
    error('myApp:argChk', 'At least 2 points are required');
end

%% Algorithm

if ( size(p,1)==3 && size(p,2)==size(q,2) && size(q,1)==1 ) || ...
   ( size(p,1)==3 && size(p,2)==5 && size(q,2)==3 && size(q,1)==1 )

    if n==2 % This has to be equal to 4A^2
       
        p1 = p(:,1);
        p2 = p(:,2);

        a = norm( p2 - p1 );
        b = l(1);
        c = l(2);
        
        % From wikipedia, the area of a triangle knowing the length of the
        % sides is:
        %
        % A = 1/4 sqrt((a+b+c)(-a+b+c)(a-b+c)(a+b-c))
        %
        % So:
        %
        % 4*A^2 = 1/4(a+b+c)(-a+b+c)(a-b+c)(a+b-c)
        % 
        % So, the paper is wrong!!!! :/

        D = (   a + b + c )*...
            ( - a + b + c )*...
            (   a - b + c )*...
            (   a + b - c )/4;
        
    elseif n==3 % This has to be equal to 36V^2

        p1 = p(:,1);
        p2 = p(:,2);
        p3 = p(:,3);

        l1 = l(1);
        l2 = l(2);
        l3 = l(3);

        a = norm( p2 - p1 );
        b = norm( p3 - p1 );
        c = norm( p3 - p2 );

        X = ( l1 - a + l2 ) * ( a + l1 + l2 );
        Y = ( l3 - b + l1 ) * ( b + l1 + l3 );
        Z = ( l2 - c + l3 ) * ( c + l3 + l2 );

        x = ( a - l2 + l1 ) * ( l2 - l1 + a );
        y = ( b - l1 + l3 ) * ( l1 - l3 + b );
        z = ( c - l3 + l2 ) * ( l3 - l2 + c );

        XI_var     = sqrt(x*Y*Z);
        eta_var    = sqrt(y*Z*X);
        xi_var     = sqrt(z*X*Y);
        lambda_var = sqrt(x*y*z);

        D = 1 / ( 1024 * l1^2 * l2^2 * l3^2 ) * ...
            (  XI_var + eta_var + xi_var - lambda_var ) * ...
            (  XI_var + eta_var - xi_var + lambda_var ) * ...
            (  XI_var - eta_var + xi_var + lambda_var ) * ...
            ( -XI_var + eta_var + xi_var + lambda_var );
        
%         D = 1 / 8 * det( ...
% [0                     1                     1                     1    1;
%  1                     0 caley_menger([p1,p2]) caley_menger([p1,p3]) l1^2;
%  1 caley_menger([p1,p2])                     0 caley_menger([p2,p3]) l2^2;
%  1 caley_menger([p1,p3]) caley_menger([p2,p3])                     0 l3^2;
%  1                  l1^2                  l2^2                  l3^2  0] );

    elseif n==5

        % i
        p1 = p(:,1);
        p2 = p(:,2);
        p3 = p(:,3);
        % j
        p1_ = p(:,4);
        p2_ = p(:,5);

        l1 = l(1);
        l2 = l(2);
        l3 = l(3);

        D = - 1 / 4 * det( ...
                  [0                    1                    1    1;
                   1 caley_menger([p1,p1_]) caley_menger([p1,p2_]) l1^2;
                   1 caley_menger([p2,p1_]) caley_menger([p2,p2_]) l2^2;
                   1 caley_menger([p3,p1_]) caley_menger([p3,p2_]) l3^2] );
    else
        error('myApp:argChk', 'No implementation for more than 3 points');
    end

else
    
    % If n=2 and p=q the formula reduces to the squared Euclidean distance
    if n==2 && nargin==1
        D = norm(p(:,2)-p(:,1))^2;
        return
    end

    % If symbolic variables are included
    %D = sym(zeros(n));
    D = zeros(n);

    for i=1:n
        for j=1:n
            D(i,j) = norm(q(:,j)-p(:,i))^2;
        end
    end

    D = [0                 ones(1,size(D,2));
         ones(size(D,1),1) D];

    D = 2*(-1/2)^n*det(D);

    %if size(symvar(D),2)==0
    %    D = double(D);
    %end

end

end
