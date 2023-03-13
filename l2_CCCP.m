function [x_est_mtx, fail, Y_BFGS] = l2_CCCP(Tx, Rx, Rg, b, Nmax, epsilon)
global Y;
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','OutputFcn',@outfun);
%Concave-Convex Procedure (CCCP) for l2-minimization

%Inputs:
%Tx - transmitter position matrix H \times M
%Rx - receiver position matrix H \times L
%Rg - M \times L TSOA-based range sum observation matrix
%Nmax - user-defined maximum number of CCCP iterations
%epsilon - threshold that defines the termination criterion
%b - balancing parameter

%Outputs:
%x_est_mtx - location estimate matrix H \times Ncccp
%fail - indicator for the feasibility of the run

fail = false;

Y_BFGS = zeros(4,6,200);

[H, M] = size(Tx);
[~, L] = size(Rx);

w_mtx = zeros(M,L);

for mm = 1:M
    for ll = 1:L
        w_mtx(mm,ll) = 1 - (Rg(mm,ll)/sum(sum(Rg)));
    end
end

cnt = 0;

x_ini = zeros(H,1);

x_est_mtx = [];

x_est_mtx = [x_est_mtx,x_ini];

x_kth = x_ini;

%main loop
while cnt <= Nmax

    cnt = cnt + 1;
    
    [y_est,fval,exitflag,output] = fminunc(@obj_fun,x_kth,options);

    x = y_est;

    [sizex, sizey] = size(Y);

    Y_BFGS(1:sizex,1:sizey,cnt) = Y;
    
    if (isnan(sum(x))) || (sum(x) == +Inf) || (sum(x) == -Inf)
        fail = true;
        return
    end
    
    x_est_mtx = [x_est_mtx,x];
    
    x_kth = x;
    
    if (norm(x_est_mtx(:,end) - x_est_mtx(:,end-1))/(min(norm(x_est_mtx(:,end)), norm(x_est_mtx(:,end-1)))) <= epsilon)
        return
    end
    
end



    function obj = obj_fun(x_vec)
        
        obj = 0;
        
        for m = 1:M
            for l = 1:L
                
                obj = obj + w_mtx(m,l)*( Rg(m,l)^2 + norm(x_vec(1:H) - Tx(:,m))^2 + norm(x_vec(1:H) - Rx(:,l))^2 + b^2 ...
                    + 2*norm(x_vec(1:H) - Tx(:,m))*norm(x_vec(1:H) - Rx(:,l)) + (x_vec(1:H) - Tx(:,m))'*(x_vec(1:H) - Rx(:,l)) ...
                    + 2*norm(x_vec(1:H) - Tx(:,m))*b + 2*norm(x_vec(1:H) - Rx(:,l))*b ...
                    ...
                    - 2*Rg(m,l)*norm(x_kth(1:H) - Tx(:,m)) - 2*Rg(m,l)*norm(x_kth(1:H) - Rx(:,l)) - 2*Rg(m,l)*b ...
                    - (x_kth(1:H) - Tx(:,m))'*(x_kth(1:H) - Rx(:,l)) ...
                    ...
                    - ( [(2*Rg(m,l)*(x_kth(1:H)-Tx(:,m))/norm(x_kth(1:H)-Tx(:,m)))']*(x_vec - x_kth) + [(2*Rg(m,l)*(x_kth(1:H)-Rx(:,l))/norm(x_kth(1:H)-Rx(:,l)))']*(x_vec - x_kth) ...
                       + [(2*x_kth(1:H) - (Tx(:,m)+Rx(:,l)))']*(x_vec - x_kth) ) );
                
            end
        end
        
        
    end



end

