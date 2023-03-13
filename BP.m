function [x_est] = BP(Tx, Rx, Rg, Nmax, epsilon)
%Code for Error-Reduced Elliptic Positioning via Joint Estimation of Location and a Balancing Parameter
%implemented by Wenxin Xiong
%Feel free to drop me an email if you have any question
%w.x.xiong@outlook.com

%Inputs:
%Tx - transmitter position matrix H \times M
%Rx - receiver position matrix H \times L
%Rg - M \times L TSOA-based range sum observation matrix
%Nmax - user-defined maximum number of CCCP iterations
%epsilon - threshold that defines the termination criterion

%Outputs:
%x_est - location estimate


[~, M] = size(Tx);
[~, L] = size(Rx);

[x_est_mtx_CCCP, ~] = l2_CCCP(Tx, Rx, Rg, 0, Nmax, epsilon);
b_est = 0;
for m = 1:M
    for l = 1:L
        b_est = b_est + (Rg(m,l)-norm(x_est_mtx_CCCP(:,end)-Tx(:,m))-norm(x_est_mtx_CCCP(:,end)-Rx(:,l)))/(M*L);
    end
end
[x_est_mtx_CCCP, ~] = l2_CCCP(Tx, Rx, Rg, b_est, Nmax, epsilon);
x_est = x_est_mtx_CCCP(:,end);

end

