function [Rs_est, a_est] = GEVD(Rx, Rn)
[U,S]=eig(Rx, Rn, 'vector');
[MAX, Index]=max(S);
 U=(U')^(-1);
u1 = U(:, Index);

Rs_est = u1*(MAX-1)*u1';
[Q, Lam, ~]=svd(Rs_est+eps);
a_est = Q(:,1)*sqrt(Lam(1,1)/Rs_est(1,1));
end