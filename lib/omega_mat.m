function o = omega_mat(w)
%OMEGA_MAT  Omega(w) of Eq. (1), q' = Omega(w) q with q scalar first.
o = 0.5*[0    -w(1) -w(2) -w(3);
         w(1)  0     w(3) -w(2);
         w(2) -w(3)  0     w(1);
         w(3)  w(2) -w(1)  0];
end
