function r = q2r(q)
%Q2R  Body-to-inertial rotation matrix of a scalar-first quaternion, R' = R w^x.
s = q(1); v = q(2:4);
r = (s^2 - v'*v)*eye(3) + 2*(v*v') + 2*s*skew(v);
end
