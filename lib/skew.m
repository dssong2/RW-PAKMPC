function s = skew(a)
%SKEW  (a)^x, so that skew(a)*b = cross(a, b).
s = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
