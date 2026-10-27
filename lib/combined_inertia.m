function [jc, jt] = combined_inertia(p, s)
%COMBINED_INERTIA  Combined inertia with the target's mass and inertia jointly scaled by s, Sec. V-D.
%   [jc, jt] = combined_inertia(p, s) returns Jc(s) about the combined center of mass and the
%   scaled target inertia jt = s Jt. s = 1 gives p.jc exactly.
%   [guess] The paper gives Jc but not how it is assembled. The geometry follows
%   build_pakmpc_multibody: servicer cube (side 1.5 m), target box matching Jt, docked along the
%   servicer y axis with a 0.2 m ring, and the rest of Jc (arms, panels, ring) lumped at the
%   nominal combined center of mass and held fixed while the target is scaled.

ls = sqrt(6*p.js(1, 1)/p.ms);
lt = sqrt(6/p.mt*(trace(p.jt) - 2*diag(p.jt)'));
r = [0; ls/2 + 0.2 + lt(2)/2; 0];
pa = @(m, d) m*((d'*d)*eye(3) - d*d');
parts = @(mt, jt) p.js + pa(p.ms, -mt*r/(p.ms + mt)) + jt + pa(mt, r - mt*r/(p.ms + mt));
jres = p.jc - parts(p.mt, p.jt);
jt = s*p.jt;
jc = jres + parts(s*p.mt, jt);
end
