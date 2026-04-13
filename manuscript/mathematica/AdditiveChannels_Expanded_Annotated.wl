(* Additive Channels — expanded Wolfram Language companion *)

ClearAll["Global`*"];

(* 2-trait worked example *)
G2 = {{1.0, 0.3}, {0.3, 0.6}};
b2 = {0.8, 0.4};
H2 = {{-0.4, 0}, {0, -0.4}};

Vlin2 = b2.G2.b2;
Vquad2 = 1/2 Tr[(H2.G2).(H2.G2)];
Ag2 = N[Vlin2/(Vlin2 + Vquad2)];

Print[{"Vlin" -> Vlin2, "Vquad" -> Vquad2, "Ag" -> Ag2}];

(* Logit identity *)
ClearAll[VlinT, VquadT, t];
AgT = VlinT[t]/(VlinT[t] + VquadT[t]);
LlogitT = Log[VlinT[t]/VquadT[t]];
rT = D[LlogitT, t];

Print["Identity check:"];
Print[FullSimplify[D[AgT, t] - AgT (1 - AgT) rT]];

(* One-trait moving optimum *)
ClearAll[s, d, Gsc, thetaDot, lambdaRec, Gstar, m, Ne];
Vlin1D = (s d)^2 Gsc;
Vquad1D = (1/2) s^2 Gsc^2;
Ag1D = FullSimplify[Vlin1D/(Vlin1D + Vquad1D)];
r1D = FullSimplify[2 thetaDot/d - s Gsc + lambdaRec (1 - Gstar/Gsc) - m/Gsc + 1/(2 Ne)];

Print[{"Ag1D" -> Ag1D, "r1D" -> r1D}];

(* Two-trait gradient rotation *)
ClearAll[lambda1, lambda2, gamma, theta, bnorm];
VlinRot = bnorm^2 (lambda1 Cos[theta]^2 + lambda2 Sin[theta]^2);
VquadRot = 1/2 gamma^2 (lambda1^2 + lambda2^2);
AgRot = FullSimplify[VlinRot/(VlinRot + VquadRot)];

Print["AgRot(theta):"];
Print[AgRot];