NDOS[m_, n_] := Module[{x, xp, prec = Floor[1.5 n m Log[2]/Log[10]], b, a, c2, s2, c0, s0, cn, sn},
  b = 2 x (1 - x^2);
  a[k_] := (1 + x^2)^2 - b Cos[Pi k/n];
  c0 = (1 - x)^m + x^m (1 + x)^m;
  s0 = (1 - x)^m - x^m (1 + x)^m;
  cn = (1 + x)^m + x^m (1 - x)^m;
  sn = (1 + x)^m - x^m (1 - x)^m;
  c2[k_] := (Sum[ m!/(2 j)!/(m - 2 j)! (a[k]^2 - b^2)^j a[k]^(m - 2 j), {j, 0, IntegerPart[m/2]}] + b^m)/2^(m - 1);
  s2[k_] := (Sum[ m!/(2 j)!/(m - 2 j)! (a[k]^2 - b^2)^j a[k]^(m - 2 j), {j, 0, IntegerPart[m/2]}] - b^m)/2^(m - 1);
  xp = Expand[ N[ (1/2) If[Mod[n, 2] == 0,
    Product[c2[2 k + 1], {k, 0, n/2 - 1}]
  + Product[s2[2 k + 1], {k, 0, n/2 - 1}]
  + c0 cn Product[c2[2 k], {k, 1, n/2 - 1}]
  + s0 sn Product[s2[2 k], {k, 1, n/2 - 1}],
    cn Product[c2[2 k + 1], {k, 0, (n - 3)/2}]
  + sn Product[s2[2 k + 1], {k, 0, (n - 3)/2}]
  + c0 Product[c2[2 k], {k, 1, (n - 1)/2}]
  + s0 Product[s2[2 k], {k, 1, (n - 1)/2}]], prec]];
  Take[Round[Chop[CoefficientList[xp, x]]], {1, -1, 2}]];

savels[fn_, ls_] := Module[{fp = OpenWrite[fn], i},
  For[i = 1, i <= Length[ls], i++, Write[fp, ls[[i]]]]; Close[fp]];

easydos[n_, m_] := Module[{dos = NDOS[n, m], logdos = Table[0, {n m + 1}], i},
  savels["IsingDOS" <> ToString[n] <> "x" <> ToString[m] <> ".txt", dos];
  For[i = 1, i <= n m + 1, i++,
    logdos[[i]] = If[dos[[i]] == 0, -10000, N[Log[dos[[i]]], 17]]];
  savels["islogdos" <> ToString[n] <> "x" <> ToString[m] <> ".txt", logdos]];
