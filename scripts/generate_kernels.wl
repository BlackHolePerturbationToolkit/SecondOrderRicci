(* ::Package:: *)
(* Generate the C++ mode-coupling kernels from the PerturbationEquations package.

   Same operator and C conversion rules as notebooks/GenerateCouplingCode.nb, but every
   component d2R_i, which is linear in the coupling integrals CInt (the only place where
   l3, m3 enter), is split by coupling:
       d2R_i = sum_k c_k(l3, m3, l1, m1, l2, m2) G_ik(r; l1, m1, l2, m2),
   and written as
     src/R2G_<i>.cc : R2G_<i>(...) adds G_ik(r_j) to G[k][j] for jb <= j < je, taking the
                      first-order modes as ModeBlocks (include/modeblock.h);
                      R2C_<i>(l3, m3, l1, m1, l2, m2, c) evaluates the c_k;
     include/R2G.h  : declarations, R2_NC[i] (number of couplings), R2_NC_MAX.
   A term with coupling spins (s3; s2, s1) requires l1 >= |s1|, l2 >= |s2| and every
   mu[l, n] in its denominator to be non-zero (l >= n); otherwise it is skipped.
   The output is scaled (in R2C) to keep the normalisation of the original
   (pre-PerturbationEquations) coupling code: d2R from PerturbationEquations 0.3.0 is 2x the
   original, and 0.4.0 normalises it 4x smaller than 0.3.0 (verified numerically against the
   original code: all components and modes agree to ~1e-14 after scaling).

   Code forms (optional argument, default "real"):
     direct   : each G_ik written as the (nested) expression from PerturbationEquations;
     factored : each G_ik expanded as  sum_{p,q} r^p f^q sum_{a,b} K_{k,a,b,p,q} X_a Y_b,  with
                X_a, Y_b the fields of the two modes (the operator is bilinear in them), the
                constants K (functions of l, m only) computed once per call, and each distinct
                product X_a Y_b and power r^p f^q computed once per grid point.

     real     : as factored, but written in explicit double arithmetic on the separate real
                and imaginary planes of the fields (ModeBlock) and of the outputs G_k,
                branch-free (constants of guarded-off couplings are zero), so that the
                grid-point loop is vectorised across grid points.
   Only the real form matches the current storage; the other two are kept in the git history.

   Usage: WolframKernel -script scripts/generate_kernels.wl [direct|factored|real] *)

root = DirectoryName[$InputFileName] // ParentDirectory;
form = SelectFirst[Join[$ScriptCommandLine, $CommandLine], MemberQ[{"direct", "factored", "real"}, #] &, "real"];
Print["code form: ", form];
If[form =!= "real", Print["only the real form matches the current storage (separate real/imaginary planes); the direct and complex-factored forms are in the git history (commit 0bed18f)"]; Abort[]];

Block[{Print}, << xAct`PerturbationEquations`];
Print["PerturbationEquations ", xAct`PerturbationEquations`$Version];
(* overall factor applied to every coupling coefficient, to match the original normalisation *)
scale = Switch[First[xAct`PerturbationEquations`$Version], "0.3.0", "0.5", "0.4.0", "2.0",
  _, Print["Unverified PerturbationEquations version: check the normalisation against the original code and add it here."]; Abort[]];
Print["output scale: ", scale];
d2R = SchwarzschildQuadraticOperator["d2R", "Lorenz", "BLS", "BLS"];
d2R = Map[Collect[# //. FrequencyDomainConversion, CInt[__], Simplify] &, d2R];

(* --- conversion to C, as in GenerateCouplingCode.nb ------------------------------- *)
hRules = Sequence[PD[{1, -BL}][PD[{1, -BL}][hBS[LI[i_], LI[l1], LI[m1]]]] :> ddh1[i], PD[{1, -BL}][hBS[LI[i_], LI[l1], LI[m1]]] :> dh1[i], hBS[LI[i_], LI[l1], LI[m1]] :> h1[i], PD[{1, -BL}][PD[{1, -BL}][hBS[LI[i_], LI[l2], LI[m2]]]] :> ddh2[i], PD[{1, -BL}][hBS[LI[i_], LI[l2], LI[m2]]] :> dh2[i], hBS[LI[i_], LI[l2], LI[m2]] :> h2[i]];
LRules = Sequence[mu[LI[l1], -LI[1]] -> mu1l1, mu[LI[l1], -LI[2]] -> mu2l1, mu[LI[l2], -LI[1]] -> mu1l2, mu[LI[l2], -LI[2]] -> mu2l2, mu[LI[l1], -LI[3]] -> mu3l1, mu[LI[l2], -LI[3]] -> mu3l2];
gridRules = Sequence[r[] -> r, f[r[]] -> f];
toC[expr_] := StringReplace[StringReplace[ToString[CForm[expr /. {\[Omega]1 -> omega1, \[Omega]2 -> omega2} /. {LRules} /. {gridRules} /. {x : (omega1 + omega2)^_ :> Expand[x]} /. {hRules, p : Power[(r | f | M | omega1 | omega2 | mu1l1 | mu1l2 | mu2l1 | mu2l2 | mu3l1 | mu3l2), _?(# != -1 &)] :> p, i : _?NumericQ /; Abs[i] != 1 :> N[i]}]],
    {"Power(" ~~ x : ("r" | "f" | "M" | "omega1" | "omega2" | "mu1l1" | "mu1l2" | "mu2l1" | "mu2l2" | "mu3l1" | "mu3l2") ~~ "," ~~ n : NumberString ~~ ")" :> x <> "_" <> n,
     h : ("h1" | "dh1" | "ddh1" | "h2" | "dh2" | "ddh2") ~~ "(" ~~ i : NumberString ~~ ")" :> h <> "_" <> i,
     "Complex(0,1)" :> "I",
     "Complex(0," ~~ n : NumberString ~~ ")" :> ToString[N[ToExpression[n]], FormatType -> CForm] <> "*I",
     "Complex(0.," ~~ n : NumberString ~~ ")" :> ToString[N[ToExpression[n]], FormatType -> CForm] <> "*I"}],
   {"Power(" ~~ h : (("h1" | "dh1" | "ddh1" | "h2" | "dh2" | "ddh2") ~~ "_" ~~ NumberString) ~~ ",2)" :> h <> "*" <> h}];

spinOf[arg_] := First[Cases[arg, LI[x_] :> x, {0, Infinity}]];

(* --- factored form ------------------------------------------------------------------- *)
xyRules = {PD[{1, -BL}][PD[{1, -BL}][hBS[LI[i_], LI[l1], LI[m1]]]] :> X[2, i], PD[{1, -BL}][hBS[LI[i_], LI[l1], LI[m1]]] :> X[1, i], hBS[LI[i_], LI[l1], LI[m1]] :> X[0, i],
   PD[{1, -BL}][PD[{1, -BL}][hBS[LI[i_], LI[l2], LI[m2]]]] :> Y[2, i], PD[{1, -BL}][hBS[LI[i_], LI[l2], LI[m2]]] :> Y[1, i], hBS[LI[i_], LI[l2], LI[m2]] :> Y[0, i]};
fname[X[s_, i_]] := {"h1", "dh1", "ddh1"}[[s + 1]] <> "_" <> ToString[i];
fname[Y[s_, i_]] := {"h2", "dh2", "ddh2"}[[s + 1]] <> "_" <> ToString[i];
pname[x_, y_] := "P_" <> fname[x] <> "_" <> fname[y];
sgn[n_] := If[n < 0, "m" <> ToString[-n], ToString[n]];
rname[{p_, q_}] := "R_r" <> sgn[p] <> "_f" <> sgn[q];
powExpr[v_, n_] := Which[n == 0, "1.0", n > 0, StringRiffle[Table[v, {n}], "*"], True, StringRiffle[Table[v <> "inv", {-n}], "*"]];
rdecl[{p_, q_}] := "    const double " <> rname[{p, q}] <> " = " <> powExpr["r", p] <> "*" <> powExpr["f", q] <> ";\n";
(* expand one coefficient into {p, q, X, Y, constant} records *)
expandTerms[coef_] := Module[{e = Expand[coef /. xyRules /. {r[] -> rr, f[r[]] -> ff}], ts},
   ts = If[Head[e] === Plus, List @@ e, {e}];
   Map[Function[t, Module[{x = First[Cases[t, _X, {0, Infinity}]], y = First[Cases[t, _Y, {0, Infinity}]]},
       {Exponent[t, rr], Exponent[t, ff], x, y, t /. {_X -> 1, _Y -> 1, rr -> 1, ff -> 1}}]], ts]];

(* --- per component ------------------------------------------------------------------ *)
nc = <||>;
Do[
  expr = d2R[i];
  cints = SortBy[Union[Cases[expr, _CInt, Infinity]], spinOf /@ (List @@ #)[[{3, 6, 9}]] &];
  If[Simplify[expr - Total[(Coefficient[expr, #] # &) /@ cints]] =!= 0, Print["component ", i, " is not linear in CInt"]; Abort[]];
  terms = Table[
    Module[{coef = Coefficient[expr, c], s = spinOf /@ (List @@ c)[[{3, 6, 9}]], den, req1, req2, code},
      den = Denominator[Together[coef]];
      req1 = Max[Abs[s[[3]]], Cases[den, mu[LI[l1], -LI[n_]] :> n, {0, Infinity}]];
      req2 = Max[Abs[s[[2]]], Cases[den, mu[LI[l2], -LI[n_]] :> n, {0, Infinity}]];
      code = toC[coef];
      If[StringContainsQ[code, "CInt" | "LI" | "mu("], Print["unconverted symbol in component ", i]; Abort[]];
      <|"s" -> s, "req1" -> req1, "req2" -> req2, "code" -> code|>],
    {c, cints}];
  nc[i] = Length[terms];

  If[form === "factored" || form === "real",
    (* records per coupling, combined over identical (p, q, X, Y) *)
    recs = Table[Normal[GroupBy[expandTerms[Coefficient[expr, cints[[k]]]], (#[[1 ;; 4]] &) -> (#[[5]] &), Total]], {k, Length[cints]}];
    (* split every constant into real and imaginary parts (all symbols are real), so that
       each term is a real x complex multiplication; the I is applied once per sum *)
    reim[c_] := Module[{e = Expand[c]}, {e /. Complex[a_, b_] :> a, (e - (e /. Complex[a_, b_] :> a))/I // Expand}];
    recs = Table[Flatten[(Function[rule, MapThread[If[#2 === 0, Nothing, Append[rule[[1]], #1] -> #2] &, {{"re", "im"}, reim[rule[[2]]]}]] /@ recs[[k]])], {k, Length[recs]}];
    consts = DeleteDuplicates[Flatten[Table[Values[recs[[k]]], {k, Length[recs]}]]];
    consts = DeleteCases[consts, 0];
    cIndex = AssociationThread[consts -> Range[Length[consts]]];
    If[! AllTrue[consts, FreeQ[N[# /. {\[Omega]1 -> 1.1, \[Omega]2 -> 1.3, M -> 1.7, mu[__] -> 1.9}], _Complex] &], Print["complex constant left in component ", i]; Abort[]];
    cDecl = StringJoin[Table["  const double K" <> ToString[n] <> " = " <> toC[consts[[n]]] <> ";\n", {n, Length[consts]}]];
    rmons = Union[Flatten[Table[#[[1 ;; 2]] & /@ Keys[recs[[k]]], {k, Length[recs]}], 1]];
    prods = Union[Flatten[Table[#[[3 ;; 4]] & /@ Keys[recs[[k]]], {k, Length[recs]}], 1]];
    fieldsF = Union[Flatten[prods]];
    term[rule_] := "K" <> ToString[cIndex[rule[[2]]]] <> "*" <> pname @@ rule[[1, 3 ;; 4]];
    body = StringJoin[Table[
       Module[{byR = GroupBy[Normal[recs[[k]]], (#[[1, 1 ;; 2]] &)]},
        "    if (use" <> ToString[k - 1] <> ") G[" <> ToString[k - 1] <> "][j] +=\n      " <>
         StringRiffle[KeyValueMap[Function[{pq, rs},
             Module[{re = Select[rs, #[[1, 5]] === "re" &], im = Select[rs, #[[1, 5]] === "im" &], sre, sim},
              sre = If[re === {}, "", "(" <> StringRiffle[term /@ re, " + "] <> ")"];
              sim = If[im === {}, "", "I*(" <> StringRiffle[term /@ im, " + "] <> ")"];
              rname[pq] <> "*(" <> StringRiffle[DeleteCases[{sre, sim}, ""], " + "] <> ")"]], byR], "\n      + "] <> ";\n"],
       {k, Length[recs]}]];
    If[form === "real",
      kconsts = Table[DeleteCases[DeleteDuplicates[Values[recs[[k]]]], 0], {k, Length[recs]}];
      kIndex = Table[AssociationThread[kconsts[[k]] -> Range[Length[kconsts[[k]]]]], {k, Length[recs]}];
      cDecl = StringJoin[Table["  const double K" <> ToString[k - 1] <> "_" <> ToString[n] <> " = use" <> ToString[k - 1] <> " ? (" <> toC[kconsts[[k, n]]] <> ") : 0.0;\n",
          {k, Length[recs]}, {n, Length[kconsts[[k]]]}]];
      rterm[k_, rule_, part_] := "K" <> ToString[k - 1] <> "_" <> ToString[kIndex[[k]][rule[[2]]]] <> "*" <> pname @@ rule[[1, 3 ;; 4]] <> part;
      rsum[k_, rs_, part_] := If[rs === {}, "0.0", StringRiffle[rterm[k, #, part] & /@ rs, " + "]];
      body = StringJoin[Table[
         Module[{byR = GroupBy[Normal[recs[[k]]], (#[[1, 1 ;; 2]] &)], gr, gi},
          "    {\n      double gr = 0.0, gi = 0.0;\n" <>
           StringJoin[KeyValueMap[Function[{pq, rs},
               Module[{re = Select[rs, #[[1, 5]] === "re" &], im = Select[rs, #[[1, 5]] === "im" &]},
                "      { const double ar = " <> rsum[k, re, "_r"] <> ";\n        const double ai = " <> rsum[k, re, "_i"] <> ";\n" <>
                "        const double br = " <> rsum[k, im, "_r"] <> ";\n        const double bi = " <> rsum[k, im, "_i"] <> ";\n" <>
                "        gr += " <> rname[pq] <> "*(ar - bi); gi += " <> rname[pq] <> "*(ai + br); }\n"]], byR]] <>
           "      Gr" <> ToString[k - 1] <> "[j] += gr; Gi" <> ToString[k - 1] <> "[j] += gi;\n    }\n"],
         {k, Length[recs]}]];
    ];
    allcode = cDecl;
  ,
    allcode = StringRiffle[#["code"] & /@ terms, " "];
  ];
  used[pat_] := Union[StringCases[allcode, pat]];
  fields = used[RegularExpression["\\b(ddh|dh|h)[12]_\\d+\\b"]];
  rpow = used[RegularExpression["\\b[rf]_\\d+\\b"]];
  cpow = used[RegularExpression["\\b(M|omega1|omega2|mu1l1|mu1l2|mu2l1|mu2l2|mu3l1|mu3l2)_\\d+\\b"]];
  powDecl[name_] := Module[{base = StringReplace[name, "_" ~~ __ ~~ EndOfString -> ""], n = ToExpression[StringReplace[name, __ ~~ "_" -> ""]]},
    base <> StringJoin[Table["*" <> base, {n - 1}]]];
  fieldDecl[name_] := Module[{p = StringSplit[name, "_"]}, "const complex<double> " <> name <> " = " <> p[[1]] <> "[" <> p[[2]] <> "][j];"];

  src = StringJoin[
    "/* Generated by scripts/generate_kernels.wl from the PerturbationEquations package\n",
    " * (SchwarzschildQuadraticOperator[\"d2R\", \"Lorenz\", \"BLS\", \"BLS\"]), component i = ", ToString[i], ". */\n",
    "#include <complex>\n#include <cmath>\n#include <vector>\n#include \"Coupling.h\"\n#include \"R2G.h\"\n\nusing namespace std;\n\n",
    "void R2G_", ToString[i], "(const double M, const double r0,\n",
    "    const vector<double> &r_grid, const vector<double> &f_grid, const vector<double> &fp_grid,\n",
    "    const ModeBlock &h1, const ModeBlock &dh1, const ModeBlock &ddh1,\n",
    "    const ModeBlock &h2, const ModeBlock &dh2, const ModeBlock &ddh2,\n",
    "    int l1, int m1, int l2, int m2, double *const *Gr, double *const *Gi, size_t jb, size_t je)\n{\n",
    "  const complex<double> I(0.0, 1.0);\n",
    "  const double Omega = sqrt(M/(r0*r0*r0));\n",
    "  const double omega1 = m1*Omega, omega2 = m2*Omega;\n",
    "  const double mu1l1 = sqrt(l1*(l1+1.0)), mu1l2 = sqrt(l2*(l2+1.0));\n",
    "  const double mu2l1 = sqrt((l1-1.0)*(l1+2.0)), mu2l2 = sqrt((l2-1.0)*(l2+2.0));\n",
    "  const double mu3l1 = sqrt((l1-2.0)*(l1+3.0)), mu3l2 = sqrt((l2-2.0)*(l2+3.0));\n",
    StringJoin[("  const double " <> # <> " = " <> powDecl[#] <> ";\n") & /@ cpow],
    If[form === "factored", cDecl, ""],
    If[form === "real",
      StringJoin[Table["  const bool use" <> ToString[k - 1] <> " = (l1 >= " <> ToString[terms[[k, "req1"]]] <> " && l2 >= " <> ToString[terms[[k, "req2"]]] <> ");  /* spins (s3, s2, s1) = " <> ToString[terms[[k, "s"]]] <> " */\n", {k, Length[terms]}]] <>
      cDecl <>
      StringJoin[Module[{b = StringSplit[fname[#], "_"]}, "  const double *__restrict " <> fname[#] <> "re = " <> b[[1]] <> ".re[" <> b[[2]] <> "], *__restrict " <> fname[#] <> "im = " <> b[[1]] <> ".im[" <> b[[2]] <> "];\n"] & /@ fieldsF] <>
      StringJoin[Table["  double *__restrict Gr" <> ToString[k - 1] <> " = Gr[" <> ToString[k - 1] <> "], *__restrict Gi" <> ToString[k - 1] <> " = Gi[" <> ToString[k - 1] <> "];\n", {k, Length[terms]}]] <>
      "  (void) fp_grid; (void) Omega;\n#pragma omp simd\n  for (size_t j = jb; j < je; ++j) {\n" <>
      "    const double r = r_grid[j];\n    const double f = f_grid[j];\n    const double rinv = 1.0/r, finv = 1.0/f;\n    (void) rinv; (void) finv;\n" <>
      StringJoin[rdecl /@ rmons] <>
      StringJoin[("    const double " <> fname[#] <> "_r = " <> fname[#] <> "re[j], " <> fname[#] <> "_i = " <> fname[#] <> "im[j];\n") & /@ fieldsF] <>
      StringJoin[("    const double " <> pname @@ # <> "_r = " <> fname[#[[1]]] <> "_r*" <> fname[#[[2]]] <> "_r - " <> fname[#[[1]]] <> "_i*" <> fname[#[[2]]] <> "_i, " <>
          pname @@ # <> "_i = " <> fname[#[[1]]] <> "_r*" <> fname[#[[2]]] <> "_i + " <> fname[#[[1]]] <> "_i*" <> fname[#[[2]]] <> "_r;\n") & /@ prods] <>
      body <> "  }\n}\n\n",
    StringJoin[
    StringJoin[Table["  const bool use" <> ToString[k - 1] <> " = (l1 >= " <> ToString[terms[[k, "req1"]]] <> " && l2 >= " <> ToString[terms[[k, "req2"]]] <> ");  /* spins (s3, s2, s1) = " <> ToString[terms[[k, "s"]]] <> " */\n", {k, Length[terms]}]],
    "  (void) fp_grid; (void) Omega;\n",
    "  for (size_t j = jb; j < je; ++j) {\n",
    "    const double r = r_grid[j];\n    const double f = f_grid[j];\n",
    If[form === "factored",
      StringJoin[
        "    const double rinv = 1.0/r, finv = 1.0/f;\n    (void) rinv; (void) finv;\n",
        StringJoin[rdecl /@ rmons],
        StringJoin[("    const complex<double> " <> fname[#] <> " = " <> StringReplace[fname[#], "_" -> "["] <> "][j];\n") & /@ fieldsF],
        StringJoin[("    const complex<double> " <> pname @@ # <> " = " <> fname[#[[1]]] <> "*" <> fname[#[[2]]] <> ";\n") & /@ prods],
        body],
      StringJoin[
        StringJoin[("    const double " <> # <> " = " <> powDecl[#] <> ";\n") & /@ rpow],
        StringJoin[("    " <> fieldDecl[#] <> "\n") & /@ fields],
        StringJoin[Table["    if (use" <> ToString[k - 1] <> ") G[" <> ToString[k - 1] <> "][j] += " <> terms[[k, "code"]] <> ";\n", {k, Length[terms]}]]]],
    "  }\n}\n\n"]],
    "void R2C_", ToString[i], "(int l3, int m3, int l1, int m1, int l2, int m2, double *c)\n{\n",
    StringJoin[Table["  c[" <> ToString[k - 1] <> "] = " <> scale <> "*Coupling(l3,m3," <> ToString[terms[[k, "s"]][[1]]] <> ",l2,m2," <> ToString[terms[[k, "s"]][[2]]] <> ",l1,m1," <> ToString[terms[[k, "s"]][[3]]] <> ");\n", {k, Length[terms]}]],
    "}\n"];
  Export[FileNameJoin[{root, "src", "R2G_" <> ToString[i] <> ".cc"}], src, "Text"];
  Print["component ", i, " (", form, "): ", Length[terms], " couplings, ", Length[fields], " field arrays, ",
    If[form === "factored", ToString[Length[consts]] <> " constants, " <> ToString[Length[prods]] <> " products, " <> ToString[Length[rmons]] <> " r/f monomials, ", ""], "requirements ",
    ({#["req1"], #["req2"]} & /@ terms)],
  {i, 1, 10}];

header = StringJoin[
  "/* Generated by scripts/generate_kernels.wl: the d2R components split by coupling coefficient,\n",
  " *   d2R_i = sum_k R2C_i(l3, m3, ...)[k] * G_ik,  G_ik from R2G_i (independent of l3, m3),
 *   returned as real (Gr) and imaginary (Gi) planes. */\n",
  "#pragma once\n#include <complex>\n#include <cstddef>\n#include <vector>\n#include \"modeblock.h\"\n\n",
  StringJoin[Table["void R2G_" <> ToString[i] <> "(const double M, const double r0,\n" <>
     "  const std::vector<double> &r_grid, const std::vector<double> &f_grid, const std::vector<double> &fp_grid,\n" <>
     "  const ModeBlock &h1, const ModeBlock &dh1, const ModeBlock &ddh1,\n" <>
     "  const ModeBlock &h2, const ModeBlock &dh2, const ModeBlock &ddh2,\n" <>
     "  int l1, int m1, int l2, int m2, double *const *Gr, double *const *Gi, std::size_t jb, std::size_t je);\n" <>
     "void R2C_" <> ToString[i] <> "(int l3, int m3, int l1, int m1, int l2, int m2, double *c);\n", {i, 1, 10}]],
  "\n/* number of coupling coefficients of each component (index 0 unused) */\n",
  "static const int R2_NC[11] = {0, ", StringRiffle[ToString /@ Table[nc[i], {i, 1, 10}], ", "], "};\n",
  "#define R2_NC_MAX ", ToString[Max[Values[nc]]], "\n"];
Export[FileNameJoin[{root, "include", "R2G.h"}], header, "Text"];
Print["done"];
