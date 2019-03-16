(* ::Package:: *)

BeginPackage["FastRecurrenceFilters`"]

ClearAll[FastRecurrenceFilter, Biquad, BiquadReverse, BiquadUp2, BiquadDown2];
FastRecurrenceFilter::usage="\!\(\*RowBox[{\"FastRecurrenceFilter\"}]\) is a \
compiled implementation of some numerical recurrence filters."
Biquad::usage="Biquad implementation."
BiquadReverse::usage="Biquad implementation with reversed output."
BiquadUp2::usage="2x upsampling biquad implementation."
BiquadDown2::usage="2x downsampling biquad implementation."

Begin["`Private`"]

ClearAll[compilationTarget, runtimeOptions];
compilationTarget = "C";
runtimeOptions = {"Speed",
  "CatchMachineOverflow" -> False,
  "CatchMachineIntegerOverflow" -> False,
  "CompareWithTolerance" -> False,
  "EvaluateSymbolically" -> False};

ClearAll[biquad2t];
With[{ro = runtimeOptions}, biquad2t = 
  Compile[{{x, _Real, 
     1}, {b0, _Real}, {b1, _Real}, {b2, _Real}, {a0, _Real}, {a1, \
_Real}, {a2, _Real}}, 
   With[{b0n = b0/a0, b1n = b1/a0, b2n = b2/a0, a1n = a1/a0, 
     a2n = a2/a0}, 
    Module[{z1 = 0., z2 = 0.}, 
     Table[With[{y = z1 + x[[i]] b0n}, z1 = z2 + x[[i]] b1n - y a1n; 
       z2 = x[[i]] b2n - y a2n; y], {i, Length[x]}]]], 
   CompilationTarget -> compilationTarget, RuntimeOptions -> ro];

ClearAll[biquad2tUp2];
With[{ro = runtimeOptions}, biquad2tUp2 = 
  Compile[{{x0, _Real, 
     1}, {b0, _Real}, {b1, _Real}, {b2, _Real}, {a0, _Real}, {a1, \
_Real}, {a2, _Real}}, 
   With[{a = b0/a0, a1n = -a1/a0, a2n = -a2/a0}, 
    Module[{z0 = 0., z1 = 0.}, 
     a Table[With[{xi = x0[[i]]}, 
        With[{yi0 = xi + z0}, 
         With[{yi1 = xi + xi + yi0 a1n + z1 a2n}, 
          z0 = xi + yi0 a2n + yi1 a1n; z1 = yi1; 
          If[j == 1, yi0, yi1]]]], {j, 2}, {i, Length[x0]}]]], 
   CompilationTarget -> compilationTarget, RuntimeOptions -> ro]];

ClearAll[biquad2tDown2];
With[{ro = runtimeOptions}, biquad2Down2 = 
  Compile[{{x0, _Real, 1}, {x1, _Real, 
     1}, {b0, _Real}, {b1, _Real}, {b2, _Real}, {a0, _Real}, {a1, \
_Real}, {a2, _Real}}, 
   With[{a = b0/a0, a1n = a1/a0, a2n = a2/a0}, 
    Module[{z0 = 0., z1 = 0.}, 
     Table[With[{x0i = x0[[i]], x1i = x1[[i]]}, 
       With[{z0i = x0i - (z1 a1n + z0 a2n)}, 
        With[{z1i = x1i - (z0i a1n + z1 a2n)}, 
         With[{yi = z0i + z0i + z1i + z1}, z0 = z0i; z1 = z1i; 
          yi a]]]], {i, Length[x0]}]]], CompilationTarget -> compilationTarget,
   RuntimeOptions -> ro]];

ClearAll[biquad2tr];
With[{ro = runtimeOptions}, biquad2tr =
  Compile[{{x, _Real,
     1}, {b0, _Real}, {b1, _Real}, {b2, _Real}, {a0, _Real}, {a1, \
_Real}, {a2, _Real}},
   With[{b0n = b0/a0, b1n = b1/a0, b2n = b2/a0, a1n = a1/a0,
     a2n = a2/a0},
    Module[{z1 = 0., z2 = 0.},
     Table[With[{y = z1 + x[[i]] b0n}, z1 = z2 + x[[i]] b1n - y a1n;
       z2 = x[[i]] b2n - y a2n; y], {i, Length[x], 1, -1}]]],
   CompilationTarget -> compilationTarget, RuntimeOptions -> ro]];

FastRecurrenceFilter::noimpl =
  "No matching implementation found for coefficients `1` and ratio `2`, falling\
 back to \!\(\*RowBox[{\"RecurrenceFilter\"}]\).";

ClearAll[fastRecurrenceFilter];
fastRecurrenceFilter[c_, x_, r_:1] := With[{
  a = First[c],
  b = Last[c]},
  Message[FastRecurrenceFilter::noimpl, c, r];
  Switch[r,
    1, RecurrenceFilter[{a, b}, x],
    (*I need to revise the following two cases if I want to use this*)
    ratio_/;ratio > 1, Transpose@Partition[RecurrenceFilter[{a, b},
      Upsample[x, r, 1, 0.]], r],
    ratio_/;ratio < 1, 
      Downsample[RecurrenceFilter[{a, b}, Catenate[Transpose@x]], 1/r, 1/r]]
  ];
(*fastRecurrenceFilter[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, x_, 1] :=
  biquad2t[x, b0, b1, b2, a0, a1, a2];*)
Biquad[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, x_] :=
  biquad2t[x, b0, b1, b2, a0, a1, a2];

BiquadReverse[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, x_] :=
  biquad2tr[x, b0, b1, b2, a0, a1, a2];

(*For some reason, biquad2tUp2 and the fallback in fastRecurrenceFilter do not
  yield the same output*)
(*fastRecurrenceFilter[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, x_, 2] :=
  biquad2tUp2[x, b0, b1, b2, a0, a1, a2];*)
BiquadUp2[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, x_] :=
  biquad2tUp2[x, b0, b1, b2, a0, a1, a2];

(*For some reason, biquad2tDown2 and the fallback in fastRecurrenceFilter do not
  yield the same output*)
(*fastRecurrenceFilter[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, x_, 1/2] :=
  biquad2Down2[First[x], Last[x], b0, b1, b2, a0, a1, a2];*)
BiquadDown2[{{a0_, a1_, a2_}, {b0_, b1_, b2_}}, {x0_, x1_}] :=
  biquad2Down2[x0, x1, b0, b1, b2, a0, a1, a2];

FastRecurrenceFilter = fastRecurrenceFilter;

End[]

EndPackage[]
