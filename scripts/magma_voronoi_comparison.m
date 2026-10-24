// Run in an isolated VMH checkout with BasicData.m and AVComputeComplex set.
AttachSpec("Magma/VMH-spec");
SetSeed(1);
major, minor, patch := GetVersion();
Write("magma_version.json", Sprint([major,minor,patch]));
t := Cputime();
V := VoronoiAlgorithm();
forms := [[[[Sprintf("\"%o\"", c) : c in Eltseq(F`Matrix[i,j])]
    : j in [1..NumberOfColumns(F`Matrix)]]
    : i in [1..NumberOfRows(F`Matrix)]] : F in V`PerfectList];
Write("forms.json", Sprint(forms));
printf "Perfect forms: %o; enumeration CPU seconds: %o\n", #forms, Cputime(t);
if AVComputeComplex then
    G, phi := Presentation(V);
    Write("abelianization.json", Sprint(Invariants(AbelianQuotient(G))));
    ComputeComplexGL("complex_gl.g", V);
    K := BaseRing(V`PerfectList[1]`Matrix);
    if Degree(K) eq 2 and Discriminant(Integers(K)) lt 0 then
        D := Discriminant(Integers(K));
        m := D eq -4 select 4 else (D eq -3 select 6 else 2);
        z := D eq -4 select K.1 else (D eq -3 select (1+K.1)/2 else K!(-1));
        n := NumberOfRows(V`PerfectList[1]`Matrix);
        reps := [DiagonalMatrix([z^j] cat [K!1 : i in [2..n]]) : j in [0..m-1]];
        ComputeComplexLowIndexSubgroup("complex_sl.g", reps, IsInSL, V);
    end if;
end if;
printf "Total CPU seconds: %o\n", Cputime(t);
Write("complete.txt", "complete");
quit;
