// Run in an isolated VMH checkout, with the field and rank in BasicData.m.
AttachSpec("Magma/VMH-spec");
SetSeed(1);
V := VoronoiAlgorithm();
G, phi := Presentation(V);
matrices := [[[[Sprintf("\"%o\"", c) : c in Eltseq(phi(G.k)[i,j])]
    : j in [1..NumberOfColumns(phi(G.k))]]
    : i in [1..NumberOfRows(phi(G.k))]] : k in [1..Ngens(G)]];
forms := [[[[Sprintf("\"%o\"", c) : c in Eltseq(F`Matrix[i,j])]
    : j in [1..NumberOfColumns(F`Matrix)]]
    : i in [1..NumberOfRows(F`Matrix)]] : F in V`PerfectList];
relations := [Eltseq(LHS(r)*RHS(r)^-1) : r in Relations(G)];
major, minor, patch := GetVersion();
Write("reference.json", Sprintf(
    "{\"magma_version\":%o,\"perfect_forms\":%o,\"forms\":%o,\"matrices\":%o,\"relators\":%o,\"abelianization\":%o}",
    [major,minor,patch], #V`PerfectList, forms, matrices, relations,
    Invariants(AbelianQuotient(G))));
Write("complete.txt", "complete");
quit;
