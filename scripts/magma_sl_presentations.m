// Independent VMH SL presentation, followed by its full scalar quotient.
// parameters.m is written by the Julia validation runner.
AttachSpec("Magma/VMH-spec");
load "parameters.m";
SetSeed(1);
function ScalarCenter(K, dimension)
    R<x> := PolynomialRing(K);
    roots := [r[1] : r in Roots(x^dimension-1)];
    d := #roots;
    return [a : a in roots | #{a^j : j in [0..dimension-1]} eq d][1], d;
end function;
if FieldDiscriminant lt 0 then
    V := VoronoiAlgorithm(: SL:=true);
    S, phi := Presentation(V: SL:=true);
    K := BaseRing(phi(S.1));
    z, centralOrder := ScalarCenter(K,Dimension);
    central := SolveWordProblem(GL(Dimension,K)!(z*IdentityMatrix(K,Dimension)),V);
elif FieldDiscriminant eq 1 then
    // Real VMH's presentation routine assumes -I belongs to its group.
    // Compute GL afresh, then let Magma rewrite its determinant kernel.
    V := VoronoiAlgorithm();
    G, glMatrices := Presentation(V);
    K := BaseRing(glMatrices(G.1));
    U := SymmetricGroup(2);
    determinant := hom<G -> U | [Determinant(glMatrices(G.i)) eq 1
        select U!1 else U!(1,2) : i in [1..Ngens(G)]]>;
    H := Kernel(determinant);
    assert Index(G,H) eq 2;
    z, centralOrder := ScalarCenter(K,Dimension);
    // Give the central word explicitly to the subgroup constructor. Direct
    // coercion of an arbitrary GL word into the coset-table kernel does not
    // provide the required change of alphabet.
    if centralOrder gt 1 then
        raw := [p[2] : p in V`StabilizerPreimages
            | p[1] eq z*IdentityMatrix(K,Dimension)][1];
        glCentral := V`GNonsimplifiedToSimplifiedHom(raw);
        assert glMatrices(glCentral) eq z*IdentityMatrix(K,Dimension);
        H := sub<G | GeneratingWords(G,H), glCentral>;
        assert Index(G,H) eq 2;
        hCentral := H.Ngens(H);
        assert glMatrices(G!hCentral) eq z*IdentityMatrix(K,Dimension);
    else
        hCentral := H!1;
    end if;
    S0, rewrite := Rewrite(G,H);
    S, toS := Simplify(S0);
    // Rewrite and Simplify return subgroups, so coercion gives the inclusion.
    // Their returned maps go from the old alphabet to the new one.
    phi := hom<S -> GL(Dimension,K) | [glMatrices(G!S.i) : i in [1..Ngens(S)]]>;
    central := toS(rewrite(hCentral));
else
    error "This comparison runner covers rational and imaginary quadratic fields";
end if;
assert &and[Determinant(phi(S.i)) eq 1 : i in [1..Ngens(S)]];
assert phi(central) eq z*IdentityMatrix(K,Dimension);
Q, quotient := quo<S | central>;
P := Simplify(Q);
assert Ngens(Q) eq Ngens(S);
projectiveMatrices := [phi(S!Eltseq(Q!P.i)) : i in [1..Ngens(P)]];

function Data(G, images)
    matrices := [[[[Sprintf("\"%o\"",c) : c in Eltseq(M[i,j])]
        : j in [1..NumberOfColumns(M)]] : i in [1..NumberOfRows(M)]] : M in images];
    relations := [Eltseq(LHS(r)*RHS(r)^-1) : r in Relations(G)];
    return Sprintf("{\"matrices\":%o,\"relators\":%o,\"abelianization\":%o}",
        matrices,relations,Invariants(AbelianQuotient(G)));
end function;
major,minor,patch := GetVersion();
Write("reference_sl.json", Sprintf("{\"magma_version\":%o,\"center_order\":%o,\"SL\":%o,\"PSL\":%o}",
    [major,minor,patch],centralOrder,
    Data(S,[phi(S.i) : i in [1..Ngens(S)]]),Data(P,projectiveMatrices)));
Write("sl_complete.txt","complete");
quit;
