// Input contains two finite presentations and words defining maps both ways.
load "presentation_maps.m";
FN := FreeGroup(#NativeToReference);
FR := FreeGroup(#ReferenceToNative);
Nfp := quo<FN | [FN!w : w in NativeRelations]>;
Rfp := quo<FR | [FR!w : w in ReferenceRelations]>;
function EvaluateWord(w, images, group)
    g := group!1;
    for a in w do
        g *:= images[Abs(a)]^Sign(a);
    end for;
    return g;
end function;
function RewriteWord(w,Q)
    // Magma's sequence-to-word coercion has a separate input-length limit.
    // Reduce chunks before multiplying them; this uses the same derived rules.
    value := Q!1;
    for first in [1..#w by 1024] do
        value *:= Q!w[first..Min(first+1023,#w)];
    end for;
    return value;
end function;
nToR := [FR!w : w in NativeToReference];
rToN := [FN!w : w in ReferenceToNative];
inN := [Eltseq(EvaluateWord(w,rToN,FN)) : w in ReferenceRelations] cat
    [Eltseq(EvaluateWord(NativeToReference[i],rToN,FN)*FN.i^-1) : i in [1..Ngens(Nfp)]];
inR := [Eltseq(EvaluateWord(w,nToR,FR)) : w in NativeRelations] cat
    [Eltseq(EvaluateWord(ReferenceToNative[i],nToR,FR)*FR.i^-1) : i in [1..Ngens(Rfp)]];

procedure ProveIdentities(G, words, label)
    // GrpFPElt equality is only free reduction. Every positive RWS reduction
    // is a consequence of the original relators, even without confluence.
    // An unreduced word is inconclusive: enlarge the completion search and
    // retry, rather than reporting a false mathematical counterexample.
    capacity := 32767; n := Ngens(G);
    relations := [Eltseq(LHS(r)*RHS(r)^-1) : r in Relations(G)];
    inputLength := &+[#r : r in relations];
    wordLength := Max([1] cat [#w : w in words]);
    orders := [];
    for k in [0..n-1] do
        order := [(j+k-1) mod n+1 : j in [1..n]];
        for p in [order,Reverse(order)] do
            if not p in orders then Append(~orders,p); end if;
        end for;
    end for;
    proved := [false : w in words];
    repeat
        for order in orders do
            // Relabel BOTH the defining relators and the queried words.
            // Positive proofs from different orderings are all consequences
            // of the original presentation and may therefore be combined.
            F := FreeGroup(n);
            H := quo<F | [F![Sign(i)*order[Abs(i)] : i in r] : r in relations]>;
            printf "%o: relation rewriting capacity %o, order %o\n",label,capacity,order;
            Q := RWSGroup(H : MaxRelations := capacity,
                MaxStates := Max(16*capacity,2*inputLength),
                MaxReduceLen := Max(32767,2*wordLength));
            for j in [1..#words] do
                if not proved[j] then
                    w := [Sign(i)*order[Abs(i)] : i in words[j]];
                    proved[j] := RewriteWord(w,Q) eq Q!1;
                end if;
            end for;
            printf "%o: proved %o/%o required identities\n",label,#[b : b in proved | b],#proved;
            if &and proved then break; end if;
        end for;
        capacity *:= 2;
    until &and proved;
end procedure;

print "Checking both homomorphisms and inverse maps on every abstract generator";
ProveIdentities(Nfp,inN,"Native presentation");
ProveIdentities(Rfp,inR,"Reference presentation");
Write("maps_complete.txt", "complete");
print "Abstract presentation isomorphism verified";
quit;
