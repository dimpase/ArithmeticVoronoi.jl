// Run with Magma from a directory containing Page's unmodified
// KleinianGroups-1.0 archive. See docs/imaginary_quadratic.md.
// Override the JSON destination with magma -b OutputFile:=/path/results.json ...
AttachSpec("KleinianGroups-1.0/klngpspec");
import "KleinianGroups-1.0/bianchi.m": QuatToMatrix;
SetVerbose("Kleinian", 0);
SetColumns(0); // Magma's display continuations are not valid JSON escapes.
if not assigned OutputFile then
    OutputFile := "page_presentations.json";
end if;

function JSONArray(items)
    result := "[";
    for i := 1 to #items do
        if i gt 1 then result cat:= ","; end if;
        result cat:= items[i];
    end for;
    return result cat "]";
end function;

function MatrixJSON(M, D)
    rows := [];
    for i := 1 to 2 do
        entries := [];
        for j := 1 to 2 do
            a, b := Explode(Eltseq(M[i,j]));
            // Page uses alpha = sqrt(-d); our bundle uses w = (s+sqrt(D))/2.
            if D mod 4 eq 1 then
                a -:= b;
                b *:= 2;
            end if;
            Append(~entries, Sprintf("[\"%o\",\"%o\"]", a, b));
        end for;
        Append(~rows, JSONArray(entries));
    end for;
    return JSONArray(rows);
end function;

function GroupJSON(D, name, G, generators, seconds)
    matrices := JSONArray([MatrixJSON(QuatToMatrix(g), D) : g in generators]);
    relators := [Eltseq(LHS(r)*RHS(r)^-1) : r in Relations(G)];
    return Sprintf("{\"discriminant\":%o,\"group\":\"%o\",\"cpu_seconds\":%o,\"abelian_invariants\":%o,\"generator_matrices\":%o,\"relators\":%o}",
        D, name, seconds, Invariants(AbelianQuotient(G)), matrices, relators);
end function;

records := [];

for field in [<-3,3>, <-4,1>, <-7,7>, <-15,15>, <-23,23>] do
    D, d := Explode(field);
    // Each computation starts independently: do not inherit random state
    // or cached quaternion-algebra data from a previous group computation.
    SetSeed(1);
    O := BianchiOrder(d); // parameter d means Q(sqrt(-d)), not discriminant
    t := Cputime();
    _, faces, edges := NormalizedBasis(O);
    PSL, pgens := Presentation(faces, edges, O);
    psl_seconds := Cputime(t);
    SL, sgens := LiftPresentation(PSL, pgens, O, true);
    sl_seconds := Cputime(t);
    print "D, PSL/SL seconds:", D, sl_seconds;
    print "PSL abelian quotient:", AbelianQuotient(PSL);
    print "SL abelian quotient:", AbelianQuotient(SL);
    print "SL presentation:", SL;
    print "SL generator matrices:", [QuatToMatrix(g): g in sgens];
    Append(~records, GroupJSON(D, "PSL", PSL, pgens, psl_seconds));
    Append(~records, GroupJSON(D, "SL", SL, sgens, sl_seconds));

    SetSeed(1);
    O := BianchiOrder(d);
    t := Cputime();
    _, faces, edges := NormalizedBasis(O: GroupType := "Units");
    PGL, pgens := Presentation(faces, edges, O);
    GL, ggens := LiftPresentation(PGL, pgens, O, false);
    gl_seconds := Cputime(t);
    print "D, PGL/GL seconds:", D, gl_seconds;
    print "GL abelian quotient:", AbelianQuotient(GL);
    print "GL presentation:", GL;
    print "GL generator matrices:", [QuatToMatrix(g): g in ggens];
    Append(~records, GroupJSON(D, "GL", GL, ggens, gl_seconds));
end for;
major, minor, patch := GetVersion();
result := Sprintf("{\"schema\":\"av_page_presentations_v1\",\"magma_version\":\"%o.%o.%o\",\"kleinian_groups_version\":\"1.0\",\"random_seed_per_computation\":1,\"timing_note\":\"CPU seconds for domain and presentation, including the lift for SL/GL; excludes order construction, export and abelianization. SL includes the PSL computation.\",\"groups\":%o}",
    major, minor, patch, JSONArray(records));
PrintFile(OutputFile, result : Overwrite := true);
print "Saved exact presentations:", OutputFile;
quit;
