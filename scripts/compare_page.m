// Run with Magma from a directory containing Page's unmodified
// KleinianGroups-1.0 archive. See docs/imaginary_quadratic.md.
AttachSpec("KleinianGroups-1.0/klngpspec");
import "KleinianGroups-1.0/bianchi.m": QuatToMatrix;
SetVerbose("Kleinian", 0);
SetSeed(1);

for field in [<-3,3>, <-4,1>, <-7,7>, <-15,15>, <-23,23>] do
    D, d := Explode(field);
    O := BianchiOrder(d); // parameter d means Q(sqrt(-d)), not discriminant
    t := Cputime();
    _, faces, edges := NormalizedBasis(O);
    PSL, pgens := Presentation(faces, edges, O);
    SL, sgens := LiftPresentation(PSL, pgens, O, true);
    print "D, PSL/SL seconds:", D, Cputime(t);
    print "PSL abelian quotient:", AbelianQuotient(PSL);
    print "SL abelian quotient:", AbelianQuotient(SL);
    print "SL presentation:", SL;
    print "SL generator matrices:", [QuatToMatrix(g): g in sgens];

    t := Cputime();
    _, faces, edges := NormalizedBasis(O: GroupType := "Units");
    PGL, pgens := Presentation(faces, edges, O);
    GL, ggens := LiftPresentation(PGL, pgens, O, false);
    print "D, PGL/GL seconds:", D, Cputime(t);
    print "GL abelian quotient:", AbelianQuotient(GL);
    print "GL presentation:", GL;
    print "GL generator matrices:", [QuatToMatrix(g): g in ggens];
end for;
quit;
