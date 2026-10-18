using ArithmeticVoronoi, Oscar, Test, JSON3
import Hecke
const AV=ArithmeticVoronoi
@testset "maximal-order trace model and exact K-linear action" begin
    for D in (5,8)
        P=arithmetic_problem(D);K,w=P.field,P.generator
        F=matrix(K,[3+w//4 1//3+w//7;1//3+w//7 4-w//5])
        A=trace_gram(P,F)
        @test positive_definite(A)
        @test form_from_trace(P,A)==F
        @test A*P.multiplication==transpose(P.multiplication)*A
        @test rank(matrix(QQ,[B[i,j] for B in P.form_basis,(i,j) in [(i,j) for i in 1:4 for j in 1:4]]))==6
        for z in ([1,0,0,0],[0,1,0,0],[1,-2,3,1],[-4,3,2,-5])
            v=field_vector(P,z)
            @test coordinate_vector(P,v)==QQ.(z)
            @test AV.quadratic_value(A,z)==tr(sum(v[i]*F[i,j]*v[j] for i in 1:2,j in 1:2))
        end
        eps=D==5 ? w : 1+w
        matrices=[matrix(K,[1 1;0 1]),matrix(K,[1 w;0 1]),matrix(K,[0 -1;1 0]),matrix(K,[eps 0;0 1]),-identity_matrix(K,2)]
        for g in matrices
            R=encode_matrix(P,g)
            @test check_group_matrix(P,R)
            @test decode_matrix(P,R)==g
            @test trace_gram(P,transpose(g)*F*g)==transpose(R)*A*R
            bg=backend_gram(P,transpose(R))
            @test internal_matrix(P,bg)==g
            @test internal_gram(P,bg*backend_gram(P,A)*transpose(bg))==transpose(R)*A*R
            @test equivalence(P,F,transpose(g)*F*g)!==nothing
        end
        C=block_diagonal_matrix([matrix(QQ,[1 P.root_sum;0 -1]),matrix(QQ,[1 P.root_sum;0 -1])])
        AI=trace_gram(P,identity_matrix(K,2))
        @test transpose(C)*AI*C==AI  # a genuine integral isometry which is not K-linear
        @test_throws GroupConstraintError decode_matrix(P,C)
        @test_throws GroupConstraintError decode_matrix(P,2identity_matrix(QQ,4))
        if D==8
            @test_throws GroupConstraintError decode_matrix(P,transpose(encode_matrix(P,matrices[2])))
        end
        bad=deepcopy(A);bad[1,1]+=1
        @test_throws ArgumentError form_from_trace(P,bad)
        stab=finite_stabilizer(P,identity_matrix(K,2))
        @test stab.order>0
        @test -identity_matrix(K,2) in stab.generators
    end
end
@testset "Hecke short vectors use lattice coordinates and squared lengths" begin
    B=matrix(QQ,[2 1;0 3]);L=integer_lattice(B)
    @test Set(AV.signed_key(v) for v in shortest_vectors(L))==Set([(BigInt(1),BigInt(0))])
    @test minimum(L)==5
    @test all(AV.quadratic_value(gram_matrix(L),v)==minimum(L) for v in shortest_vectors(L))
    A=matrix(QQ,[2//3 1//3;1//3 2//3]);data=shortest_data(A)
    @test data.minimum==QQ(2,3)
    @test length(data.vectors)==3
    @test Set(AV.signed_key(v) for (v,n) in short_vectors(data.lattice,data.minimum,data.minimum))==Set(Tuple(v) for v in data.vectors)
end
println("Exact coordinate, constrained-isometry, stabilizer and short-vector model tests passed.")
