using ArithmeticVoronoi, Oscar, Test
const IAV=ArithmeticVoronoi

@testset "imaginary quadratic Hermitian model" begin
    for (D,m) in ((-3,6),(-4,4),(-7,2),(-15,2),(-23,2))
        P=arithmetic_problem(D);K=P.field;w=P.generator
        @test length(P.form_basis)==4
        @test discriminant(P.order)==D
        F=matrix(K,[3 w;P.root_sum-w 4P.root_product])
        A=trace_gram(P,F)
        @test positive_definite(A)
        @test form_from_trace(P,A)==F
        g=matrix(K,[1 w;0 1]);M=encode_matrix(P,g)
        @test trace_gram(P,IAV.act_form(P,F,g))==transpose(M)*A*M
        @test decode_matrix(P,M)==g
        @test_throws ArgumentError trace_gram(P,transpose(g)*F*g)
        # Conjugation preserves the identity Hermitian trace form but is not
        # an allowed K-linear lattice transformation.
        c=matrix(QQ,[1 P.root_sum;0 -1]);C=block_diagonal_matrix([c,c])
        B=trace_gram(P,identity_matrix(K,2))
        @test transpose(C)*B*C==B
        @test_throws GroupConstraintError decode_matrix(P,C)
        units=IAV.verified_units(P)
        @test units.rank==0 && units.torsion_order==m
        for k in 0:m-1
            @test IAV.finite_unit_exponent(units,units.torsion_generator^k)==k
        end
        H=IAV.act_form(P,F,g)
        h=equivalence(P,F,H)
        @test !isnothing(h) && IAV.act_form(P,F,h)==H
    end
end

@testset "imaginary presentations and published polytope types" begin
    # Yasaki, arXiv:0908.1762v2, Tables 2--4. Tuples are (vertices,edges,faces).
    expected=Dict(-3=>[(4,6,4)],-4=>[(6,12,8)],-7=>[(6,9,5)],
        -15=>[(4,6,4),(6,12,8)],-23=>[(5,8,5),(6,9,5),(6,12,8)])
    for D in (-3,-4,-7,-15,-23)
        mktempdir() do directory
            R=withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/backend/config.json") do
                run_field_pipeline(D;backend=:native,directory,reuse_baseline=false)
            end
            P=R.enumeration.problem;K=P.field;I=identity_matrix(K,2)
            shapes=sort!([(length(unique(IAV.primitive_ray(v) for v in IAV.matrix_rows(r.evaluation_rays))),
                length(R.context.ridges[i]),length(r.facets)) for (i,r) in enumerate(R.enumeration.representatives)])
            @test shapes==expected[D]
            @test verify_presentation(R.gl) && verify_sl_presentation(R.sl)
            @test length(R.sl.transversal_words)==R.sl.units.torsion_order
            @test isempty(R.sl.scalar_word) && isnothing(R.sl.scalar_reduction)
            @test all(iszero,R.sl.scalar_corrections)
            @test IAV.evaluate_word(R.sl.generator_matrices,R.sl.minus_identity_word)==-I
            e=R.sl.units.torsion_generator
            matrices=[-I,matrix(K,[1 P.generator;0 1]),matrix(K,[e 0;0 inv(e)])]
            for g in matrices,target in (:GL,:SL,:PSL),simplified in (false,true)
                c=matrix_to_word(R,g;target,simplified)
                @test verify_word_certificate(R,g,c)
            end
            # Non-self-inverse determinant images exercise signed cyclic
            # Schreier transitions (not the xor quotient of the real case).
            u=matrix(K,[e 0;0 1]);v=matrix(K,[1 1;0 1])
            g=u*v*inv(u)*inv(v)
            @test det(g)==1
            c=matrix_to_word(R,g;target=:SL)
            @test verify_word_certificate(R,g,c)
        end
    end
end
