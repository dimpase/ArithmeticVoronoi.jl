using ArithmeticVoronoi, JSON, SHA, Test

@testset "serialized certificate SHA-256 compatibility" begin
    @test ArithmeticVoronoi.text_sha("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    payload=JSON.json((;field="Q(√−3)",relators=repeat([[1,-2,3]],100_000)))
    for text in ("", "Q(√−3)", SubString("abcde",2,4), payload)
        # An independent streaming path checks the exact UTF-8 bytes, including
        # a large presentation payload, without Julia 1.10's CodeUnits issue.
        @test ArithmeticVoronoi.text_sha(text)==bytes2hex(sha256(IOBuffer(text)))
    end
end
