# The input uses integral column-action matrices and signed, oriented boundaries.
# All data are captured per object: constructing another complex cannot change it.
if not IsBound(HAPArithmeticVoronoiComplex) then
BindGlobal("HAPArithmeticVoronoiComplex", function(data, length)
    local dims, elts, stabs, evens, boundaries, group, dimension, stabilizer,
          boundary, action;
    dims := ShallowCopy(data.dimensions);
    elts := StructuralCopy(data.elts);
    stabs := List(data.stabilizers, layer -> List(layer, Group));
    evens := List(data.evenStabilizers, layer -> List(layer, Group));
    boundaries := StructuralCopy(data.boundaries);
    group := Group(data.generators);
    SetName(group, data.name);
    dimension := function(n)
        if n < 0 or n >= Length(dims) then return 0; fi;
        return dims[n+1];
    end;
    stabilizer := function(n,k) return stabs[n+1][AbsInt(k)]; end;
    boundary := function(n,k)
        if n=0 then return []; fi;
        if k>0 then return StructuralCopy(boundaries[n+1][k]); fi;
        return List(boundaries[n+1][-k], x -> [-x[1],x[2]]);
    end;
    action := function(n,k,g)
        local H, r, u;
        H := stabilizer(n,k);
        r := CanonicalRightCosetElement(H,One(H))^-1
             * CanonicalRightCosetElement(H,elts[g]^-1);
        u := r*elts[g];
        if u in evens[n+1][AbsInt(k)] then return 1; fi;
        return -1;
    end;
    return Objectify(HapNonFreeResolution,rec(
        dimension:=dimension, elts:=elts, group:=group, stabilizer:=stabilizer,
        action:=action, boundary:=boundary, homotopy:=fail, hasse:=fail,
        arithmeticVoronoi:=rec(discriminant:=data.discriminant,
                              groupKind:=data.groupKind, backend:="julia"),
        properties:=[["type","resolution"],["reduced",true],
                    ["length",length],["characteristic",0]]));
end);
fi;
