# Compatibility loader for HAP releases predating the optional backend patch.
# It only changes calls carrying an explicit backend option. Patched HAP needs
# no wrapper, and re-reading this file leaves the registered providers intact.
# Some releases do not load the legacy entry points by default. Reading these
# function definitions does not execute Magma or access a VMH checkout.
if not IsBound(VoronoiComplexGL) or not IsBound(VoronoiComplexSL) then
    ReadPackage("hap", "lib/Voronoi/seb.gi");
fi;
if not IsBound(RegisterVoronoiBackend) then
BindGlobal("HAP_VORONOI_BACKENDS", rec());
BindGlobal("RegisterVoronoiBackend", function(name, provider)
    if not IsString(name) or Length(name)=0 or name="magma"
       or not IsFunction(provider) then
        Error("expected a nonempty backend name other than magma and a function");
    fi;
    HAP_VORONOI_BACKENDS.(name) := provider;
end);
BindGlobal("CallVoronoiBackend", function(name, kind, n, d)
    if not IsString(name) or not IsBound(HAP_VORONOI_BACKENDS.(name)) then
        Error("unknown Voronoi backend; register it with RegisterVoronoiBackend");
    fi;
    return HAP_VORONOI_BACKENDS.(name)(kind,n,d);
end);
BindGlobal("AV_HAP_MAGMA_GL", VoronoiComplexGL);
BindGlobal("AV_HAP_MAGMA_SL", VoronoiComplexSL);
VoronoiComplexGL := function(n,d)
    local backend;
    backend := ValueOption("backend");
    if backend<>fail and backend<>"magma" then
        return CallVoronoiBackend(backend,"GL",n,d);
    fi;
    return AV_HAP_MAGMA_GL(n,d);
end;
VoronoiComplexSL := function(n,d)
    local backend;
    backend := ValueOption("backend");
    if backend<>fail and backend<>"magma" then
        return CallVoronoiBackend(backend,"SL",n,d);
    fi;
    return AV_HAP_MAGMA_SL(n,d);
end;
fi;
