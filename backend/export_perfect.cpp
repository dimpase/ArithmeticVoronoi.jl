// Narrow exporter linked against GPL-2.0 polyhedral_common (Mathieu Dutour Sikiric).
// Uses upstream initial forms, neighbours, constrained equivalences and stabilizers.
// It deliberately constructs no lower-dimensional complex or boundary matrices.
#include "NumberTheory.h"
#include "Permutation.h"
#include "Group.h"
#include "perfect_tspace.h"
#include <filesystem>
#include <fstream>
#include <cstdlib>

using Q=mpq_class;
using Z=mpz_class;
using Perm=permutalib::SingleSidedPerm<uint32_t>;
using Group=permutalib::Group<Perm,Z>;
using BaseData=DataPerfectTspaceFunc<Q,Z,Group>;

struct ExportFailure:std::runtime_error {
  std::string code;
  ExportFailure(std::string c,std::string message):std::runtime_error(message),code(std::move(c)){}
};
void quoted(std::ostream& os,std::string const& value) {
  char const* hex="0123456789abcdef";os << '"';
  for(unsigned char c:value) {
    if(c=='"' || c=='\\')os << '\\' << c;
    else if(c<32)os << "\\u00" << hex[c>>4] << hex[c&15];
    else os << c;
  }
  os << '"';
}
void report_error(std::string const& code,std::string const& message) {
  std::cerr << "AV_ERROR {\"schema\":\"av_error_v1\",\"code\":";
  quoted(std::cerr,code);std::cerr << ",\"message\":";quoted(std::cerr,message);std::cerr << "}\n";
}

template<class T> void matrix_json(std::ostream& os,MyMatrix<T> const& A) {
  os << '[';
  for(int i=0;i<A.rows();++i) {
    if(i)os << ',';
    os << '[';
    for(int j=0;j<A.cols();++j) {
      if(j)os << ',';
      os << '"' << A(i,j) << '"';
    }
    os << ']';
  }
  os << ']';
}
template<class T> void matrices_json(std::ostream& os,std::vector<MyMatrix<T>> const& list) {
  os << '[';
  for(size_t i=0;i<list.size();++i) {if(i)os << ',';matrix_json(os,list[i]);}
  os << ']';
}

struct LabelledAdjacency {
  Face incidence;
  MyMatrix<Z> transporter_row;
  MyMatrix<Q> neighbour_gram;
  MyMatrix<Z> neighbour_vectors;
  Q neighbour_minimum;
};
struct ExportData {
  BaseData base;
  std::string journal;
  size_t completed_objects=0;
  using Tobj=BaseData::Tobj;
  using TadjI=BaseData::TadjI;
  using TadjO=LabelledAdjacency;
  std::ostream& get_os(){return base.get_os();}
  Tobj f_init(){return base.f_init();}
  size_t f_hash(size_t const& seed,Tobj const& x){return base.f_hash(seed,x);}
  Tobj f_adji_obj(TadjI const& x){return base.f_adji_obj(x);}
  LabelledAdjacency augment(PerfectTspace_AdjO<Z> const& label,TadjI const& y) {
    return {label.eInc,label.eBigMat,y.Gram,conversion_and_duplication<Z,Z>(y.tsp.rec_shv.SHV),y.tsp.rec_shv.min};
  }
  std::optional<TadjO> f_repr(Tobj const& x,TadjI const& y) {
    auto result=base.f_repr(x,y);
    if(!result)return {};
    return augment(*result,y);
  }
  std::pair<Tobj,TadjO> f_spann(TadjI const& y) {
    auto result=base.f_spann(y);
    return {result.first,augment(result.second,y)};
  }
  std::optional<std::vector<TadjI>> f_adj(Tobj& x) {
    auto result=base.f_adj(x);
    if(result) {
      ++completed_objects;
      std::ofstream out(journal,std::ios::app);
      out.exceptions(std::ios::failbit|std::ios::badbit);
      out << "{\"schema\":\"av_progress_v1\",\"status\":\"partial\",\"completed_objects\":"
          << completed_objects << ",\"facet_orbits\":" << result->size() << ",\"gram\":";
      matrix_json(out,x.Gram);out << "}\n";
    }
    return result;
  }
};

int main(int argc,char** argv) {
  std::string phase="input";
  try {
    if(argc!=4)throw std::runtime_error("usage: AV_ExportPerfectOrbits input.nml output.json request-token");
    std::string token=argv[3];
    if(token.size()!=64 || token.find_first_not_of("0123456789abcdef")!=std::string::npos)
      throw std::runtime_error("request token must contain 64 lowercase hexadecimal digits");
    std::srand(1);srandom(1);set_random_seed(1);
    FullNamelist input=NAMELIST_GetStandard_ENUMERATE_PERFECT_COMPLEX_TSPACE();
    NAMELIST_ReadNamelistFile(argv[1],input);
    auto const& config=input.get_block("TSPACE");
    auto const& data=input.get_block("DATA");
    if(data.get_string("arithmetic")!="gmp")throw ExportFailure("unsupported_arithmetic","only gmp is compiled into this exporter");
    if(data.get_bool("ComputeBoundary") || data.get_bool("ComputePreImage") || data.get_bool("ComputeComplex"))
      throw ExportFailure("unsupported_operation","exporter enumerates perfect forms and facet orbits only");
    auto kind=config.get_string("TypeTspace");
    if(kind!="Classic" && kind!="RealQuad")throw std::runtime_error("unsupported form space");
    if(kind=="Classic" && config.get_int("ClassicDim")!=2)throw std::runtime_error("only classical rank two is initially supported");
    if(kind=="RealQuad") {
      int s=config.get_int("RealImagSum"),p=config.get_int("RealImagProd");
      if(config.get_int("RealImagDim")!=2 || !((s==1 && (p==-1 || p==-3))||(s==0 && (p==-2 || p==-3 || p==-10))))
        throw std::runtime_error("only maximal-order discriminants 5, 8, 12, 13 and 40, rank two, are supported");
      if(config.get_string("ListComm")!="Use_realimag")throw std::runtime_error("K-linear commuting constraint is mandatory");
    }
    auto space=ReadTspace<Q,Z,Group>(config,std::cerr);
    phase="enumeration";
    auto heuristics=AllStandardHeuristicSerial<Q,Z>(space.ListMat.size(),std::cerr);
    RecordDualDescOperation<Q,Group> dual(heuristics,std::cerr);
    DataPerfectTspace<Q,Z,Group> kernel{space,OnlineHierarchicalMatrixReduction<Z>(space.n,std::cerr),false,false,std::move(dual)};
    ExportData engine{BaseData{std::move(kernel)},std::string(argv[2])+".progress.jsonl"};
    auto invalid=[](BaseData::Tobj const&){return false;};
    // No internal time cap: upstream's serial routine otherwise returns a partial
    // database without a completion flag. The adapter enforces an outer deadline.
    auto result=EnumerateAndStore_Serial(engine,invalid,0);
    if(!result)throw std::runtime_error("enumeration rejected an object");
    if(engine.completed_objects!=result->size())throw std::runtime_error("unprocessed orbit in enumeration queue");
    phase="output";
    std::string temporary=std::string(argv[2])+".writing";
    std::ofstream out(temporary);
    out.exceptions(std::ios::failbit|std::ios::badbit);
    out << "{\"schema\":\"av_polyhedral_export_v1\",\"request_token\":\"" << token
        << "\",\"backend_commit\":\"cfcd79dbddf67a494a1c71b776174b190da5209a\","
        << "\"arithmetic\":\"gmp\",\"seed\":1,\"coordinates\":\"grouped rows (x1,x2,y1,y2); classical rows (x1,x2)\","
        << "\"action\":\"H=M*F*transpose(M); row vectors act by right multiplication\","
        << "\"computed_levels\":[\"perfect_forms\",\"facet_orbits\",\"labelled_neighbours\"],"
        << "\"boundary_square\":\"not_computed\",\"minimal_vector_convention\":\"both signs in conversion_and_duplication order\",\"form_basis\":";
    matrices_json(out,space.ListMat);out << ",\"commuting_matrices\":";matrices_json(out,space.ListComm);
    out << ",\"pointwise_stabilizer\":";matrices_json(out,space.PtStabGens);
    out << ",\"representatives\":[";
    for(size_t i=0;i<result->size();++i) {
      if(i)out << ',';
      auto const& entry=(*result)[i];auto const& obj=entry.x;
      out << "{\"id\":" << i+1 << ",\"gram\":";matrix_json(out,obj.Gram);
      out << ",\"minimum\":\"" << obj.tsp.rec_shv.min << "\",\"minimal_vectors\":";matrix_json(out,conversion_and_duplication<Z,Z>(obj.tsp.rec_shv.SHV));
      out << ",\"stabilizer_order\":\"" << obj.GRP.size() << "\",\"stabilizer_generators_row\":";matrices_json(out,obj.GRP_matr);
      out << ",\"neighbours\":[";
      for(size_t j=0;j<entry.ListAdj.size();++j) {
        if(j)out << ',';
        auto const& a=entry.ListAdj[j];
        out << "{\"target\":" << a.iOrb+1 << ",\"facet_vector_indices_zero_based\":[";
        bool first=true;
        for(size_t k=0;k<a.x.incidence.size();++k)if(a.x.incidence[k]){if(!first)out << ',';out << k;first=false;}
        out << "],\"transporter_row\":";matrix_json(out,a.x.transporter_row);
        out << ",\"neighbour_gram\":";matrix_json(out,a.x.neighbour_gram);
        out << ",\"neighbour_minimum\":\"" << a.x.neighbour_minimum << "\",\"neighbour_vectors\":";matrix_json(out,a.x.neighbour_vectors);
        out << '}';
      }
      out << "]}";
    }
    out << "],\"pending_facets\":[],\"completed\":true}\n";
    out.flush();out.close();
    std::filesystem::rename(temporary,argv[2]);
    std::cerr << "AV: completed " << result->size() << " perfect-form orbits; no lower complex computed\n";
    return 0;
  } catch(ExportFailure const& e) {
    report_error(e.code,e.what());return 64;
  } catch(TerminalException const& e) {
    report_error(phase=="input"?"malformed_input":"backend_failure","upstream TerminalException "+std::to_string(e.eVal));return e.eVal?e.eVal:1;
  } catch(std::exception const& e) {
    report_error(phase=="input"?"malformed_input":phase=="output"?"unwritable_output":"backend_failure",e.what());return 70;
  } catch(...) {
    report_error("backend_failure","unknown nonstandard upstream exception");return 70;
  }
}
