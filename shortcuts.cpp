#include "DGtal/base/Common.h"

#include "DGtal/helpers/Stddefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/shapes/MeshVoxelizer.h"
#include "DGtal/topology/DigitalSetBoundary.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <igl/writeDMAT.h>
#include <igl/readDMAT.h>
#include <igl/writeOBJ.h>

#include <igl/readSTL.h>

#include <Eigen/Core>

using namespace std;
using namespace DGtal;
// using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

typedef int Integer;                 // choose your digital line here.
typedef SpaceND<3,Integer> Z3;       // Z^3
typedef MetricAdjacency<Z3,1> Adj6;  // 6-adjacency type
typedef Z3::Point Point;
typedef HyperRectDomain< Z3 > Domain;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet;

typedef MetricAdjacency< Z3, 2 > Adj18;
typedef DigitalTopology< Adj6, Adj18 > DT6_18;
Adj6 adj6;
Adj18 adj18;
DT6_18 dt6_18( adj6, adj18, JORDAN_DT );
//typedef DigitalTopology< Adj6, Adj18 > DT6_18;
typedef Object<DT6_18, DigitalSet> ObjectType;

typedef DigitalSetBoundary< Z3i::KSpace, DigitalSet > Boundary;

namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("inputFilename,f", po::value<std::string>(), "input matrix of visibilities (.dmat)" )
  ("outputFilename,o", po::value<std::string>(), "output obj file (.obj)" )
  ("isovalue,i",  po::value<float>()->default_value(1.0), "isovalue (not normalized) of surface to extract in [1, num_views]" );

  bool parseOK=true;
  po::variables_map vm;
  try
  {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }
  catch(const std::exception& ex)
  {
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);

  if( !parseOK || vm.count("help")||argc<=1)
  {
    std::cout << "Usage: " << argv[0] << " -f [input file] -o [output]\n"
              << "Export the boundary of a volume to OBJ format."<< endl
              << general_opt << "\n";
    return 0;
  }

  std::string inputFilename;
  std::string outputFilename;

  if(! vm.count("inputFilename"))
  {
    trace.error() << " The input file name was not defined" << endl;
    return 0;
  }

  if(! vm.count("outputFilename"))
  {
    outputFilename = "test.obj";
    return 0;
  }

  inputFilename = vm["inputFilename"].as<std::string>();
  outputFilename = vm["outputFilename"].as<string>();
  float isovalue = vm["isovalue"].as<float>();

  auto params = SH3::defaultParameters()
            | SHG3::defaultParameters();

    // auto al_capone = SH3::makeBinaryImage( "../../DGtal_source/examples/samples/Al.100.vol", params );
    // auto K         = SH3::getKSpace( al_capone );
    // auto surface   = SH3::makeLightDigitalSurface( al_capone, K, params );
    // trace.info() << "#surfels=" << surface->size() << std::endl;

    // auto canonical_scene = SH3::makeBinaryImage( "../hidden-supports/gpu/build/output_vol_001.pgm3d", params );
  Eigen::MatrixXf S;
  igl::readDMAT(inputFilename, S);
  // int isovalue = 61;
  Eigen::MatrixXf S_iso = (S.array() >= isovalue).select(isovalue, S.array()-S.array());
  std::cout << S_iso.rows() << " " << S_iso.cols() << std::endl;

  Point p1( 1, 1, 1 );
  Point p2( 200,200,145 );
  Domain domain( p1, p2 );

  // canonical scene
  DigitalSet scene_set( domain );

  // MeshVoxelizer<DigitalSet, 6> voxelizer;
  // Mesh<Point> a3DMesh;
  // bool importOK = a3DMesh << "../../gp-cli/build/birds.off";

  trace.beginBlock( "Making digital set from the isovalue matrix. " ); 

  int i = 0;
  for ( DomainConstIterator it = domain.begin(); it != domain.end(); ++it )
  {
    if ( S_iso.row(i)(0) == isovalue ) 
    {
      scene_set.insertNew( *it );
    }
    i++;
  }
  trace.endBlock();

  trace.beginBlock( "Getting the border of the digital set. " ); 
  ObjectType sc( dt6_18, scene_set );

  ObjectType b_scene = sc.border();

  std::vector< ObjectType > objects;
  std::back_insert_iterator< std::vector< ObjectType > > inserter( objects );
  unsigned int nbc = b_scene.writeComponents( inserter );
  std::cout << "number of components " << nbc << std::endl;
  trace.endBlock();

  Board3D<> board;
  board << SetMode3D(domain.className(), "Paving");
  board << CustomColors3D(Color(250, 0,0),Color(250, 0,0));

  // typedef KhalimskySpaceND<3,int> KSpace;
  // KSpace ks;
  board << b_scene;
  // DigitalSet border = b_scene.pointSet();
  // std::cout << "border made" << std::endl;
  // Boundary boundary(ks, border);
  // std::cout << boundary.nbSurfels() << std::endl;

  // std::cout << typeid(boundary).name() << std::endl;


  // for( SurfelConstIterator it = boundary.begin(); it != boundary.end(); ++it )
  //   std::cout << typeid(it).name() << std::endl;

  board.saveOBJ("border.obj");


  trace.beginBlock( "Making a binary image from the iso value matrix. " ); 
  CountedPtr<SH3::BinaryImage> b_image = CountedPtr<SH3::BinaryImage>(new SH3::BinaryImage( domain ));

  int j = 0;
  for(SH3::BinaryImage::Domain::ConstIterator it = b_image->domain().begin();
                      it != b_image->domain().end(); it++)
  {
    b_image->setValue(*it, S_iso.row(j)(0));
    j++;
  }
  trace.endBlock();

  // trace.beginBlock( "Making a grayscale image from the visibilities matrix. " );   
  // CountedPtr<SH3::GrayScaleImage> g_image = CountedPtr<SH3::GrayScaleImage>(new SH3::GrayScaleImage( domain ));

  // int k = 0;
  // for(SH3::GrayScaleImage::Domain::ConstIterator it = g_image->domain().begin();
  //                     it != g_image->domain().end(); it++)
  // {
  //   g_image->setValue(*it, S.row(k)(0));
  //   k++;
  // }
  // trace.endBlock();

  trace.beginBlock( "Making digital surface from binary image. " );
  auto K         = SH3::getKSpace( b_image );
  auto surface   = SH3::makeDigitalSurface( b_image, K, params );
  trace.info() << "#surfels=" << surface->size() << std::endl;
  trace.endBlock();

  trace.beginBlock( "Making polygonal surface. " );
  SH3::Cell2Index c2i;
  SH3::Surfel2Index s2i;

  // auto polysurf = SH3::makePrimalPolygonalSurface( c2i, surface );
  //auto mesh   = SH3::makeMesh( trisurf );
  //trace.info() << "#mesh=" << mesh->nbVertex() << std::endl;
  auto surfels   = SH3::getSurfelRange( surface, params );
  auto pointels  = SH3::getPointelRange( c2i, surface );

  Eigen::MatrixXf V(pointels.size(),3);
  Eigen::MatrixXi F(surfels.size(),4);

  std::cout << "V size " << V.rows() << " " << V.cols() << std::endl;
  std::cout << "F size " << F.rows() << " " << F.cols() << std::endl;

  const SH3::KSpace& space = SH3::refKSpace( surface );
  auto embedder = SH3::getCellEmbedder( surface );
  
  int num_p = 0;
  for(auto&& pointel : pointels)
  {
    SH3::RealPoint p = embedder( pointel );
    V.row(num_p) = Eigen::RowVector3f(p[0],p[1],p[2]);
    num_p++;
  }
  int num_f = 0;
  for(auto s : surfels)
  {
    Eigen::RowVector4i face;
    auto primal_vtcs = SH3::getPointelRange( space, s );
    int f_i = 0;
    for ( auto&& primal_vtx : primal_vtcs )
    {
      face(f_i) = c2i[ primal_vtx ];
      f_i ++;
    }
    F.row(num_f) = face;
    num_f++;
  }

  igl::writeOBJ("igl-test.obj",V,F);
  auto ok = SH3::saveOBJ( surface, SH3::RealVectors(), SH3::Colors(),
                                 "test.obj");
  
  trace.endBlock();

  // trace.beginBlock( "Wrapping a digital set around image. " );
  // typedef functors::IntervalForegroundPredicate<SH3::GrayScaleImage> ThresholdedImage;
  // ThresholdedImage thresholdedImage( *g_image, isovalue, num_views );
  // trace.endBlock();

  // trace.beginBlock( "Extracting boundary by scanning the space. " );
  
  // typedef KhalimskySpaceND<3,int> KSpace;
  // KSpace ks;
  // typedef KSpace::SurfelSet SurfelSet;
  // typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  // typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
  // typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
  // MySurfelAdjacency surfAdj( true ); // interior in all directions.
  // MySetOfSurfels theSetOfSurfels( ks, surfAdj );
  // Surfaces<KSpace>::sMakeBoundary( theSetOfSurfels.surfelSet(),
  //                                    ks, thresholdedImage,
  //                                     domain.lowerBound(),
  //                                     domain.upperBound() );
  // MyDigitalSurface digSurf( theSetOfSurfels );
  // trace.info() << "Digital surface has " << digSurf.size() << " surfels."
  //                << std::endl;
  
  // trace.endBlock();

  // trace.beginBlock( "Displaying everything. " );
  // Board3D<Z3,Z3i::KSpace> board(ks);

  // board << SetMode3D(  ks.unsigns( *digSurf.begin() ).className(), "Basic" );

  // std::string mode = "INNER";
  // typedef MyDigitalSurface::ConstIterator ConstIterator;
  //      if ( mode == "BDRY" )
  //       for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
  //         board << ks.unsigns( *it );
  //      else if ( mode == "INNER" )
  //        for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
  //          board << ks.sCoords( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) );
  //      else if ( mode == "OUTER" )
  //        for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
  //          board << ks.sCoords( ks.sIndirectIncident( *it, ks.sOrthDir( *it ) ) );
    
  // board.saveOBJ("test.obj", false);
  // trace.endBlock();

  return 0;
}
