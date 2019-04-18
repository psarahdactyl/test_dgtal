#include "DGtal/helpers/Stddefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <igl/writeDMAT.h>
#include <igl/readDMAT.h>

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
typedef DigitalTopology< Adj6, Adj18 > DT6_18;
typedef Object<DT6_18, DigitalSet> ObjectType;

namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("inputFilename,f", po::value<std::string>(), "input matrix of visibilities (.dmat)" )
  ("outputFilename,o", po::value<std::string>(), "output obj file (.obj)" )
  ("isovalue,i",  po::value<int>()->default_value(1), "isovalue (not normalized) of surface to extract in [1, num_views]" );

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
  int isovalue = vm["isovalue"].as<int>();

  auto params = SH3::defaultParameters()
            | SHG3::defaultParameters();

    // auto al_capone = SH3::makeBinaryImage( "../../DGtal_source/examples/samples/Al.100.vol", params );
    // auto K         = SH3::getKSpace( al_capone );
    // auto surface   = SH3::makeLightDigitalSurface( al_capone, K, params );
    // trace.info() << "#surfels=" << surface->size() << std::endl;

    // auto canonical_scene = SH3::makeBinaryImage( "../hidden-supports/gpu/build/output_vol_001.pgm3d", params );
  Eigen::MatrixXi S;
  igl::readDMAT(inputFilename, S);
  // int isovalue = 61;
  int num_views = 500;
  Eigen::MatrixXi S_iso = (S.array() >= isovalue).select(isovalue, S.array()-S.array());
  std::cout << S_iso.rows() << " " << S_iso.cols() << std::endl;

  Point p1( 1, 1, 1 );
  Point p2( 100,90,29 );
  Domain domain( p1, p2 );
  // canonical scene
  DigitalSet scene_set( domain );

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

  ObjectType b_scene = sc.border(); // one component

  std::vector< ObjectType > objects;
  std::back_insert_iterator< std::vector< ObjectType > > inserter( objects );
  unsigned int nbc = b_scene.writeComponents( inserter );
  std::cout << "number of components " << nbc << std::endl;
  trace.endBlock();

  Board3D<> board;
  board << SetMode3D(domain.className(), "Paving");
  board << CustomColors3D(Color(250, 0,0),Color(250, 0,0));

  for(auto it : b_scene.pointSet())
    board << it;

  board.saveOBJ(outputFilename);


  // trace.beginBlock( "Making a binary image from the iso value matrix. " ); 
  // CountedPtr<SH3::BinaryImage> b_image = CountedPtr<SH3::BinaryImage>(new SH3::BinaryImage( domain ));

  // int j = 0;
  // for(SH3::BinaryImage::Domain::ConstIterator it = b_image->domain().begin();
  //                     it != b_image->domain().end(); it++)
  // {
  //   b_image->setValue(*it, S_iso.row(j)(0));
  //   j++;
  // }
  // trace.endBlock();

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

  // trace.beginBlock( "Making digital surface from binary image. " );
  // SH3::Surfel2Index s2i;
  // auto K         = SH3::getKSpace( b_image );
  // auto surface   = SH3::makeDigitalSurface( b_image, K, params );
  // trace.info() << "#surfels=" << surface->size() << std::endl;

  // //auto trisurf = SH3::makeTriangulatedSurface( s2i, surface );
  // //auto mesh   = SH3::makeMesh( trisurf );
  // //trace.info() << "#mesh=" << mesh->nbVertex() << std::endl;
  
  // auto ok      = SH3::saveOBJ( surface, SH3::RealVectors(), SH3::Colors(),
  //                                "test.obj");
  // trace.endBlock();

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
