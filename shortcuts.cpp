#include "DGtal/helpers/Stddefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"

#include <igl/writeDMAT.h>
#include <igl/readDMAT.h>

// #include "DGtal/io/viewers/Viewer3D.h"


using namespace std;
using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;


int main( int argc, char** argv )
{
    auto params = SH3::defaultParameters()
              | SHG3::defaultParameters();

    auto al_capone = SH3::makeBinaryImage( "../DGtal_source/examples/samples/Al.100.vol", params );
    auto K         = SH3::getKSpace( al_capone );
    auto surface   = SH3::makeLightDigitalSurface( al_capone, K, params );
    trace.info() << "#surfels=" << surface->size() << std::endl;

    // params( "faceSubdivision", "Centroid" )( "surfelAdjacency", 1);
    // auto gimage    = SH3::makeGrayScaleImage( "../DGtal_source/examples/samples/lobster.vol" );
    // auto trisurf150= SH3::makeTriangulatedSurface( gimage, params( "thresholdMin", 150 ) );
    // auto trisurf40 = SH3::makeTriangulatedSurface( gimage, params( "thresholdMin", 40 ) );
    // auto mesh150   = SH3::makeMesh( trisurf150 );
    // auto mesh40    = SH3::makeMesh( trisurf40 );
    // trace.info() << "#mesh150=" << mesh150->nbVertex()
    //              << " #mesh40=" << mesh40->nbVertex() << std::endl;


    // typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char> Image3D;
    // Image3D image = DGtal::GenericReader<Image3D>::import(inputFilename);
    // trace.info() << "Image loaded: "<<image<< std::endl;

    // auto canonical_scene = SH3::makeBinaryImage( "../hidden-supports/gpu/build/output_vol_001.pgm3d", params );
    
    Point p1( 0, 0, 0 );
    Point p2( 300,300,215 );
    Domain domain( p1, p2 );
    Point c( 0, 0 );
    // canonical scene
    DigitalSet scene( domain );
    int i = 0;
    for ( DomainConstIterator it = domain.begin(); it != domain.end(); ++it )
    {
        if ( S.row(i) == isovalue ) scene.insertNew( *it );
        i++;
    }
    ObjectType sc( dt6_18, scene );
    // The following line takes almost no time.
    ObjectType sc_clone( sc);

    auto K         = SH3::getKSpace( scene );
    auto surface   = SH3::makeLightDigitalSurface( scene, K, params );
    trace.info() << "#surfels=" << surface->size() << std::endl;

    return 0;
}
