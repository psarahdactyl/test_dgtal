#include <DGtal/base/Common.h>
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/helpers/StdDefs.h"

using namespace DGtal;

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

int main(int argc, char** argv)
{

  Point p1( -50, -50, -50 );
  Point p2( 50, 50, 50 );
  Domain domain( p1, p2 );
  Point c( 0, 0 );
  // diamond of radius 30
  DigitalSet diamond_set( domain );
  for ( DomainConstIterator it = domain.begin(); it != domain.end(); ++it )
  {
    if ( (*it - c ).norm1() <= 30 ) diamond_set.insertNew( *it );
  }
  ObjectType diamond( dt6_18, diamond_set );
  // The following line takes almost no time.
  ObjectType diamond_clone( diamond );
  // Since one of the objects is modified, the set is duplicated at the following line
  diamond_clone.pointSet().erase( c );
  ObjectType bdiamond = diamond.border(); // one component
  std::vector< ObjectType > objects;
  std::back_insert_iterator< std::vector< ObjectType > > inserter( objects );
  // nbc == 1 since the boundary of the diamond is connected.
  unsigned int nbc = bdiamond.writeComponents( inserter );
  std::cout << nbc << std::endl;

  ObjectType bdiamond_clone = diamond_clone.border(); // two components

  DGtal::trace.info() << "Helloworld from DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;
  
  return 0;
}
