#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/Color.h"
#include "DGtal/images/ImageSelector.h"
// #include "DGtal/io/viewers/Viewer3D.h"

#include "DGtal/helpers/Stddefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"


#include <boost/pending/disjoint_sets.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
 @page volCComponentCounter volCComponentCounter
 
 @brief Counts the number of connected component (same values) in a  volume (Vol) file image.
 @b Usage:  volCComponentCounter [input]
 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                  display this message
  -c [ --connectivity ] arg (=6) object connectivity (6,18,26) (default: 6 )
  -i [ --input ] arg             volume file (Vol) (default: standard input)
 @endcode
 @b Example: 
 @code
 $ volCComponentCounter -i $DGtal/examples/samples/Al.100.vol 
 @endcode
 You should obtain:
@verbatim
New Block [Initial disjoint sets construction]
EndBlock [Initial disjoint sets construction] (515.602 ms)
New Block [Merging neighboring sets]
EndBlock [Merging neighboring sets] (1672.21 ms)
Number of disjoint 6-components = 2
@endverbatim
 
 @see
 @ref volCComponentCounter.cpp
 */


template <typename Rank, typename Parent, typename Image>
void CCCounter(Rank& r, Parent& p, const Image& elements, const unsigned int connectivity)
{
  
  boost::disjoint_sets<Rank,Parent> dsets(r, p);
  trace.beginBlock("Initial disjoint sets construction");
  for(typename Image::Domain::ConstIterator e = elements.domain().begin();
      e != elements.domain().end(); ++e)
    dsets.make_set(*e);
  trace.endBlock();

  trace.beginBlock("Merging neighboring sets");
  typename Image::Point decx(1,0,0);
  typename Image::Point decy(0,1,0);
  typename Image::Point decz(0,0,1);
  
  //Merging process
  for ( typename Image::Domain::ConstIterator e = elements.domain().begin();
        e != elements.domain().end(); ++e )
  {
    if ( elements.domain().isInside( *e + decx ) &&
         ( elements( *e ) == elements( *e + decx ) ) )
      dsets.union_set( *e, *e + decx );

    if ( elements.domain().isInside( *e + decy ) &&
         ( elements( *e ) == elements( *e + decy ) ) )
      dsets.union_set( *e, *e + decy );

    if ( elements.domain().isInside( *e + decz ) &&
         ( elements( *e ) == elements( *e + decz ) ) )
      dsets.union_set( *e, *e + decz );

    if ( connectivity > 6 )
    {
      if ( elements.domain().isInside( *e + decx + decy ) &&
           ( elements( *e ) == elements( *e + decx + decy ) ) )
        dsets.union_set( *e, *e + decx + decy );

      if ( elements.domain().isInside( *e + decx + decz ) &&
           ( elements( *e ) == elements( *e + decx + decz ) ) )
        dsets.union_set( *e, *e + decx + decz );

      if ( elements.domain().isInside( *e + decy + decz ) &&
           ( elements( *e ) == elements( *e + decy + decz ) ) )
        dsets.union_set( *e, *e + decy + decz );

      if ( connectivity == 26 )
        if ( elements.domain().isInside( *e + decy + decz + decx ) &&
             ( elements( *e ) == elements( *e + decy + decz + decx ) ) )
          dsets.union_set( *e, *e + decy + decz + decx );
    }
    }
  trace.endBlock();
  std::cout << "Number of disjoint "<<connectivity<<"-components = "
            <<dsets.count_sets(elements.domain().begin(),
                               elements.domain().end())
            << std::endl;
}



int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("connectivity,c", po::value<unsigned int>()->default_value(6), "object connectivity (6,18,26)"    " (default: 6 )")
    ("input,i", po::value<std::string>(), "volume file (Vol)"    " (default: standard input)");
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }  
  po::notify(vm);    
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input]\n"
                << "Count the number of connected component (same values) in a  volume (Vol) file image\n"
                << general_opt << "\n"
                << "Example : \n \t volCComponentCounter -i $DGtal/examples/samples/Al.100.vol ";
      return 0;
    }
  string inputFilename = vm["input"].as<std::string>();
  unsigned int connectivity = vm["connectivity"].as<unsigned int>();
 
  if ((connectivity != 6) && (connectivity != 18) && (connectivity != 26))
    {
      trace.error() << "Bad connectivity value.";
      trace.info() << std::endl;
      exit(1);
    }

  typedef ImageSelector<Domain, unsigned char>::Type Image;
  // Image image = VolReader<Image>::importVol( inputFilename );
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char> Image3D;
  Image3D image = DGtal::GenericReader<Image3D>::import(inputFilename);

  trace.info() << "Image loaded: "<<image<< std::endl;

  typedef std::map<Point,std::size_t> rank_t; // => order on Element
  typedef std::map<Point,Point> parent_t;
  rank_t rank_map;
  parent_t parent_map;
 
  boost::associative_property_map<rank_t>   rank_pmap(rank_map);
  boost::associative_property_map<parent_t> parent_pmap(parent_map);
 
  CCCounter(rank_pmap, parent_pmap, image, connectivity);
 

  return 0;
}
