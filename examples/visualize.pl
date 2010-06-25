#!/usr/bin/perl -w

use lib '../blib/lib','../blib/arch';

use Algorithm::KMeans qw(kmeans visualize cluster_data_generator);


# IMPORTANT: visualize() will do its job correctly only when
#            is called on the same datafile that you last
#            clustered with a call to kmeans()

# In most cases, you would not change the value of the mask
# between clustering and visualization.  But, if you are
# clustering multi-dimensional data and you wish to
# visualize the projection of of the data on each plane
# separately, you can do so by changing the value of the
# mask as shown below:

$mask = "I111";

visualize( datafile => "mydatafile.dat", 
           mask     =>  $mask );

