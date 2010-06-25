#!/usr/bin/perl -w

use lib '../blib/lib','../blib/arch';

use Algorithm::KMeans qw(kmeans visualize cluster_data_generator);


my $parameter_file = "param.txt";
my $datafile = "mydatafile.dat";

# The Parameter File:

# How the synthetic data is generated for clustering is
# controlled entirely by the input_parameter_file keyword in
# the function call shown below.  The mean vector and
# covariance matrix entries in file must be according to the
# syntax shown in the example param.txt file.  It is best to
# edit this file as needed for the purpose of data
# generation.

cluster_data_generator( input_parameter_file => $parameter_file,
                        output_datafile => $datafile,
                        number_data_points_per_cluster => 60 );

# Mask:

# The mask tells the module which columns of the data file
# are are to be used for clustering, which columns are to be
# ignored and which column contains the symbolic ID tag for
# a data point.  If the ID is in column 1 and you are
# clustering 3D data, the mast would be "I111".  Note the
# first character in the mask in this case is `I' for "Id".
# If, on the other hand, you wanted to ignore the first data
# coordinate for clustering, the mask would be "I011".  The
# symbolic ID can be in any column --- you just have to
# place the character `I' at the right place:

my $mask = "I111";

kmeans( datafile => $datafile,
        mask     =>  $mask,
        terminal_output => 1,  
        K => 3 );

# In most cases, you would not change the value of the mask
# between clustering and visualization.  But, if you are
# clustering multi-dimensional data and you wish to
# visualize the projection of of the data on each plane
# separately, you can do so by changing the value of the
# mask as shown below:

$mask = "I111";

visualize( datafile => $datafile,
           mask     =>  $mask );

