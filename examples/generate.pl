#!/usr/bin/perl -w

use lib '../blib/lib','../blib/arch';

use Algorithm::KMeans qw(kmeans visualize cluster_data_generator);


# The Parameter File:

# How the synthetic data is generated for clustering is
# controlled entirely by the input_parameter_file keyword in
# the function call shown below.  The mean vector and
# covariance matrix entries in file must be according to the
# syntax shown in the example param.txt file.  It is best to
# edit this file as needed for the purpose of data
# generation.

cluster_data_generator( input_parameter_file => "param.txt",
                        output_datafile => "mynewdatafile.dat",
                        number_data_points_per_cluster => 10 );

