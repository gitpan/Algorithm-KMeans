#!/usr/bin/perl -w

#use lib '../blib/lib', '../blib/arch';

use strict;
use Algorithm::KMeans;

my $datafile = "mydatafile1.dat";

# Mask:

# See the comment block in 
#
#          cluster_and_vsualize.pl
#
# for how to set up the mask and what it means.

my $mask = "N111";
my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                        mask     => $mask,
                                        K        => 3,
                                        terminal_output => 1,
                                        do_variance_normalization => 1,
#                                        write_clusters_to_files => 1,
    );

$clusterer->read_data_from_file();

$clusterer->kmeans();

# CLUSTER VISUALIZATION:

# Visualization mask:

# See the comment block in 
#
#          cluster_and_vsualize.pl
#
# for how to set up the mask and what it means.

my $visualization_mask = "111";

$clusterer->visualize_clusters($visualization_mask);

