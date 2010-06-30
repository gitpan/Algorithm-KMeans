package Algorithm::KMeans;

#---------------------------------------------------------------------------
# Copyright (c) 2010 Avinash Kak. All rights reserved.
# This program is free software.  You may modify and/or
# distribute it under the same terms as Perl itself.
# This copyright notice must remain attached to the file.
#
# Algorithm::KMeans is a pure Perl implementation for
# clustering multi-dimensional data.
#---------------------------------------------------------------------------

use 5.10.0;
use strict;
use warnings;
use Carp;
use File::Basename;
use Math::Random;
use Graphics::GnuplotIF;

our $VERSION = '1.1.1';

# from perl docs:
my $_num_regex =  '^[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$'; 

# Constructor:
sub new { 
    my ($class, %args) = @_;
    bless {
        _datafile         =>   $args{datafile} || croak("datafile required"),
        _mask             =>   $args{mask}     || croak("mask required"),
        _K                =>   $args{K}        || 0,
        _K_best           =>   'unknown',
        _K_min            =>   $args{Kmin} || 'unknown',
        _K_max            =>   $args{Kmax} || 'unknown',
        _terminal_output  =>   $args{terminal_output} || 0,
        _data             =>   {},
        _data_id_tags     =>   [],
        _N                =>   0,
        _QoC_values       =>   {},
        _clusters         =>   [],
        _cluster_centers  =>   [],
        _data_dimensions  =>   0,
        _clusters_2_files =>  $args{write_clusters_to_files} || 0,
    }, $class;
}

sub read_data_from_file {
    my $self = shift;
    my $datafile = $self->{_datafile};
    my $mask = $self->{_mask};
    my @mask = split //, $mask;

    $self->{_data_dimensions} = scalar grep {$_ eq '1'} @mask;
    print "data dimensionality:  $self->{_data_dimensions} \n"
	if $self->{_terminal_output};

    open INPUT, $datafile
        or die "unable to open file $datafile: $!\n";
    chomp( my @raw_data = <INPUT> );
    close INPUT;

    # Transform strings into number data
    foreach my $record (@raw_data) {
        next if $record =~ /^#/;
        my @data_fields;
        my @fields = split /\s+/, $record;
        die "\nABORTED: Mask size does not correspond to row record size\n" 
            if $#fields != $#mask;
        my $record_id;
        foreach my $i (0..@fields-1) {
            if ($mask[$i] eq '0') {
                next;
            } elsif ($mask[$i] eq 'N') {
                $record_id = $fields[$i];
            } elsif ($mask[$i] eq '1') {
                push @data_fields, $fields[$i];
            } else {
                die "misformed mask for reading the data file\n";
            }
        }
        my @nums = map {/$_num_regex/;$_} @data_fields;
        $self->{_data}->{ $record_id } = \@nums;
    }
    my @all_data_ids = keys %{$self->{_data}};
    $self->{_data_id_tags} = \@all_data_ids;
    $self->{_N} = scalar @all_data_ids;
    die "You need at least 8 data samples. The number of data points " .
        "must satisfy the relation N = 2xK**2 where K is the " .
        "number of clusters.  The smallest value for K is 2.\n"
        if $self->{_N} <= 8;
}

sub kmeans {
    my $self = shift;
    my $K = $self->{_K};
    if ( ($K == 0) ||
              ( ($self->{_K_min} ne 'unknown') &&
                ($self->{_K_max} ne 'unknown') ) ) {
        $self->iterate_through_K();
    } elsif ( $K =~ /\d+/) {
        $self->cluster_for_fixed_K_multiple_tries($K);
    } else {
        die "Incorrect call syntax used.  See documentation.\n";
    }
    if ((defined $self->{_clusters}) && (defined $self->{_cluster_centers})){
        $self->write_clusters_to_files( $self->{_clusters} )
            if $self->{_clusters_2_files};
        return ($self->{_clusters}, $self->{_cluster_centers});
    } else {
        die "Clustering failed.  Try again.";
    }
}

sub cluster_for_fixed_K_multiple_tries {
    my $self = shift;
    my $K = shift;
    my @all_data_ids = @{$self->{_data_id_tags}};
    my $QoC;
    my @QoC_values;
    my @array_of_clusters;
    my @array_of_cluster_centers;
    my $clusters;
    my $new_clusters;
    my $cluster_centers;
    my $new_cluster_centers;
    print "Clustering for K = $K\n" if $self->{_terminal_output};
    foreach my $trial (1..20) {
        my ($new_clusters, $new_cluster_centers) =
                              $self->cluster_for_given_K($K);

        my $newQoC = $self->cluster_quality( $new_clusters, 
                                             $new_cluster_centers );
        if ( (!defined $QoC) || ($newQoC < $QoC) ) {
            $QoC = $newQoC;
            $clusters = deep_copy( $new_clusters );
            $cluster_centers = deep_copy( $new_cluster_centers );
        } 
    }
    $self->{_clusters} = $clusters;
    $self->{_cluster_centers} = $cluster_centers;  
    $self->{_QoC_values}->{"$K"} = $QoC; 
    if ($self->{_terminal_output}) {
        print "\nDisplaying final clusters for best K (= $K) :\n";
        display_clusters( $clusters );
        display_cluster_centers( $cluster_centers );
        print "QoC value: $QoC\n";
#       print "QoC values array for different K: " . 
#               "@{[ map {my $x = sprintf('%.4f', $_); $x} @QoC_values ]}\n";
    }
}

# The following subroutine is the top-level routine to call
# if you want the system to figure out on its own what value
# to use for K, the number of clusters.  It will try every
# possible value of K between 2 and the maximum possible
# depending on the number of data points available. For
# example, if the number of data points is 10,000, it will
# try all possible values of K between 2 and 70. For how the
# maximum value is set for K, see the comments made under
# Description.  Note also how this function makes 20
# different tries for each value of K as a defense against
# the problem of the final result corresponding to some
# local minimum in the values of the QoC metric.  Out of
# these 20 tries for each K, it retains the clusters and the
# cluster centers for only that try that yields the smallest
# value for the QoC metric.  After estimating the "best" QoC
# values for all possible K in this manner, it then finds
# the K for which the QoC is the minimum.  This is taken to
# be the best value for K.  Finally, the output clusters are
# written out to separate files.
sub iterate_through_K {
    my $self = shift;
    my @all_data_ids = @{$self->{_data_id_tags}};
    my $N = $self->{_N};
    die "You need at least 8 data samples. The number of data points " .
        "must satisfy the relation N = 2xK**2 where K is the " .
        "number of clusters.  The smallest value for K is 2.\n"
        if $N <= 8;
    my $K_statistical_max = int( sqrt( $N / 2.0 ) );
    my $Kmin = $self->{_K_min} eq 'unknown' 
                          ? 2
                          : $self->{_K_min};
    my $Kmax = $self->{_K_max} eq 'unknown' 
                          ? int( sqrt( $N / 2.0 ) )
                          : $self->{_K_max};
    die  "\n\nYour Kmax value is too high.  Ideally, it should " .
         "not exceed sqrt(N/2) where N is the number of data points"
          if $Kmax > $K_statistical_max;

    print "Value of Kmax is: $Kmax\n" if $self->{_terminal_output};
    my @QoC_values;
    my @array_of_clusters;
    my @array_of_cluster_centers;
    foreach my $K ($Kmin..$Kmax) {
        my $QoC;
        my $clusters;
        my $cluster_centers;
        print "Clustering for K = $K\n" if $self->{_terminal_output};
        foreach my $trial (1..20) {
            my ($new_clusters, $new_cluster_centers) = 
                           $self->cluster_for_given_K($K);
            my $newQoC = $self->cluster_quality( $new_clusters, 
                                                 $new_cluster_centers );
            if ( (!defined $QoC) || ($newQoC < $QoC) ) {
                $QoC = $newQoC;
                $clusters = deep_copy( $new_clusters );
                $cluster_centers = deep_copy( $new_cluster_centers );
            } 
        }
        push @QoC_values, $QoC;
        push @array_of_clusters, $clusters;
        push @array_of_cluster_centers, $cluster_centers;
    }
    my ($min, $max) = minmax( \@QoC_values );
    die "Unsuccessful. Try again.\n" if ($max - $min ) < 0.00001;
    my $K_best_relative_to_Kmin = get_index_at_value($min, \@QoC_values );
    my $K_best = $K_best_relative_to_Kmin + $Kmin;

    if ($self->{_terminal_output}) {
        print "\nDisplaying final clusters for best K (= $K_best) :\n";
        display_clusters( $array_of_clusters[$K_best_relative_to_Kmin] );
        display_cluster_centers( 
                   $array_of_cluster_centers[$K_best_relative_to_Kmin] );
        print "\nBest clustering achieved for K=$K_best with QoC = $min\n";
        print "QoC values array for different K starting " .
                            "with K=$Kmin:  @QoC_values\n";
    }
    $self->{_K_best} = $K_best;
    foreach my $i (0..@QoC_values-1) {
        my $k = $i + $Kmin;
        $self->{_QoC_values}->{"$k"} = $QoC_values[$i]; 
    }
    $self->{_clusters} = $array_of_clusters[$K_best_relative_to_Kmin];
    $self->{_cluster_centers} =  
                $array_of_cluster_centers[$K_best_relative_to_Kmin];
}

# This is the function to call if you already know what
# value you want to use for K, the number of expected
# clusters.  The purpose of this function is to do the
# initialization of the cluster centers and to carry out the
# initial assignment of the data to the clusters with the
# initial cluster centers.  The initialization consists of 3
# steps: Construct a random sequence of K integers between 0
# and N-1 where N is the number of data points to be
# clustered; 2) Call get_initial_cluster_centers() to index
# into the data array with the random integers to get a list
# of K data points that would serve as the initial cluster
# centers; and (3) Call assign_data_to_clusters_initial() to
# assign the rest of the data to each of the K clusters on
# the basis of the proximity to the cluster centers.
sub cluster_for_given_K {
    my $self = shift;
    my $K = shift;
    my @all_data_ids = @{$self->{_data_id_tags}};
    my @cluster_center_indices =  $self->initialize_cluster_centers($K);
    my $cluster_centers = $self->get_initial_cluster_centers( 
                                         \@cluster_center_indices );

    my $clusters = $self->assign_data_to_clusters_initial($cluster_centers);  

    my $cluster_nonexistant_flag = 0;
    foreach my $trial (0..2) {
        ($clusters, $cluster_centers) =
                         $self->assign_data_to_clusters( $clusters, $K );
        my $num_of_clusters_returned = @$clusters;
        foreach my $cluster (@$clusters) {
            $cluster_nonexistant_flag = 1 if ((!defined $cluster) 
                                             ||  (@$cluster == 0));
        }
        last unless $cluster_nonexistant_flag;
    }
    return ($clusters, $cluster_centers);
}

# Returns a set of K random integers.  These serve as
# indices to reach into the data array.  A data element
# whose index is one of the random numbers returned by this
# routine serves as an initial cluster center.  Note the
# quality check it runs on the list of K random integers
# constructed.  We first make sure that all K random
# integers are different.  Subsequently, we carry out a
# quality assessment of the K random integers constructed.
# This quality measure consists of the ratio of the values
# spanned by the random integers to the value of N, the
# total number of data points to be clustered.  Currently,
# if this ratio is less than 0.3, we discard the K integers
# and try again.
sub initialize_cluster_centers {
    my $self = shift;
    my $K = shift;
    my $data_store_size = $self->{_N};
    my @cluster_center_indices;
    while (1) {
        foreach my $i (0..$K-1) {
            $cluster_center_indices[$i] = int rand( $data_store_size );
            next if $i == 0;
            foreach my $j (0..$i-1) {
                while ( $cluster_center_indices[$j] == 
                        $cluster_center_indices[$i] ) {
                    my $old = $cluster_center_indices[$i];
                    $cluster_center_indices[$i] = int rand($data_store_size);
                }
            }
        }
        my ($min,$max) = minmax(\@cluster_center_indices );
        my $quality = ($max - $min) / $data_store_size;
        last if $quality > 0.3;
    }
    return @cluster_center_indices;
}

# This routine merely reaches into the data array with the
# random integers, as constructed by the previous routine,
# serving as indices and fetching values corresponding to
# those indices.  The fetched data samples serve as the
# initial cluster centers.
sub get_initial_cluster_centers {
    my $self = shift;
    my @cluster_center_indices = @{shift @_};
    my @result;
    foreach my $i (@cluster_center_indices) {    
        my $tag = $self->{_data_id_tags}[$i];     
        push @result, $self->{_data}->{$tag};
    }
    return \@result;
}

# The purpose of this routine is to form initial clusters by
# assigning the data samples to the initial clusters formed
# by the previous routine on the basis of the best proximity
# of the data samples to the different cluster centers.
sub assign_data_to_clusters_initial {
    my $self = shift;
    my @cluster_centers = @{ shift @_ };
    my @clusters;
    foreach my $ele (@{$self->{_data_id_tags}}) {
        my $best_cluster;
        my @dist_from_clust_centers;
        foreach my $center (@cluster_centers) {
            push @dist_from_clust_centers, $self->distance($ele, $center);
        }
        my ($min, $best_center_index) = minimum( \@dist_from_clust_centers );
        push @{$clusters[$best_center_index]}, $ele;
    }
    return \@clusters;
}    

# This is the main routine that along with the
# update_cluster_centers() routine constitute the two key
# steps of the K-Means algorithm.  In most cases, the
# infinite while() loop will terminate automatically when
# the cluster assignments of the data points remain
# unchanged. For the sake of safety, we keep track of the
# number of iterations. If this number reaches 100, we exit
# the while() loop anyway.  In most cases, this limit will
# not be reached.
sub assign_data_to_clusters {
    my $self = shift;
    my $clusters = shift;
    my $K = shift;
    my $final_cluster_centers;
    my $iteration_index = 0;
    while (1) {
        my $new_clusters;
        my $assignment_changed_flag = 0;
        my $current_cluster_center_index = 0;
        my $cluster_size_zero_condition = 0;
        my $how_many = @$clusters;
        my $cluster_centers = $self->update_cluster_centers( $clusters );
        $iteration_index++;
        foreach my $cluster (@$clusters) {
            my $current_cluster_center = 
                          $cluster_centers->[$current_cluster_center_index];
            foreach my $ele (@$cluster) {
                my @dist_from_clust_centers;
                foreach my $center (@$cluster_centers) {
                    push @dist_from_clust_centers, 
                               $self->distance($ele, $center);
                }
                my ($min, $best_center_index) = 
                              minimum( \@dist_from_clust_centers );
                my $best_cluster_center = 
                                 $cluster_centers->[$best_center_index];
                if (vector_equal($current_cluster_center, 
                                         $best_cluster_center)){
                    push @{$new_clusters->[$current_cluster_center_index]}, 
                                  $ele;
                } else {
                    $assignment_changed_flag = 1;             
                    push @{$new_clusters->[$best_center_index]}, $ele;
                }
            }
            $current_cluster_center_index++;
        }
        # Now make sure that we still have K clusters since K is fixed:
        next if ((@$new_clusters != @$clusters) && ($iteration_index < 100));
        # Now make sure that none of the K clusters is an empty cluster:
        foreach my $newcluster (@$new_clusters) {
            $cluster_size_zero_condition = 1 if ((!defined $newcluster) 
                                             or  (@$newcluster == 0));
        }
        next if (($cluster_size_zero_condition) && ($iteration_index < 100));
        last if $iteration_index == 100;
        # Now do a deep copy of new_clusters into clusters
	$clusters = deep_copy( $new_clusters );
        last if $assignment_changed_flag == 0;
    }
    $final_cluster_centers = $self->update_cluster_centers( $clusters );
    return ($clusters, $final_cluster_centers);
}

# After each new assignment of the data points to the
# clusters on the basis of the current values for the
# cluster centers, we call the routine shown here for
# updating the values of the cluster centers.
sub update_cluster_centers {
    my $self = shift;
    my @clusters = @{ shift @_ };
    my @new_cluster_centers;
    foreach my $cluster (@clusters) {
        die "Cluster became empty --- untenable condition " .
            "for a given K.  Try again. \n" if !defined $cluster;
        my $cluster_size = @$cluster;
        die "Cluster size is zero --- untenable.\n" if $cluster_size == 0;
        my @new_cluster_center = @{$self->add_point_coords( $cluster )};
        @new_cluster_center = map {my $x = $_/$cluster_size; $x} 
                                  @new_cluster_center;
        push @new_cluster_centers, \@new_cluster_center;
    }        
    return \@new_cluster_centers;
}

# The following function returns the value of QoC for a
# given partitioning of the data into K clusters.  It
# calculates two things: the average value for the distance
# between a data point and the center of the cluster in
# which the data point resides, and the average value for
# the distances between the cluster centers.  We obviously
# want to minimize the former and maximize the latter.  All
# of the "from center" distances within each cluster are
# stored in the variable $sum_of_distances_for_one_cluster.
# When this variable, after it is divided by the number of
# data elements in the cluster, is summed over all the
# clusters, we get a value that is stored in
# $avg_dist_for_cluster.  The inter-cluster-center distances
# are stored in the variable $inter_cluster_center_dist.
sub cluster_quality {
    my $self = shift;
    my $clusters = shift;
    my $cluster_centers = shift;
    my $K = @$cluster_centers;          # Number of clusters
    my $cluster_radius = 0;
    foreach my $i (0..@$clusters-1) {
        my $sum_of_distances_for_one_cluster = 0;
        foreach my $ele (@{$clusters->[$i]}) {
            $sum_of_distances_for_one_cluster += 
                $self->distance( $ele, $cluster_centers->[$i] );
        }
       $cluster_radius += 
           $sum_of_distances_for_one_cluster / @{$clusters->[$i]};
    }
    my $inter_cluster_center_dist = 0;
    foreach my $i (0..@$cluster_centers-1) {
        foreach my $j (0..@$cluster_centers-1) {
            $inter_cluster_center_dist += 
              $self->distance2( $cluster_centers->[$i], 
                                $cluster_centers->[$j] );
        }
    }
    my $avg_inter_cluster_center_dist = $inter_cluster_center_dist /
                    ( $K * ($K-1) / 2.0 );
    return $cluster_radius / $avg_inter_cluster_center_dist;
}

# The following routine is for computing the distance
# between a data point specified by its symbolic name in the
# master datafile and a point (such as the center of a
# cluster) expressed as a vector of coordinates:
sub distance {
    my $self = shift;
    my $ele1_id = shift @_;            # symbolic name of data sample
    my @ele1 = @{$self->{_data}->{$ele1_id}};
    my @ele2 = @{shift @_};
    die "wrong data types for distance calculation\n" if @ele1 != @ele2;
    my $how_many = @ele1;
    my $squared_sum = 0;
    foreach my $i (0..$how_many-1) {
        $squared_sum += ($ele1[$i] - $ele2[$i])**2;
    }    
    my $dist = sqrt $squared_sum;
    return $dist;
}

# The following routine does the same as above but now both
# arguments are expected to be arrays of numbers:
sub distance2 {
    my $self = shift;
    my @ele1 = @{shift @_};
    my @ele2 = @{shift @_};
    die "wrong data types for distance calculation\n" if @ele1 != @ele2;
    my $how_many = @ele1;
    my $squared_sum = 0;
    foreach my $i (0..$how_many-1) {
        $squared_sum += ($ele1[$i] - $ele2[$i])**2;
    }    
    return sqrt $squared_sum;
}

sub write_clusters_to_files {
    my $self = shift;
#    my @clusters = @{shift @_};        
    my @clusters = @{$self->{_clusters}};
    unlink glob "Cluster*.dat";
    foreach my $i (1..@clusters) {
        my $filename = "Cluster" . $i . ".dat";
        print "Writing cluster $i to file $filename\n"
                            if $self->{_terminal_output};
        open FILEHANDLE, "| sort > $filename"
            or die "Unable to pen file: $!";
        foreach my $ele (@{$clusters[$i-1]}) {        
            print FILEHANDLE "$ele\n";
        }
        close FILEHANDLE;
    }
}

sub get_K_best {
    my $self = shift;
    croak "You need to run the clusterer with K=0 option " .
          "before you can call this method" 
                            if $self->{_K_best} eq 'unknown';
    print "The best value of K: $self->{_K_best}\n"
                     if $self->{_terminal_output};
    return $self->{_K_best};
}

sub show_QoC_values {
    my $self = shift;
    croak "You need to run the clusterer with K=0 option " .
          "before you can call this method" 
                            if $self->{_K_best} eq 'unknown';
    print "Show below are K on the left and the QoC values on the right\n";
    foreach my $key (sort keys %{$self->{_QoC_values}} ) {
        print " $key  =>  $self->{_QoC_values}->{$key}\n";
    }
}

sub DESTROY {
    unlink "__temp_" . basename($_[0]->{_datafile});
}

##################  Cluster Visualization Code #################

#  Note that it makes sense to call visualize() only AFTER
#  you have called kmeans() The latter function deposits the
#  symbolic labels of the data points in the different
#  clusters in different files that are named ClusterX.dat
#  for values of X starting with X=0.
#
#
#  IMPORTANT:  visualize() can only be run on the same 
#              source data file on which you last ran
#              kmeans()
#
#  The visualize() implementation automatically figures out
#  whether it should do a 2D plot or a 3D plot.  If the
#  number of on bits in the mask that is supplied as one of
#  the arguments is greater than 2, it does a 3D plot for
#  the first three data coordinates.  That is, the clusters
#  will be displayed in the 3D space formed by the first
#  three data coordinates. On the other hand, if the number
#  of on bits in the mask is exactly 2, it does a 2D plot.
#  Should it happen that only one on bit is specified for
#  the mask, visualize() aborts.
#
#  The visualization code consists of reading the original
#  source data file (which is presumed to have the same
#  structure as in the call to the clustering routine) and
#  creating the hash %all_data.  Subsequently, we read each
#  ClusterX.dat file that was created by the kmeans()
#  subroutine.  Note that the cluster files contain only the
#  symbolic names for the individual records in the source
#  data file.  We therefore next reach into the %all_data
#  hash table and get the data coordinates associated with
#  each symbolic label in the ClusterX.dat file.  The
#  numerical data thus generated is then written out to a
#  temp file.  When doing so we must remember to insert TWO
#  BLANK LINES between the data blocks corresponding to the
#  different ClusterX.dat files.  This constraint is imposed
#  on us by Gnuplot when plotting data from the same file.
#  Note that we want to use different point styles for the
#  data points in different cluster files.
#
#  Subsequently, we call upon the Perl interface provided by
#  the Graphics::GnuplotIF module to plot the data clusters.
sub visualize {
    my $self = shift;
    my $v_mask;
    my $pause_time;
    if (@_ == 1) {
        $v_mask = shift || croak "visualization mask missing";
    } elsif (@_ == 2) {
        $v_mask = shift || croak "visualization mask missing";    
        $pause_time = shift;
    } else {
        croak "visualize() called with wrong args";
    }
    my $master_datafile = $self->{_datafile};

    my @v_mask = split //, $v_mask;
    my $visualization_mask_width = @v_mask;
    my $original_data_mask = $self->{_mask};
    my @mask = split //, $original_data_mask;
    my $data_field_width = scalar grep {$_ eq '1'} @mask;    

    croak "\n\nABORTED: The width of the visualization mask (including " .
          "all its 1s and 0s) must equal the width of the original mask " .
          "used for reading the data file (counting only the 1's)"
          if $visualization_mask_width != $data_field_width;

    my $visualization_data_field_width = scalar grep {$_ eq '1'} @v_mask;

    my %visualization_data;

    while ( my ($record_id, $data) = each %{$self->{_data}} ) {
        my @fields = @$data;
        croak "\nABORTED: Visulization mask size exceeds data record size\n" 
            if $#v_mask > $#fields;
        my @data_fields;
        foreach my $i (0..@fields-1) {
            if ($v_mask[$i] eq '0') {
                next;
            } elsif ($v_mask[$i] eq '1') {
                push @data_fields, $fields[$i];
            } else {
                croak "Misformed visualization mask. It can only have 1s and 0s\n";
            }
        }
        $visualization_data{ $record_id } = \@data_fields;
    }

    my @all_data_ids = @{$self->{_data_id_tags}};

    my $K = scalar @{$self->{_clusters}};

    my $filename = basename($master_datafile);
    my $temp_file = "__temp_" . $filename;
    unlink $temp_file if -e $temp_file;
    open OUTPUT, ">$temp_file"
           or die "Unable to open a temp file in this directory: $!\n";
    foreach my $cluster (@{$self->{_clusters}}) {
        foreach my $item (@$cluster) {
            print OUTPUT "@{$visualization_data{$item}}";
            print OUTPUT "\n";
        }
        print OUTPUT "\n\n";
    }
    close OUTPUT;

    my $plot;
    if (!defined $pause_time) {
        $plot = Graphics::GnuplotIF->new( persist => 1 );
    } else {
        $plot = Graphics::GnuplotIF->new();
    }

    $plot->gnuplot_cmd( "set noclip" );
    $plot->gnuplot_cmd( "set pointsize 2" );

    my $arg_string = "";
    if ($visualization_data_field_width > 2) {
        foreach my $i (0..$K-1) {
            my $j = $i + 1;
            $arg_string .= "\"$temp_file\" index $i using 1:2:3 notitle with points lt $j pt $j, ";
        }
    } elsif ($visualization_data_field_width == 2) {
        foreach my $i (0..$K-1) {
            my $j = $i + 1;
            $arg_string .= "\"$temp_file\" index $i using 1:2 notitle with points lt $j pt $j, ";
        }
    } elsif ($visualization_data_field_width == 1 ) {
        foreach my $i (0..$K-1) {
            my $j = $i + 1;
            $arg_string .= "\"$temp_file\" index $i using 1 notitle with points lt $j pt $j, ";
        }
    }

    $arg_string = $arg_string =~ /^(.*),[ ]+$/;
    $arg_string = $1;

    if ($visualization_data_field_width > 2) {
        $plot->gnuplot_cmd( "splot $arg_string" );
        $plot->gnuplot_pause( $pause_time ) if defined $pause_time;
    } elsif ($visualization_data_field_width == 2) {
        $plot->gnuplot_cmd( "plot $arg_string" );
        $plot->gnuplot_pause( $pause_time ) if defined $pause_time;
    } elsif ($visualization_data_field_width == 1) {
        die "No provision for plotting 1-D data\n";
    }
}

###########  Generating Synthetic Data for Clustering  ############

#  The data generated corresponds to a multivariate
#  distribution.  The mean and the covariance of each
#  Gaussian in the distribution are specified individually
#  in a parameter file.  See the example parameter file
#  param.txt in the examples directory.  Just edit this
#  file for your own needs.
#
#  The multivariate random numbers are generated by calling
#  the Math::Random module.  As you would expect, that
#  module will insist that the covariance matrix you
#  specify be symmetric and positive definite.
sub cluster_data_generator {
    my $class = shift;
    croak "illegal call of a class method" 
        unless $class eq 'Algorithm::KMeans';
    my %args = @_;
    my $input_parameter_file = $args{input_parameter_file};
    my $output_file = $args{output_datafile};
    my $N = $args{number_data_points_per_cluster};

    my @all_params;
    my $param_string;
    if (defined $input_parameter_file) {
        open INPUT, $input_parameter_file
            || "unable to open parameter file: $!";
        @all_params = <INPUT>;
        @all_params = grep { $_ !~ /^[ ]*#/ } @all_params;
        chomp @all_params;
        $param_string = join ' ', @all_params;
    } else {
        # Just for testing. Unsed in t/test.t
        $param_string = "cluster 5 0 0  1 0 0 0 1 0 0 0 1 " .
                        "cluster 0 5 0  1 0 0 0 1 0 0 0 1 " .
                        "cluster 0 0 5  1 0 0 0 1 0 0 0 1";
    }

    my @cluster_strings = split /[ ]*cluster[ ]*/, $param_string;
    @cluster_strings = grep  $_, @cluster_strings;

    my $K = @cluster_strings;
    die "Too many clusters requested" if $K > 12;
    my @point_labels = ('a'..'z');

    print "Number of Gaussians used for the synthetic data: $K\n";

    my @means;
    my @covariances;
    my $data_dimension;
    foreach my $i (0..$K-1) {
        my @num_strings = split /  /, $cluster_strings[$i];
        my @cluster_mean = map {/$_num_regex/;$_} split / /, $num_strings[0];
        $data_dimension = @cluster_mean;
        push @means, \@cluster_mean;
        my @covariance_nums = map {/$_num_regex/;$_} split / /, $num_strings[1];
        die "dimensionality error" if @covariance_nums != 
                                      ($data_dimension ** 2);
        
        my $cluster_covariance;
        foreach my $j (0..$data_dimension-1) {
            foreach my $k (0..$data_dimension-1) {        
                $cluster_covariance->[$j]->[$k] = 
                         $covariance_nums[$j*$data_dimension + $k];
            }
        }
        push @covariances, $cluster_covariance;
    }

    random_seed_from_phrase( 'hellojello' );

    my @data_dump;
    foreach my $i (0..$K-1) {
        my @m = @{shift @means};
        my @covar = @{shift @covariances};
        my @new_data = Math::Random::random_multivariate_normal( $N, 
                                                           @m, @covar );
        my $p = 0;
        my $label = $point_labels[$i];
        @new_data = map {unshift @$_, $label.$i; $i++; $_} @new_data;
        push @data_dump, @new_data;     
    }

    fisher_yates_shuffle( \@data_dump );

    open OUTPUT, ">$output_file";
    foreach my $ele (@data_dump) {
        foreach my $coord ( @$ele ) {
            print OUTPUT "$coord ";
        }
        print OUTPUT "\n";
    }
    print "Data written out to file $output_file\n";
    close OUTPUT;
}

sub add_point_coords {
    my $self = shift;
    my @arr_of_ids = @{shift @_};      # array of data element names
    my @result;
    my $data_dimensionality = $self->{_data_dimensions};
    foreach my $i (0..$data_dimensionality-1) {
        $result[$i] = 0.0;
    }
    foreach my $id (@arr_of_ids) {
        my $ele = $self->{_data}->{$id};
        my $i = 0;
        foreach my $component (@$ele) {
            $result[$i] += $component;
            $i++;
        }
    }
    return \@result;
}

######################   Support Routines  ########################

sub get_index_at_value {
    my $value = shift;
    my @array = @{shift @_};
    foreach my $i (0..@array-1) {
        return $i if $value == $array[$i];
    }
}

# This routine is really not necessary in light of the new
# `~~' operator in Perl.  Will use the new operator in the
# next version.
sub vector_equal {
    my $vec1 = shift;
    my $vec2 = shift;
    die "wrong data types for distance calculation\n" if @$vec1 != @$vec2;
    foreach my $i (0..@$vec1-1){
        return 0 if $vec1->[$i] != $vec2->[$i];
    }
    return 1;
}

# Returns the minimum value and its positional index in an array
sub minimum {
    my $arr = shift;
    my $min;
    my $index;
    foreach my $i (0..@{$arr}-1) {
        if ( (!defined $min) || ($arr->[$i] < $min) ) {
            $index = $i;
            $min = $arr->[$i];
        }
    }
    return ($min, $index);
}

sub minmax {
    my $arr = shift;
    my $min;
    my $max;
    foreach my $i (0..@{$arr}-1) {
        if ( (!defined $min) && (!defined $max) ) {
            $min = $arr->[$i];
            $max = $arr->[$i];
        } elsif ( $arr->[$i] < $min ) {
            $min = $arr->[$i];
        } elsif ( $arr->[$i] > $max ) {
            $max = $arr->[$i];
        }
    }
    return ($min, $max);
}

# Meant only for constructing a deep copy of an array of arrays
sub deep_copy {
    my $ref_in = shift;
    my $ref_out;
    foreach my $i (0..@{$ref_in}-1) {
        foreach my $j (0..@{$ref_in->[$i]}-1) {
            $ref_out->[$i]->[$j] = $ref_in->[$i]->[$j];
        }
    }
    return $ref_out;
}

sub display_cluster_centers {
    my @cluster_center_arr = @{shift @_};
    my $i = 1;
    foreach my $ele (@cluster_center_arr) {
        print "Cluster center $i: " .
               "@{[ map {my $x = sprintf('%.4f', $_); $x} @$ele ]}\n";
        $i++;
    }
}

# For displaying the individual clusters on a terminal
# screen.  Each cluster is displayed through the symbolic
# names associated with the data points.
sub display_clusters {
    my @clusters = @{shift @_};
    my $i = 1;
    foreach my $cluster (@clusters) {
        @$cluster = sort @$cluster;
        my $cluster_size = @$cluster;
        print "\n\nCluster $i ($cluster_size records):\n";
        foreach my $ele (@$cluster) {
            print "  $ele";
        }
        $i++
    }
    print "\n\n";
}

# from perl docs:
sub fisher_yates_shuffle {                
    my $arr =  shift;                
    my $i = @$arr;                   
    while (--$i) {                   
        my $j = int rand( $i + 1 );  
        @$arr[$i, $j] = @$arr[$j, $i]; 
    }
}

1;

__END__

=head1 NAME

Algorithm::KMeans - Clustering multi-dimensional data with a pure-Perl implementation

=head1 SYNOPSIS

  use Algorithm::KMeans;

  #  Name the data file:

  my $datafile = "mydatafile.dat";


  #  Set the mask to tell system which columns of the datafile to use for 
  #  clustering and which column contains a symbolic ID for each data record.
  #  For example, if the symbolic name is in the first column, if you want the
  #  second column to be ignored, and you want the next three columns to be 
  #  used for 3D clustering:

  my $mask = "N0111";


  #  Now construct an instance of the clusterer.  The parameter K controls 
  #  the number of clusters.  If you know how many clusters you want (in 
  #  this case 3), call

  my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                          mask     => $mask,
                                          K        => 3,
                                          terminal_output => 1,
                                        );


  #  Set K to 0 if you want the module to figure out the optimum number of 
  #  clusters from the data (it is best to run this option with the 
  #  terminal_output set to 1 so that you can see the different value of 
  #  QoC for the different K):

  my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                          mask     => $mask,
                                          K        => 0,
                                          terminal_output => 1,
                                        );

  #  For very large data files, setting K to 0 will result in searching
  #  through too many values for K.  For such cases, you can range limit
  #  the values of K to search through by

  my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                          mask     => "N111",
                                          Kmin     => 3,
                                          Kmax     => 10,
                                          terminal_output => 1,
                                        );

  #  Use the following call if you wish for the clusters to be written out 
  #  to files. Each cluster will be deposited in a file named 'ClusterX.dat' 
  #  with X starting from 0:

  my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                          mask     => $mask,
                                          K        => $K,
                                          write_clusters_to_files => 1,
                                        );


  #  FOR ALL CASES ABOVE, YOU'D NEED TO MAKE THE FOLLOWING CALLS ON THE 
  #  CLUSTERER INSTANCE TO ACTUALLY CLUSTER THE DATA:

  $clusterer->read_data_from_file();
  $clusterer->kmeans();


  #  If you want to directly access the clusters and the cluster centers in 
  #  your top-level script:

  my ($clusters, $cluster_centers) = $clusterer->kmeans();

  #  You can now access the symbolic data names in the clusters directly, 
  #  as in:

  foreach my $cluster (@$clusters) {
      print "Cluster:   @$cluster\n\n"
  }


  # CLUSTER VISUALIZATION:

  #  You must first set the mask for cluster visualization. This mask tells 
  #  the module which 2D or 3D subspace of the original data space you wish 
  #  to visualize the clusters in:

  my $visualization_mask = "111";
  $clusterer->visualize($visualization_mask);


  # SYNTHETIC DATA GENERATION:

  #  The module has been provided with a class method for generating 
  #  multivariate data for experimenting with clustering.  The data 
  #  generation is controlled by the contents of the parameter file that 
  #  is supplied as an argument to the data generator method.  The mean 
  #  and covariance matrix entries in the parameter file must be according 
  #  to the syntax shown in the param.txt file in the examples directory.  
  #  It is best to edit this file as needed:

  my $parameter_file = "param.txt";
  my $out_datafile = "mydatafile.dat";
  Algorithm::KMeans->cluster_data_generator(
                          input_parameter_file => $parameter_file,
                          output_datafile => $out_datafile,
                          number_data_points_per_cluster => $N );

=head1 CHANGES

Version 1.1.1 allows for range limiting the values of K to
search through.  K stands for the number of clusters to
form.  This version also declares the module dependencies in
the Makefile.PL file.

Version 1.1 is a an object-oriented version of the
implementation presented in version 1.0.  The current
version should lend itself more easily to code extension.
You could, for example, create your own class by subclassing
from the class presented here and, in your subclass, use
your own criteria for the similarity distance between the
data points and for the QoC (Quality of Clustering) metric,
and, possibly a different rule to stop the iterations.
Version 1.1 also allows you to directly access the clusters
formed and the cluster centers in your calling script.

=head1 DESCRIPTION

B<Algorithm::KMeans> is a I<perl5> module for the clustering
of numerical data in multidimensional spaces.  Since the
module is entirely in Perl (in the sense that it is not a
Perl wrapper around a C library that actually does the
clustering), the code in the module can easily be modified
to experiment with several aspects of automatic clustering.
For example, one can change the criterion used to measure
the "distance" between two data points, the stopping
condition for accepting final clusters, the criterion used
for measuring the quality of the clustering achieved, etc.

A K-Means clusterer is a poor man's implementation of the EM
algorithm.  EM stands for Expectation Maximization. For the
case of isotropic Gaussian data, the results obtained with a
good K-Means implementation should match those obtained with
the EM algorithm. Clustering with K-Means takes place
iteratively and involves two steps: 1) assignment of data
samples to clusters; and 2) Recalculation of the cluster
centers.  The assignment step can be shown to be akin to the
Expectation step of the EM algorithm, and the calculation of
the cluster centers akin to the Maximization step of the EM
algorithm.

Of the two key steps of the K-Means algorithm, the
assignment step consists of assigning each data point to
that cluster from whose center the data point is the
closest.  That is, during assignment, you compute the
distance between the data point and each of the current
cluster centers.  You assign the data sample on the basis of
the minimum value of the computed distance.  The second step
consists of re-computing the cluster centers for the newly
modified clusters.

Obviously, before the two-step approach can proceed, we need
to initialize the both the cluster center values and the
clusters that can then be iteratively modified by the
two-step algorithm.  How this initialization is carried out
is very important.  The implementation here uses a random
number generator to find K random integers between 0 and N
where N is the total number of data samples that need to be
clustered and K the number of clusters you wish to form.
The K random integers are used as indices for the data
samples in the overall data array --- the data samples thus
selected are treated as seed cluster centers.  This
obviously requires a prior knowledge of K.

How to specify K is one of the most vexing issues in any
approach to clustering.  In some case, we can set K on the
basis of prior knowledge.  But, more often than not, no such
prior knowledge is available.  When the programmer does not
explicitly specify a value for K, the approach taken in the
current implementation is to try all possible values between
2 and some largest possible value that makes statistical
sense.  We then choose that value for K which yields the
best value for the QoC (Quality of Clustering) metric.  It
is generally believed that the largest value for K should
not exceed sqrt(N/2) where N is the number of data point to
be clustered.

How to set the QoC metric is obviously a critical issue unto
itself.  In the current implementation, the value of QoC is
a ratio of the average radius of the clusters and the
average distance between the cluster centers.  But note that
this is a good criterion only when the data exhibits the
same variance in all directions.  When the data variance is
different directions, but still remains the same for all
clusters, a more appropriate QoC can be formulated using
other distance metrics such as the Mahalanobis distance.

Every iterative algorithm requires a stopping criterion.
The criterion implemented here is that we stop iterations
when there is no re-assignment of the data points during the
assignment step.

Ordinarily, the output produced by a K-Means clusterer will
correspond to a local minimum for the QoC values, as opposed
to a global minimum.  The current implementation protects
against that, but only in a very small way, by trying
different randomly selected initial cluster centers and then
selecting the one that gives the best overall QoC value.

=head1 METHODS

The module provides the following methods for clustering,
for cluster visualization, and for the generation of data
for testing a clustering algorithm:

=over

=item new()

    my $clusterer = Algorithm::KMeans->new(datafile => $datafile,
                                           mask     => $mask,
                                           K        => $K,
                                           terminal_output => 1,     
                                           write_clusters_to_files => 1,
                                          );

A call to new() constructs a new instance of the
Algorithm::KMeans class.  When $K is a non-zero positive
integer, the module will construct exactly that many
clusters.  However, when $K is 0, the module will find the
best number of clusters to partition the data into.

The data file is expected to contain entries in the
following format

   c20  0  10.7087017086940  9.63528386251712  10.9512155258108  ...
   c7   0  12.8025925026787  10.6126270065785  10.5228482095349  ...
   b9   0  7.60118206283120  5.05889245193079  5.82841781759102  ...
   ....
   ....

where the first column contains the symbolic ID tag for each
data record and the rest of the columns the numerical
information.  As to which columns are actually used for
clustering is decided by the string value of the mask.  For
example, if we wanted to cluster on the basis of the entries
in just the 3rd, the 4th, and the 5th columns above, the
mask value would be `N0111' where the character `N'
indicates that the ID tag is in the first column, the
character '0' that the second column is to be ignored, and
the '1's that follow that the 3rd, the 4th, and the 5th
columns are to be used for clustering.

The parameter `terminal_output' is boolean; when not
supplied in the call to new() it defaults to 0.  When set,
this parameter determines what you will see on the terminal
screen of the window in which you make these method calls.
When set to 1, you will see on the terminal screen the
different clusters as lists of the symbolic IDs and their
cluster centers. You will also see the QoC (Quality of
Clustering) value for the clusters displayed. 

The parameter `write_clusters_to_files' is boolean; when not
supplied in the call to new(), it defaults to 0.  When set
to 1, the clusters are written out to files named
    
     Cluster0.dat 
     Cluster1.dat 
     Cluster2.dat
     ...
     ...

Before the clusters are written to these files, the module
destroys all files with such names in the directory in which
you call the module.

If you wish for the clusterer to search through a
(Kmin,Kmax) range of values for K, the constructor should be
called in the following fashion:

    my $clusterer = Algorithm::KMeans->new(datafile => $datafile,
                                           mask     => $mask,
                                           Kmin     => 3,
                                           Kmax     => 10,
                                           terminal_output => 1,     
                                          );

where obviously you can choose any reasonable values for
Kmin and Kmax.  If you choose a value for Kmax that is
statistically too large, the module will let you know.

=item read_data_from_file()

    $clusterer->read_data_from_file()

=item kmeans()

    $clusterer->kmeans();

    or 

    my ($clusters, $cluster_centers) = $clusterer->kmeans();

The first call above works solely by side-effect.  The
second call also returns the clusters and the cluster
centers.

=item get_K_best()

    $clusterer->get_K_best();

This call makes sense only if you supply either the K=0
option to the constructor, or you specify values for the
Kmin and Kmax options. The K=0 and the (Kmin,Kmax) options
cause the KMeans algorithm to figure out on its own the best
value for K.  Remember, K is the number of clusters the data
is partitioned into.

=item show_QoC_values()

    $clusterer->show_QoC_values();

presents a table with K values in the left column and the
corresponding QoC (Quality-of-Clustering) values in the
right column.  Note that this call makes sense only if you
either supply the K=0 option to the constructor, or you
specify values for the Kmin and Kmax options.

=item visualize()

    $clusterer->visualize( $visualization_mask )

The visualization mask here does not have to be identical to
the one used for clustering, but must be a subset of that
mask.  This is convenient for visualizing the clusters in
two- or three-dimensional subspaces of the original space.

=item  cluster_data_generator()

    Algorithm::KMeans->cluster_data_generator(
                            input_parameter_file => $parameter_file,
                            output_datafile => $out_datafile,
                            number_data_points_per_cluster => 20 );

for generating multivariate data for clustering if you wish
to play with synthetic data for clustering.  The input
parameter file contains the means and the variances for the
different Gaussians you wish to use for the synthetic data.
See the file param.txt provided in the examples directory.
It will be easiest for you to just edit this file for your
data generation needs.  In addition to the format of the
parameter file, the main constraint you need to observe in
specifying the parameters is that the dimensionality of the
covariance matrix must correspond to the dimensionality of
the mean vectors.  The multivariate random numbers are
generated by calling the Math::Random module.  As you would
expect, this module requires that the covariance matrices
you specify in your parameter file be symmetric and positive
definite.  Should the covariances in your parameter file not
obey this condition, the Math::Random module will let you
know.

=back

=head1 HOW ARE THE CLUSTERS OUTPUT?

When the option terminal_output is set in the call to the
constructor, the clusters are displayed on the terminal
screen.

When the option write_clusters_to_files is set in the call
to the constructor, the module dumps the clusters in files
named

    Cluster0.dat
    Cluster1.dat
    Cluster2.dat
    ...
    ...

in the directory in which you execute the module.  The
number of such files will equal the number of clusters
formed.  All such existing files in the directory are
destroyed before fresh ones created.  Each cluster file
contains the symbolic ID tags of the data points in that
cluster.

=head1 REQUIRED

This module requires the following two modules:

   Math::Random
   Graphics::GnuplotIF

the former for generating the multivariate random numbers
and the latter for the visualization of the clusters.

=head1 EXAMPLES

See the examples directory in the distribution for how to
make calls to the clustering and the visualization methods.
The examples directory also includes a parameter file,
param.txt, for generating synthetic data for clustering.
Just edit this file if you would like to generate your own
multivariate data for clustering.

=head1 EXPORT

None by design.

=head1 CAVEATS

Please note that this clustering module is not meant for
very large datafiles.  Being an all-Perl implementation, the
goal here is not the speed of execution.  On the contrary,
the goal is to make it easy to experiment with the different
facets of K-Means clustering.  If you need to process a
large data file, you'd be better off with a module like
Algorithm::Cluster.  However note that when you use a
wrapper module in which it is a C library that is actually
doing the job of clustering for you, it is more difficult to
experiment with the various aspects of clustering.  At the
least, you have to recompile the code for every change you
make to the source code of a low-level library.  You are
spared that frustration with an all-Perl implementation.

Clustering usually does not work well when the data is
highly anisotropic, that is, when the data has very
different variances along its different dimensions.  This
problem becomes particularly severe when the different
clusters you expect to see in the data have non-uniform
anisotropies.  When the anisotropies are uniform, one could
improve on the implementation shown in the current module by
first normalizing the data with respect to the different
variances along the different dimensions.  One could also
try to cluster the data in a low-dimensional space formed by
a principal components analysis of the data, especially
after the variances are normalized out.  Depending on how
the current module is received, its future versions may
include those enhancements.


=head1 BUGS

Please send me email to the author if you find any bugs.
They will be dealt with promptly. When sending email, please
place the string 'KMeans' in the subject line.

=head1 INSTALLATION

The usual

    perl Makefile.PL
    make
    make test
    make install

if you have root access.  If not, 

    perl Makefile.PL prefix=/some/other/directory/
    make
    make test
    make install

=head1 THANKS

Chad Aeschliman was kind enough to test out the interface of
this module and to give suggestions for its improvement.  His
key slogan: "If you cannot figure out how to use a module in
under 10 minutes, it's not going to be used."  That should
explain the longish Synopsis included here.

=head1 AUTHOR

Avinash Kak, kak@purdue.edu

If you send email, please place the string "KMeans" in your
subject line to get past my spam filter.

=head1 COPYRIGHT

This library is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

 Copyright 2010 Avinash Kak

=cut

