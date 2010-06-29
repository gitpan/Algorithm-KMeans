use Test::Simple tests => 3;

use lib '../blib/lib','../blib/arch';

use Algorithm::KMeans;

# Test 1 (Data Generation):

my $datafile = "__testdata.dat";
Algorithm::KMeans->cluster_data_generator( 
                        output_datafile => $datafile,
                        number_data_points_per_cluster => 20 );
open IN, $datafile;
my @data_records = <IN>;
ok( @data_records == 60,  'Data generation works' );


# Test 2 (K-Means Clustering):

my $mask = "N111";
my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                        mask     => "N111",
                                        K        => 3,
    );
$clusterer->read_data_from_file();
my ($clusters, $cluster_centers) = $clusterer->kmeans();
ok( @$clusters == 3,  'Clustering works' );

# Test 3 (Data Visualization)

eval {
    my $visualization_mask = "111";
    my $pause_time = 2;
    $clusterer->visualize($visualization_mask, $pause_time);
};
print ${$@} if ($@); 

ok( !$@,  'Visualization works' );

unlink "__testdata.dat";
