use ExtUtils::MakeMaker;
use File::Copy;

if ($^V lt v5.10) {

   die("Algorithm::KMeans has only been tested on Perl 5.10.0.\n" .
   "Your perl version is $].\n");
}

copy("perl/MANIFEST.perl","MANIFEST");

WriteMakefile(
	NAME         => 'Algorithm::KMeans',
	VERSION_FROM => 'KMeans.pm',
        AUTHOR       => 'Avinash Kak (kak@purdue.edu)',
        ABSTRACT     => 'A pure-Perl implementation of K-Means Clustering',
);
