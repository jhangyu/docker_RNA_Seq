from jhenv
run apt-get update
run conda config --add channels conda-forge
run conda config --add channels anaconda
run conda config --add channels bioconda
run conda install -y gatk4 hisat2 trimmomatic cufflinks sambamba stringtie samtools ensembl-vep stacks
run apt-get install perlbrew
run perlbrew init
run perlbrew install 5.27.5
run perlbrew switch perl-5.27.5
perlbrew install-cpanm
cpanm -v DBI
cpanm -v DBD::mysql
cpanm -v Bio::Perl 
conda install ensembl-vep
