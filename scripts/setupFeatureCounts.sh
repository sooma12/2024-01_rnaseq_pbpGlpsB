# setupFeatureCounts.sh

# Install subread package
# Downloaded subread-2.0.6-source.tar.gz from https://sourceforge.net/projects/subread/
# scp to /work/geisingerlab/Mark/software/
cd /work/geisingerlab/Mark/software/
tar zxvf subread-2.0.6-source.tar.gz
rm zxvf subread-2.0.6-source.tar.gz
cd subread-2.0.6-source/src
make -f Makefile.Linux

# Set up featureCounts as a module
# Binary files including featureCounts are in Mark/software/subread-2.0.6-source/bin
# Note, I had previously added module use --append /usr/local/usrapps/bioinfo/modulefiles to .bashrc
cd /work/geisingerlab/Mark/software/modulefiles
mkdir -p subread
cd subread
nano 2.0.6

#cat 2.0.6
##%Module
#proc ModulesHelp { } {
#   puts stderr "This module adds gffread to your path"
#}
#
#module-whatis "This module adds gffread to your path\n"
#
#set root /work/geisingerlab/Mark/software/subread-2.0.6-source/bin
#prepend-path PATH $root
#prepend-path LD_LIBRARY_PATH $root

