#1. samtools 

#2. bedtools

#3. IGV 
wget http://data.broadinstitute.org/igv/projects/downloads/IGV_2.3.89.zip
unzip IGV_2.3.89.zip
#java -jar IGV_2.3.89/igv.jar 

# 4. python scipy: https://www.scipy.org/install.html 
sudo apt-get install python-numpy python-scipy 

# 5. tabix 
git clone https://github.com/samtools/htslib
cd htslib;
make

#6. hapcut
git clone https://github.com/vibansal/hapcut.git
cd hapcut;
make all
