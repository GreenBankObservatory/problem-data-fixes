cd /stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified
echo $1
rm -rf $1/SDFITS/
mkdir $1/SDFITS
sdfits $1 -noprompt
mv $1.raw.vegas $1/SDFITS/
mv $1.raw.dcr.fits $1/SDFITS/
