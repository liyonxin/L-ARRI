DATABASE_DIR="database"
mkdir -p "$DATABASE_DIR"

######### centrifuge database #################################
echo "
Downloading Centrifuge database"
wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz --output-document 'p+b+v.tar.gz'
tar -zvxf p+b+v.tar.gz
mv centrifuge-1.0.3-beta database
wget https://github.com/liyonxin/L-ARRI/blob/main/new-species-all.txt
mv new-species-all.txt database

######### MobileOG-db #################################
echo "
Downloading MGE database"
wget "https://mobileogdb-downloads.s3.us-east-2.amazonaws.com/data-version-files/beatrix-1-6_v1_all.zip" 
unzip beatrix-1-6_v1_all.zip
mv mobileOG-db_beatrix-1.6.All.faa database

######### sARG-db #################################
echo "
Downloading sARG database"
wget https://github.com/liyonxin/L-ARRI/blob/main/SARG_20211207_14210_filter.ffn
wget https://github.com/liyonxin/L-ARRI/blob/main/structure_20181107.LIST
mv SARG_20211207_14210_filter.ffn database
mv structure_20181107.LIST database

