base_url=$1/
dir=/cellar/users/snwright/Data/RareCommon/GWASCatalog/testing_traits/

cd $dir

# get initial index
wget -nH -r -np --cut-dirs=5 --span-hosts --convert-links $base_url

cd $(basename "${base_url}")

# extract file links
grep -oP '(?<=href=")[^"]*' index.html | grep -v index | grep -v "?" | grep -v "/$"> file_links.txt

# download files
wget -i file_links.txt
rm file_links.txt
rm index.html

# get harmonised files
mkdir harmonised
cd harmonised

# get index
wget -nH -r -np --cut-dirs=7 --span-hosts --convert-links ${base_url}harmonised/

# get file links
grep -oP '(?<=href=")[^"]*' index.html | grep -v index | grep -v "?" | grep -v "/$"> file_links.txt

# download files
wget -i file_links.txt

rm file_links.txt
rm index.html
