#!/bin/bash
# Set work directory with argument. 
if [ $# -ne 1 ]; then
  echo "Please specify the working directory." 1>&2
  echo "EXAMPLE:" 1>&2
  echo "    ./setup.sh souyakuchan_library" 1>&2
  exit 1
fi
WORK_DIR=$1

# Create database.
echo 'Create database...'
docker-compose exec db createdb -U postgres souyakuRB2018
docker-compose exec db psql -U postgres -c 'create extension rdkit' souyakuRB2018
# docker-compose exec web python setup_for_souyakuRB2018/setup.py

# Create raw_data table and register compounds into it.
echo 'Create raw_data table and register compounds into it...'
docker-compose exec db psql -U postgres -c "create table raw_data (id serial, smiles text, name text, library_id integer)" souyakuRB2018
cat $WORK_DIR/Enamine_Advanced_collection_desalt.smi | sed -e 's/	/ /g' | \
    docker-compose exec -T db psql -U postgres -c "copy raw_data (smiles, name) from stdin with delimiter ' '" souyakuRB2018
docker-compose exec -T db psql -U postgres -c "update raw_data set library_id = 1  where library_id is null" souyakuRB2018

cat $WORK_DIR/Enamine_Premium_collection_desalt.smi | sed -e 's/	/ /g' | \
    docker-compose exec -T db psql -U postgres -c "copy raw_data (smiles, name) from stdin with delimiter ' '" souyakuRB2018
docker-compose exec -T db psql -U postgres -c "update raw_data set library_id = 2  where library_id is null" souyakuRB2018

cat $WORK_DIR/Enamine_HTS_collection_desalt.smi | sed -e 's/	/ /g' | \
    docker-compose exec -T db psql -U postgres -c "copy raw_data (smiles, name) from stdin with delimiter ' '" souyakuRB2018
docker-compose exec -T db psql -U postgres -c "update raw_data set library_id = 3 where library_id is null" souyakuRB2018

cat $WORK_DIR/UOS_HTS_desalt.smi | sed -e 's/	/ /g' | \
    docker-compose exec -T db psql -U postgres -c "copy raw_data (smiles, name) from stdin with delimiter ' '" souyakuRB2018
docker-compose exec -T db psql -U postgres -c "update raw_data set library_id = 4 where library_id is null" souyakuRB2018

# Create libraries table.
echo 'Create libraries table...'
docker-compose exec -T db psql -U postgres -c "create table libraries (id serial primary key, name text)" souyakuRB2018
docker-compose exec -T db psql -U postgres -c "insert into libraries (name) values ('Enamine_Advanced_collection'), ('Enamine_Premium_collection'), ('Enamine_HTS_collection'), ('UOS_HTS')" souyakuRB2018

# Convert smiles string to mol and copy them to compounds table.
# It takes about an hour on my pc.
echo 'Convert smiles string to mol and copy them to compounds table...'
docker-compose exec -T db psql -U postgres -c "select * into compounds from (select id, mol_from_smiles (smiles::cstring) structure, name, library_id, morganbv_fp(mol_from_smiles(smiles::cstring), 2) mfp2, featmorganbv_fp(mol_from_smiles(smiles::cstring), 2) ffp2, atompairbv_fp(mol_from_smiles(smiles::cstring)) atompair, torsionbv_fp(mol_from_smiles(smiles::cstring)) torsion from raw_data) tmp where structure is not null" souyakuRB2018

# Assign id to primary key.
echo 'Assign id to primary key...'
docker-compose exec -T db psql -U postgres -c "alter table compounds add primary key (id)" souyakuRB2018

# Assign library_id to foreign key.
echo 'Assign library_id to foreign key...'
docker-compose exec -T db psql -U postgres -c "alter table compounds add foreign key (library_id) references libraries(id)" souyakuRB2018

# Create index to structure column.
# It takes about an hour on my pc.
echo 'Create index to structure column...'
docker-compose exec -T db psql -U postgres -c "create index structure_idx on compounds using gist(structure)" souyakuRB2018

# Create index to fingerprint column.
# It takes about an hour on my pc.
echo 'Create index to mfp2 column...'
docker-compose exec -T db psql -U postgres -c "create index mfp2_idx on compounds using gist(mfp2)" souyakuRB2018
echo 'Create index to ffp2 column...'
docker-compose exec -T db psql -U postgres -c "create index ffp2_idx on compounds using gist(ffp2)" souyakuRB2018
echo 'Create index to atompair column...'
docker-compose exec -T db psql -U postgres -c "create index atompair_idx on compounds using gist(atompair)" souyakuRB2018
echo 'Create index to torsion column...'
docker-compose exec -T db psql -U postgres -c "create index torsion_idx on compounds using gist(torsion)" souyakuRB2018

# Drop raw_data table.
echo 'Drop raw_data table...'
docker-compose exec -T db psql -U postgres -c "drop table raw_data" souyakuRB2018

echo 'Done.'
