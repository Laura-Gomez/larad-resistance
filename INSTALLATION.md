# INSTALLATION INSTRUCTIONS

To install and run the software in this repository follow these instructions


## Docker and nextflow installation

Please sure you have docker and nextflow installed in your system



## Docker repositories 

Pull the following docker repositories:

```
docker pull laugoro/bacterial-st:public
docker pull laugoro/workshop-inmegen-resistance:public
docker pull laugoro/workshop-inmegen-assembly:public
```


## Databases

Be sure to download and install the following databases:

Resfinder databases

```
git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/
```

Virulencefinder database

```
git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db/

```

RGI database

```
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
```

RESFAM database

```
wget http://dantaslab.wustl.edu/resfams/Resfams-full.hmm.gz
gunzip Resfams-full.hmm.gz
hmmpress Resfams-full.hmm
```
