# CRAG
* Unfinished project. Version now is 0.0.1



Dependency
------------------

* R: version > 3.6.0, packages: "genio", "glmnet"
* Python: packages: "pandas-plink"


Running the tests
-----------------
* For testing, run elasnet_generation.R. For example, Rscript elasnet_generation.R test test.frq annotation_Whole_Blood.txt gamma.txt expression.txt true_evaluation.txt

  The first argument is plink file location, second is frequency file location, third is annotation file location, fourth is initialized gamma output file location, fifth is gene expression output file location, last one is expected evaluation results.

* Then run read_data_elasnet.py. For example, python read_data_elasnet.py --plink test --expression expression.txt --gamma gamma.txt --frq test.frq --annotation annotation_Whole_Blood.txt --outSNP SNP.txt --outAnno anno.txt --outPerfer pre_evaluation.txt

  The first argument is plink file location, second is gene expression file location, third is initialized gamma file location which should be the result in the first step, fourth is frequency file location, fifth is annotation file location, sixth is SNP's coefficient output file location, seventh is annotation's coefficient output file location, the last one is predicted evaluation results.

Notes
---------
This package is still on going.


Source Repository
-----------------
CRAG's source-code repository is hosted here on GitHub.


Authors
---------

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Yuzheng Dun (maintainer) | ydun2@jhu.edu | Visiting student, Department of Biostatistics  JHU |
| Wei Liu | wei.liu.vivian@yale.edu | PhD student, Department of Biostatistics Yale |
| Hongyu Zhao | hongyu.zhao@yale.edu | Professor, Department of Biostatistics  Yale |
<!--- --->
                             


Licensing
---------
CRAG



