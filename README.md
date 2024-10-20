
# ORGANOFEM, FEM applied to organoids biomechanical modelling  

<img src="https://user-images.githubusercontent.com/56252845/160375009-8c8537f5-6f84-4a2e-a21d-f0994fafde4b.png" width="300" height="300" />

# organofem_published

This repository contain the codes used in the following publication:

**Deciphering the interplay between biology and physics with a finite element method-implemented vertex organoid model: a tool for the mechanical analysis of cell behavior on a spherical organoid shell.

J. Laussu1*, D. Michel2, L. Magne1,2, S. Segonds1, S. Marguet1, D. Hamel2, M. Quaranta-Nicaise², F. Barreau2, E. Mas2,3, V. Velay4, F. Bugarin1* and A. Ferrand2*  
1/  Institut Clément Ader, Université Fédérale de Toulouse Midi-Pyrénées, Institut Clément Ader – CNRS UMR 5312 – UPS/INSA/Mines Albi/ISAE, 3, rue Caroline Aigle, Toulouse 31400, France  
2/ IRSD, Université de Toulouse, INSERM, INRAE, ENVT, UPS, U1220, CHU Purpan, CS60039, 31024 Toulouse, France  
3/ Centre de référence des maladies rares digestives et service de Gastroentérologie, Hépatologie, Nutrition, Diabétologie et Maladies Héréditaires du Métabolisme  
4/ Institut Clément Ader (ICA), Université de Toulouse, CNRS, IMT Mines Albi, INSA, ISAE-SUPAERO, UPS, Campus Jarlard, F-81013, Albi, France.  

# architecture of the project  

    //This git project is constructed of: 
    
    README.md
    ├── docs 
    │   ├── TOC.md # table of content
    │   ├── faq.md # frequently asked question
    │   ├── misc.md # miscellaneous
    │   └── usage.md # frequently used commands
    ├── lib  
    │   ├── python
    │   |   └── lib.py # python libraries     
    │   └── other
    │       └── wrap.py # other languages libraries 
    ├── test  
    │   ├── benchmarks
    │   |   └──  # benchmarks tests     
    │   ├── integration  
    │   |   └──  # integration tests
    │   ├── temp  
    │   |   └──  # for temporary resulst
    │   ├── article_rawdata  
    │   |   └──  # rawdatas from the article  
    │   └── unit
    │       └──  # unit tests 
    ├── dist 
    │   ├── io
    |   |   └──  # in/ out functions 
    │   ├── utils
    |   |   └──  # common operations and functions 
    |   ├── model
    │   |   ├── seg_organoid.tif # segmentation image to test the code step 1 
    │   |   ├── analyze_imseg_organofem.py # STEP1 python code to extract morphological data from segmented image
    │   |   ├── create_AVM_organofem.py # STEP2 python code to create a virtual organoid in an activ vertex model
    │   |   ├── node_inp_organofem.inp # input file for Abaqus that can be use in step 3 to describe organoid elemental structure
    │   |   ├── node_fil_organofem.fil # output Abaqus file with historic of nodes field information / for step 3  
    │   |   ├── elem_fil_organofem.fil # output Abaqus file with historic of elements field information / for step 3   
    |   |   └── analyze_FEM_organofem.py # STEP1 python code to extract morphological data from FEM model resolution
    |   ├── bio
    │   |   ├── mesh.py # mesh class
    │   |   ├── organoid.py # element class  
    │   |   ├── cells.py # edge class  
    │   |   ├── element.py # element class  
    │   |   ├── edge.py # edge class  
    │   |   └── node.py # node class  
    │   ├── bio  
    │   │   ├── tissu
    │   │   |   ├── behaviors.py # tissu behaviors class  
    │   │   |   └── properties.py # tissu properties class 
    │   │   └── cell  
    │   │       ├── behaviors.py # cell behaviors class
    │   │       └── properties.py # cell properties class 
    │   └── meca
    │       └── properties.py # tissu properties class 
    VERSION # Git version
    LICENCE # GNU licence v3
 

The main functions are the following:

________________________________________
## 1- analyse_imseg_organofem.py
Code developed to extract morphological cellular information like cellular surface area from segmented images. 

![step1](https://github.com/user-attachments/assets/de2eecd8-f92f-4ec1-b06d-5b346457afa5)

## 2- create_AVM_organofem.py##
Active Vertex Model for the simulation of a virtual organoid given the number of cells and the thickness of the oraganoid shell.
This part have an optional quasi static solver resolution of the oraganoid morphology with parameters of elasticity depending of a target tensions and volumes to achieve.
An other property of this code is the translation of the vertex model in a finite elment model formulation with the creation an .inp file for Abaqus.

![step2](https://github.com/user-attachments/assets/c7333a04-9eee-4801-a99c-bb09f4765b17)

##  3- analyse_FEM_organofem.py
This part of the code is an example of how to read Abaqus output files (.fil) to extract properties like strain and stress and the more important how to acces to the position/displacement of the node.
By this way, we can measure morphological properties at the cellular level after deformation.

![step3](https://github.com/user-attachments/assets/1fe37d81-f8fa-4199-9c28-ed17cfeb21dd)
________________________________________
## Installation guide  
intall anaconda  
create conda environment:  
```console
conda create -n organofem python=3.9
```
activate the environment:  
```console
conda activate organofem
```
install tyssue library using conda-forge depository:  
```console
conda install -c conda-forge tyssue
```
install libtiff library using conda-forge depository:  
```console
conda install -c conda-forge pylibtiff
```
install additionnal libraries:  
```console
pip install vedo hickle scikit-learn
```
enter in the folder containing the examples:  
```console  
cd organofem\dist\model
```
run an example:  
```console  
python "python code from example 1, 2 or 3"
```  


