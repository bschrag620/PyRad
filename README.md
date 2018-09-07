# PyRad
A python package for modeling radiative transfer through an atmsophere.

All absorption lines, intensities, and relevant info are downloaded from hitran.org and stored locally.
For a more comprehensive (ie professional package), be sure to check out HAPI (Hitran API) http://hitran.org/hapi. It is free and open source as well. This project has been more personal opportunity to test and improve my own understanding of radiative transfer through gases. The final version of this project should support modeling radiative transfer through any atmospheric composition. Currently, it serves as a gas cell simulator, similar to http://www.spectralcalc.com/calc/spectralcalc.php 

Version 1.
Currently, pyrad supports calculating radiative transmittance through a single layer of gas. The main file has an example of usage. 

#########
Thanks to contributions and guidance from the following:

HAPI Interface - R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016)

Tom Marshall of GATS-Inc.com, for providing guidance in resolving personal ignorances of units.

Eli Rabett, for his patience in helping me understand concepts that I have not received an education on previously.
