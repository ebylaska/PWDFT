![PWDFT Logo](https://raw.githubusercontent.com/ebylaska/PWDFT/gh-pages/assets/PWDFT-LOGO.png)
# PWDFT - PW-DFT development for [NWChemEx](https://www.exascaleproject.org/research-project/nwchemex/)


Based on [NWPW](https://nwchemgit.github.io/Plane-Wave-Density-Functional-Theory.html) module in [NWChem]√

https://pubs.acs.org/doi/10.1021/acs.chemrev.0c00998

We would like thank the DOE ECC program and the DOE OS OBER EMSL project for providing support that helped with the initial development of PWDFT.

Are you just learning NWChem and would like to have an easy way to generate input decks, check your output decks against a large database of calculations, perform simple thermochemistry calculations, calculate the NMR and IR spectra of modest size molecule, or just try out NWChem before installing it? EMSL Arrows scientific service can help. A Web API to EMSL Arrows is now available for Alpha testing. Click on this link.

For more information contact Eric Bylaska (eric.bylaska@pnnl.gov)

The NWChem molecular modeling software implements a robust and diverse set of molecular theories that can estimate the thermodynamics and kinetics of molecules and materials. It arguably has the most capabilities of any molecular modeling code today. The problem with NWChem and other molecular modeling codes is that:

Molecular modeling software is extremely complex, contains millions of lines of code, and takes a long time to set up and to learn how to use.
Even the most basic input for molecular modeling software requires the use of other software to generate it.
Because of this complexity people unnaturally identify with codes and molecular theories, and they are hesitant to learn new codes and new molecular simulation techniques.
TinyArrows is a software package that combines NWChem, SQL and NOSQL databases, and web applications that simplifies molecular and materials modeling and makes these modeling capabilities easier to use and more accessible to many scientists and engineers and students. TinyArrows is very simple to use. The user just enters chemical reactions into one, of serveral available web applications, and then results are posted back with thermodynamic, reaction pathway (kinetic), spectroscopy, and other results.

[PWDFT Documentation](https://ebylaska.github.io/PWDFT/Nwpw/Plane-Wave-Density-Functional-Theory.md)


Overview
The NWChemex PWDFT code uses pseudopotentials and plane-wave basis sets to perform Density Functional Theory calculations (simple introduction pw-lecture.pdf). This module complements the capabilities of the more traditional Gaussian function based approaches by having an accuracy at least as good for many applications, yet is still fast enough to treat systems containing hundreds of atoms. Another significant advantage is its ability to simulate dynamics on a ground state potential surface directly at run-time using the Car-Parrinello algorithm. This method’s efficiency and accuracy make it a desirable first principles method of simulation in the study of complex molecular, liquid, and solid state systems. Applications for this first principles method include the calculation of free energies, search for global minima, explicit simulation of solvated molecules, and simulations of complex vibrational modes that cannot be described within the harmonic approximation.

