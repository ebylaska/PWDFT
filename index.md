![PWDFT Logo](https://raw.githubusercontent.com/ebylaska/PWDFT/gh-pages/assets/PWDFT-LOGO.png)
# PWDFT - Plane-Wave DFT development for [NWChemEx](https://www.exascaleproject.org/research-project/nwchemex/)


Based on [NWPW](https://nwchemgit.github.io/Plane-Wave-Density-Functional-Theory.html) module in [NWChem]√

https://pubs.acs.org/doi/10.1021/acs.chemrev.0c00998

We would like thank the DOE ECC program and the DOE OS OBER EMSL project for providing support that helped with the initial development of PWDFT.


For more information contact Eric Bylaska (eric.bylaska@pnnl.gov)

[PWDFT Documentation](https://ebylaska.github.io/PWDFT/Nwpw/Plane-Wave-Density-Functional-Theory.html)

Overview: The NWChemex PWDFT code uses pseudopotentials and plane-wave basis sets to perform Density Functional Theory calculations (simple introduction pw-lecture.pdf). This module complements the capabilities of the more traditional Gaussian function based approaches by having an accuracy at least as good for many applications, yet is still fast enough to treat systems containing hundreds of atoms. Another significant advantage is its ability to simulate dynamics on a ground state potential surface directly at run-time using the Car-Parrinello algorithm. This method’s efficiency and accuracy make it a desirable first principles method of simulation in the study of complex molecular, liquid, and solid state systems. Applications for this first principles method include the calculation of free energies, search for global minima, explicit simulation of solvated molecules, and simulations of complex vibrational modes that cannot be described within the harmonic approximation.

