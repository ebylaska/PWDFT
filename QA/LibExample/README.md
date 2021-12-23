## Generate pwdft shared library ##
First go to the PWDFT home directory, e.g.,

https://jsfiddle.net/8ndx694g/


<img src="https://render.githubusercontent.com/render/math?math=\large E = E_{QM} %2B E_{MM} %2B E_{QM/MM}">

<img src="https://render.githubusercontent.com/render/math?math=\large %5Cbegin%7Balign*%7D%0AE_%7BMM%7D%20%26%3D%20%5Cfrac%7B1%7D%7B2%7D%20%5Csum_%7Bi%20%5Cneq%20j%7D%20%5Cfrac%7Bq_i%20q_j%7D%7Br_%7Bij%7D%7D%20%5C%5C%0A%20%20%20%20%20%20%20%26%20%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Csum_%7Bi%20%5Cneq%20j%7D%20%5C%7B%20%5Cfrac%7BA_%7Bij%7D%7D%7Br_%7Bij%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7Bij%7D%7D%7Br_%7Bij%7D%5E6%7D%20%5C%7D%0A%5Cend%7Balign*%7D">



<img src="https://render.githubusercontent.com/render/math?math=\large %5Cbegin%7Balign*%7D%0AE_%7BQM%2FMM%7D%20%26%3D%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BQ_I%20q_i%7D%7Br_%7BIi%7D%7D%20%5C%7D%5C%5C%0A%26%2B%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BA_%7BIi%7D%7D%7Br_%7BIi%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7BIi%7D%7D%7Br_%7BIi%7D%5E6%7D%20%5C%7D%0A%5Cend%7Balign*%7D">

<img src="https://render.githubusercontent.com/render/math?math=\large %5Cbegin%7Balign*%7D%0AE_%7BQM%2FMM%7D%20%26%3D%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BQ_I%20q_i%7D%7Br_%7BIi%7D%7D%20%5C%7D%5C%5C%0A%26%2B%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BA_%7BIi%7D%7D%7Br_%7BIi%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7BIi%7D%7D%7Br_%7BIi%7D%5E6%7D%20%5C%7D%5C%5C%0A%26%3D%20%5Csum_%7BI%3DQM%7D%20Q_I%20U_I%20%5C%5C%0A%26%2B%20e_%7BLJ%7D(I%2Ci)%0A%5Cend%7Balign*%7D">

```
cd /Users/bylaska/Codes/PWDFT
```

then excecute the following steps:

 1) mkdir build_library
 2) cd build_library
 3) cmake ../Nwpw -DMAKE_LIBRARY=true
 4) make


## Compiling on macOS ##

Go back to the QA directory, e.g.,

```
cd /Users/bylaska/Codes/PWDFT/QA/LibExample
```


Set location of the DYLD_LIBRARY_PATH

 1) setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared

Compile test.cpp using mpic++

 2) mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


## Compiling on LINUX ##
Add location of shared library location to LD_LIBRARY_PATH, e.g., 

1) setenv LD_LIBRARY_PATH setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib

Compile test.cpp using mpic++

 2) mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


## Running Example ##

To run the program just type 

```
a.out
```

or run using mpirun, e,g,

```
 mpirun -np 4 a.out
```

## Example Output ##

```
[Erics-MacBook-Pro:PWDFT/QA/LibExample] bylaska% mpirun -np 4 a.out
Hello world
np=4 taskid=0
dummy=Sat Dec 11 21:37:26 2021
rtdbstr={
    "current_task": "task pspw gradient",
    "dbname": "w2-pspw",
    "foundtask": false,
    "geometries": {
        "geometry": {
            "autosym": 0,
            "autoz": 0,
            "center": 1,
            "charges": [
                74.0,
                74.0
            ],
            "conv": 1.88972598858,
            "coords": [
                0.0,
                0.0,
                -2.173184886867,
                0.0,
                0.0,
                2.173184886867
            ],
            "masses": [
                183.951,
                183.951
            ],
            "nion": 2,
            "symbols": [
                "W",
                "W"
            ],
            "velocities": [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0
            ]
        }
    },
    "geometry": null,
    "nwinput_cur": 24,
    "nwinput_lines": [
        "Title \"Hellow Dolly\"",
        "",
        "memory 1900 mb",
        "start w2-pspw",
        "echo",
        "",
        "",
        "",
        "geometry noautosym noautoz center",
        "W 0 0 0    ",
        "W 0 0 2.3  ",
        "end",
        "",
        "nwpw",
        "   simulation_cell",
        "     SC 20.0",
        "   end",
        "   cutoff 20.0",
        "   xc beef",
        "end",
        "",
        "task pspw gradient",
        "",
        ""
    ],
    "nwinput_nlines": 24,
    "nwpw": {
        "cutoff": [
            20.0,
            40.0
        ],
        "simulation_cell": {
            "unita": [
                20.0,
                0.0,
                0.0,
                0.0,
                20.0,
                0.0,
                0.0,
                0.0,
                20.0
            ]
        },
        "xc": "beef",
	"apc":{"on":true}
    },
    "permanent_dir": ".",
    "psp_library_dir": "",
    "pspw": {
        "eigenvalues": [
            -0.1362203466912019,
            -0.16599015455140698,
            -0.18804690195660134,
            -0.2009022221514341,
            -0.2009114718004976,
            -0.22608901551944924
        ],
        "energies": [
            -19.386361266387933,
            -2.2363202253411827,
            12.056086535837377,
            -8.821534942920875,
            -1.740836619784234,
            13.2008277243924,
            -20.4908238457099,
            -13.590080118202735,
            24.112173071674754,
            -5.468417057495736,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
        ],
        "energy": -19.386361266387933,
        "fion": [
            0.007125932269872324,
            -0.02090412617808888,
            -0.040784071881918516,
            -0.007129528005492989,
            0.020918941207273663,
            0.04081718818417886
        ],
        "gradient": [
            -0.007125932269872324,
            0.02090412617808888,
            0.040784071881918516,
            0.007129528005492989,
            -0.020918941207273663,
            -0.04081718818417886
        ]
    },
    "scratch_dir": ".",
    "title": "Hellow Dolly"
}

          *****************************************************
          *                                                   *
          *               PWDFT PSPW Calculation              *
          *                                                   *
          *  [ (Grassmann/Stiefel manifold implementation) ]  *
          *  [              C++ implementation             ]  *
          *                                                   *
          *              version #7.00   02/27/21             *
          *                                                   *
          *    This code was developed by Eric J. Bylaska,    *
          *    Abhishek Bagusetty, David H. Bross, ...        *
          *                                                   *
          *****************************************************
          >>> job started at       Sat Dec 11 21:37:26 2021 <<<

 psp_library: /Users/bylaska/Codes/PWDFT/Nwpw/libraryps



 initializing nwpw_APC object
 ----------------------------
 nga =   3 ngs =     6
 Gc  =   2.50000
 APC gamma: 0 0.60000
 APC gamma: 1 0.90000
 APC gamma: 2 1.35000
 - not self-consistent

 reading formatted psp filename: ./W.vpp
 generating random psi from scratch
 Warning - Gram-Schmidt being performed on psi2

          ==============  summary of input  ==================

 input psi filename: ./w2-pspw.movecs

 number of processors used: 4
 processor grid           : 4 x1
 parallel mapping         : 1d-slab
 parallel mapping         : balanced

 options:
   boundary conditions  = periodic
   electron spin        = restricted
   exchange-correlation = BEEF (White and Bird) parameterization

 elements involved in the cluster:
     1 :    W   core charge:  6.0  lmax=2
           comment : Troullier-Martins pseudopotential
           pseudopotential type            :   0
           highest angular component       :   2
           local potential used            :   0
           number of non-local projections :   8
           semicore corrections included   :  1.800 (radius)  4.477 (charge)
           cutoff =    2.389   3.185   2.244

 atom composition:
   W : 2

 initial ion positions (au):
   1 W	(    0.00000    0.00000   -2.17318 ) - atomic mass = 183.951
   2 W	(    0.00000    0.00000    2.17318 ) - atomic mass = 183.951
   G.C.	(    0.00000    0.00000    0.00000 )
 C.O.M.	(    0.00000    0.00000    0.00000 )

 number of electrons: spin up=     6 (   6 per task) down=     6 (   6 per task)

 supercell:
      volume =    8000.00
      lattice:    a1=<   20.000    0.000    0.000 >
                  a2=<    0.000   20.000    0.000 >
                  a3=<    0.000    0.000   20.000 >
      reciprocal: b1=<    0.314    0.000    0.000 >
                  b2=<    0.000    0.314    0.000 >
                  b3=<    0.000    0.000    0.314 >
      lattice:    a=      20.000 b=     20.000 c=      20.000
                  alpha=  90.000 beta=  90.000 gamma=  90.000
      density cutoff=  40.000 fft=   60 x   60 x   60  (   48485 waves    12121 per task)
      wavefnc cutoff=  20.000 fft=   60 x   60 x   60  (   17133 waves     4283 per task)

 ewald parameters:
      energy cutoff=  40.000 fft=   60 x   60 x   60  (   48485 waves    12122 per task)
      Ewald summation: cut radius=    6.366 and   1
                       Mandelung Wigner-Seitz=   1.76011888 (alpha=  2.83729748 rs=12.40700982)

 technical parameters:
      time step=        5.80  ficticious mass=  400000.00
      tolerance=   1.000e-07 (energy)    1.000e-07 (density)    1.000e-04 (ion)
      max iterations =       1000 (   10 inner   100 outer)
      minimizer = Grassmann conjugate gradient



     ============ Grassmann conjugate gradient iteration ============
          >>> iteration started at Sat Dec 11 21:37:26 2021  <<<
     iter.                 Energy          DeltaE        DeltaRho
     ----------------------------------------------------------------
        - 10 steepest descent iterations performed
        10      -1.245392026983e+01   -8.933668e-01    5.785399e-04
        - 10 steepest descent iterations performed
        20      -1.807008263738e+01   -8.183934e-02    7.507453e-06
        - 10 steepest descent iterations performed
        30      -1.886354215617e+01   -1.867658e-02    1.174123e-05
        - 10 steepest descent iterations performed
        40      -1.899174285750e+01   -1.499020e-02    6.393621e-05
        - 10 steepest descent iterations performed
        50      -1.917436042361e+01   -5.114398e-03    9.831002e-07
        60      -1.927603165545e+01   -4.191019e-03    2.429206e-05
        70      -1.934484360561e+01   -2.056803e-03    6.607760e-06
        80      -1.937850862199e+01   -1.935696e-03    3.049967e-05
        90      -1.938426079709e+01   -2.523059e-04    1.030353e-06
       100      -1.938555620911e+01   -8.054072e-05    5.233555e-07
       110      -1.938607962555e+01   -3.829212e-05    2.171797e-07
       120      -1.938634874115e+01   -2.079205e-05    1.084291e-07
       130      -1.938650607456e+01   -1.364562e-05    5.682140e-08
       140      -1.938661418874e+01   -9.904439e-06    3.872361e-08
       150      -1.938669752233e+01   -8.515435e-06    3.494277e-08
       160      -1.938677238714e+01   -7.899617e-06    3.540230e-08
       170      -1.938684510741e+01   -8.285501e-06    4.697876e-08
        - 10 steepest descent iterations performed
       180      -1.938699350692e+01   -6.847641e-05    3.482651e-06
        - 10 steepest descent iterations performed
       190      -1.938717532734e+01   -9.580587e-06    5.406413e-08
       200      -1.938731163075e+01   -1.750145e-05    1.406716e-07
        - 10 steepest descent iterations performed
       210      -1.938766478997e+01   -1.848137e-04    1.466250e-05
        - 10 steepest descent iterations performed
       220      -1.938817153390e+01   -2.379111e-05    1.260890e-07
       230      -1.938855103074e+01   -4.140231e-05    3.264751e-07
        - 10 steepest descent iterations performed
       240      -1.938932116389e+01   -2.782129e-04    1.267459e-05
        - 10 steepest descent iterations performed
       250      -1.939039102745e+01   -6.382677e-05    4.383698e-07
       260      -1.939133221831e+01   -8.492176e-05    5.811988e-07
        - 10 steepest descent iterations performed
       270      -1.939319125253e+01   -7.143648e-04    4.233226e-05
        - 10 steepest descent iterations performed
       280      -1.939540696710e+01   -9.577615e-05    5.285681e-07
       290      -1.939691865235e+01   -1.272823e-04    9.957555e-07
        - 10 steepest descent iterations performed
       300      -1.939882512401e+01   -4.107051e-04    1.096077e-05
        - 10 steepest descent iterations performed
       310      -1.940047763795e+01   -8.157201e-05    4.542819e-07
       320      -1.940157261412e+01   -7.604079e-05    3.923439e-07
       330      -1.940248181333e+01   -5.383331e-05    2.394074e-07
       340      -1.940322404289e+01   -4.950076e-05    2.049163e-07
       350      -1.940381991611e+01   -3.433330e-05    1.232326e-07
       360      -1.940429773591e+01   -3.153618e-05    1.062861e-07
       370      -1.940467822734e+01   -2.188991e-05    6.499503e-08
       380      -1.940498396701e+01   -2.032339e-05    5.738865e-08
       390      -1.940522897971e+01   -1.426958e-05    3.572022e-08
       400      -1.940542815761e+01   -1.343789e-05    3.278045e-08
       410      -1.940558947860e+01   -9.517934e-06    2.086061e-08
       420      -1.940572242499e+01   -9.142691e-06    1.983736e-08
       430      -1.940583138779e+01   -6.536793e-06    1.286346e-08
       440      -1.940592254988e+01   -6.404798e-06    1.271395e-08
       450      -1.940599812531e+01   -4.607887e-06    8.343663e-09
       460      -1.940606222984e+01   -4.600482e-06    8.534901e-09
       470      -1.940611586915e+01   -3.321221e-06    5.634857e-09
       480      -1.940616191544e+01   -3.372160e-06    5.942654e-09
        - 10 steepest descent iterations performed
       490      -1.940621987945e+01   -1.872682e-05    2.314575e-07
        - 10 steepest descent iterations performed
       500      -1.940626545288e+01   -1.732387e-06    2.208313e-09
       510      -1.940628998318e+01   -1.765801e-06    3.234825e-09
        - 10 steepest descent iterations performed
       520      -1.940631870938e+01   -6.772758e-06    5.173682e-08
        - 10 steepest descent iterations performed
       530      -1.940634173834e+01   -1.061064e-06    1.628664e-09
       540      -1.940635656370e+01   -1.020112e-06    1.432302e-09
       550      -1.940636921159e+01   -8.100749e-07    1.233434e-09
       560      -1.940638053274e+01   -7.889757e-07    1.094588e-09
       570      -1.940639019274e+01   -6.249786e-07    9.425830e-10
       580      -1.940639890942e+01   -6.151006e-07    8.459941e-10
       590      -1.940640633646e+01   -4.852896e-07    7.269102e-10
       600      -1.940641308813e+01   -4.811830e-07    6.548887e-10
       610      -1.940641883443e+01   -3.816117e-07    5.581760e-10
       620      -1.940642409397e+01   -3.809830e-07    5.235669e-10
       630      -1.940642855799e+01   -2.980238e-07    4.280020e-10
       640      -1.940643267567e+01   -3.035407e-07    4.226349e-10
        - 10 steepest descent iterations performed
       650      -1.940643784408e+01   -1.647156e-06    1.616945e-08
        - 10 steepest descent iterations performed
       660      -1.940644227204e+01   -1.805499e-07    2.040816e-10
       670      -1.940644466209e+01   -1.736989e-07    3.046672e-10
       680      -1.940644693826e+01   -1.420316e-07    1.618246e-10
       690      -1.940644880630e+01   -1.352906e-07    2.241373e-10
       700      -1.940645060036e+01   -1.238282e-07    1.583741e-10
       710      -1.940645206054e+01   -9.820953e-08    1.585936e-10
     *** tolerance ok. iteration terminated
          >>> iteration ended at   Sat Dec 11 21:42:19 2021  <<<

     ==============  energy results (Molecule object)  ==============


 number of electrons: spin up=     6.00000  down=     6.00000 (real space)

 total     energy    :   -1.9406452061e+01 (   -9.70323e+00 /ion)
 total orbital energy:   -1.7122461818e+00 (   -2.85374e-01 /electron)
 hartree energy      :    1.2594417449e+01 (    2.09907e+00 /electron)
 exc-corr energy     :   -8.9802734407e+00 (   -1.49671e+00 /electron)
 ion-ion energy      :   -1.7408366198e+00 (   -8.70418e-01 /ion)

 kinetic (planewave) :    1.4597219343e+01 (    2.43287e+00 /electron)
 V_local (planewave) :   -2.0582587603e+01 (   -3.43043e+00 /electron)
 V_nl    (planewave) :   -1.5294391189e+01 (   -2.54907e+00 /electron)
 V_Coul  (planewave) :    2.5188834897e+01 (    4.19814e+00 /electron)
 V_xc    (planewave) :   -5.6213216304e+00 (   -9.36887e-01 /electron)
 Viral Coefficient   :   -1.1172994761e+00

 orbital energy:
    -1.0904280e-01 (  -2.967eV)
    -1.0907445e-01 (  -2.968eV)
    -1.4344319e-01 (  -3.903eV)
    -1.4998463e-01 (  -4.081eV)
    -1.4998519e-01 (  -4.081eV)
    -1.9459283e-01 (  -5.295eV)


 Ion Forces (au):
   1 W	(    0.00262   -0.02882    0.05993 )
   2 W	(   -0.00273    0.02878   -0.06005 )



*************************************************************
**                                                         **
**          PSPW Atomic Point Charge (APC) Analysis        **
**                                                         **
**   Point charge analysis based on paper by P.E. Blochl   **
**         (J. Chem. Phys. vol 103, page 7422, 1995)       **
**                                                         **
*************************************************************

 nwpw_APC object
 ---------------
 nga =   3 ngs =     6
 Gc  =   2.50000
 APC gamma: 0 0.60000
 APC gamma: 1 0.90000
 APC gamma: 2 1.35000
 - not self-consistent


 charge analysis on each atom
 ----------------------------

      no  atom        Qelc        Qion      Qtotal
   -----  ----     -------     -------     -------
       1    W       -6.000       6.000       0.000
       2    W       -6.000       6.000      -0.000
       Total Q     -12.000      12.000       0.000


 gaussian coefficients of model density
 --------------------------------------

      no  atom     g=0.000     g=0.600     g=0.900     g=1.350
   -----  ----     -------     -------     -------     -------
       1     W       6.000     -17.033      34.586     -23.552
       2     W       6.000     -17.036      34.590     -23.554


 output psi to filename: ./w2-pspw.movecs

 -----------------
 cputime in seconds
 prologue    : 1.04876e-01
 main loop   : 2.93376e+02
 epilogue    : 7.56090e-02
 total       : 2.93556e+02
 cputime/step: 1.17964e-01 ( 2487 evaluations, 10 linesearches)

 Time spent doing      total        step             percent
 total time            2.939356e+02 1.181888e-01     100.00%
 total FFT time        1.644260e+02 6.611421e-02      55.94%
 lagrange multipliers  4.516679e-01 1.816116e-04       0.15%
 local potentials      2.181282e-03 8.770736e-07       0.00%
 non-local potentials  1.425822e+01 5.733100e-03       4.85%
 ffm_dgemm             2.968983e+00 1.193801e-03       1.01%
 fmf_dgemm             9.132324e+00 3.672024e-03       3.11%
 m_diagonalize         1.366800e-02 5.495776e-06       0.00%
 mmm_multiply          9.537384e-03 3.834895e-06       0.00%

 >>> job completed at     Sat Dec 11 21:42:19 2021 <<<
output rtdbstr={"current_task":"task pspw gradient","dbname":"w2-pspw","foundtask":false,"geometries":{"geometry":{"autosym":0,"autoz":0,"center":1,"charges":[74.0,74.0],"conv":1.88972598858,"coords":[0.0,0.0,-2.173184886867,0.0,0.0,2.173184886867],"masses":[183.951,183.951],"nion":2,"symbols":["W","W"],"velocities":[0.0,0.0,0.0,0.0,0.0,0.0]}},"geometry":null,"nwinput_cur":24,"nwinput_lines":["Title \"Hellow Dolly\"","","memory 1900 mb","start w2-pspw","echo","","","","geometry noautosym noautoz center","W 0 0 0    ","W 0 0 2.3  ","end","","nwpw","   simulation_cell","     SC 20.0","   end","   cutoff 20.0","   xc beef","end","","task pspw gradient","",""],"nwinput_nlines":24,"nwpw":{"apc":{"on":true,"q":[5.6291123200935544e-05,-5.6291119108209386e-05]},"cutoff":[20.0,40.0],"simulation_cell":{"unita":[20.0,0.0,0.0,0.0,20.0,0.0,0.0,0.0,20.0]},"xc":"beef"},"permanent_dir":".","psp_library_dir":"","pspw":{"eigenvalues":[-0.10904279623518386,-0.1090744534904712,-0.1434431879678982,-0.14998463147416577,-0.14998518967766503,-0.1945928320749004],"energies":[-19.406452060536864,-1.712246181840569,12.59441744862369,-8.980273440714464,-1.740836619784234,14.597219343007769,-20.5825876026156,-15.294391189054037,25.18883489724738,-5.62132163042609,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],"energy":-19.406452060536864,"fion":[0.0026249595837127835,-0.028819370135758422,0.05993387407025376,-0.0027345683434777395,0.028778125118681035,-0.060047567556791644],"gradient":[-0.0026249595837127835,0.028819370135758422,-0.05993387407025376,0.0027345683434777395,-0.028778125118681035,0.060047567556791644]},"scratch_dir":".","title":"Hellow Dolly"}
[Erics-MacBook-Pro:PWDFT/QA/LibExample] bylaska% 

```
