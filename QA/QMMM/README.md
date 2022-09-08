## Generate pwdft shared library ##
First go to the PWDFT home directory, e.g.,

https://jsfiddle.net/8ndx694g/


<img src="https://render.githubusercontent.com/render/math?math=\large E = E_{QM} %2B E_{MM} %2B E_{QM/MM}">

<img src="https://render.githubusercontent.com/render/math?math=\large %5Cbegin%7Balign*%7D%0AE_%7BMM%7D%20%26%3D%20%5Cfrac%7B1%7D%7B2%7D%20%5Csum_%7Bi%20%5Cneq%20j%7D%20%5Cfrac%7Bq_i%20q_j%7D%7Br_%7Bij%7D%7D%20%5C%5C%0A%20%20%20%20%20%20%20%26%20%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Csum_%7Bi%20%5Cneq%20j%7D%20%5C%7B%20%5Cfrac%7BA_%7Bij%7D%7D%7Br_%7Bij%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7Bij%7D%7D%7Br_%7Bij%7D%5E6%7D%20%5C%7D%0A%5Cend%7Balign*%7D">



<img src="https://render.githubusercontent.com/render/math?math=\large %5Cbegin%7Balign*%7D%0AE_%7BQM%2FMM%7D%20%26%3D%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BQ_I%20q_i%7D%7Br_%7BIi%7D%7D%20%5C%7D%5C%5C%0A%26%2B%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BA_%7BIi%7D%7D%7Br_%7BIi%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7BIi%7D%7D%7Br_%7BIi%7D%5E6%7D%20%5C%7D%0A%5Cend%7Balign*%7D">

<img src="https://render.githubusercontent.com/render/math?math=\large %5Cbegin%7Balign*%7D%0AE_%7BQM%2FMM%7D%20%26%3D%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BQ_I%20q_i%7D%7Br_%7BIi%7D%7D%20%5C%7D%5C%5C%0A%26%2B%20%5Csum_%7BI%3DQM%7D%20%5Csum_%7Bi%3DMM%7D%20%5C%7B%20%5Cfrac%7BA_%7BIi%7D%7D%7Br_%7BIi%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7BIi%7D%7D%7Br_%7BIi%7D%5E6%7D%20%5C%7D%5C%5C%0A%26%3D%20%5Csum_%7BI%3DQM%7D%20Q_I%20U_I%20%5C%5C%0A%26%2B%20e_%7BLJ%7D(I%2Ci)%0A%5Cend%7Balign*%7D">


<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0AE_%7BQM%2FMM%7D%20%26%3D%20%5Csum_%7BIi%7D%20%5Cfrac%7BQ_I%20q_i%7D%7Br_%7BIi%7D%7D%5C%5C%0A%26%2B%20%5Csum_%7BIi%7D%20%5C%7B%20%5Cfrac%7BA_%7BIi%7D%7D%7Br_%7BIi%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7BIi%7D%7D%7Br_%7BIi%7D%5E6%7D%20%5C%7D%0A%5Cend%7Balign*%7D">

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0AE_%7BAPC%7D%20%26%3D%20%5Csum_I%20Q_I%20U_I%0A%5Cend%7Balign*%7D">

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0AU_%7BI%7D%20%26%3D%20%5Csum_%7Bi%7D%20%5Cfrac%7Bq_i%7D%7Br_%7BIi%7D%7D%0A%5Cend%7Balign*%7D">

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0AE_%7BQQ%7D%20%26%3D%20%5Cfrac%7B1%7D%7B2%7D%20%5Csum_%7BI%20%5Cneq%20J%7D%20%5Cfrac%7BQ_I%20Q_J%7D%7Br_%7BIJ%7D%7D%5C%5C%0A%26%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Csum_%7BI%20%5Cneq%20J%7D%20%5C%7B%20%5Cfrac%7BA_%7BIJ%7D%7D%7Br_%7BIJ%7D%5E%7B12%7D%7D%20-%20%5Cfrac%7BB_%7BIJ%7D%7D%7Br_%7BIJ%7D%5E6%7D%20%5C%7D%0A%5Cend%7Balign*%7D">
