
# Result

环上序列的截位还原问题

## cpuinfo

```sh
 % cat /proc/cpuinfo | grep 'model name'
model name      : Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz
```

## category 1

### level 1

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 18 2 --category 1 --level 1 --verbose 1 --block-size 19
SEED: 3409482587277765135

solution: (257122259, 561754033, 340567466, 370127561, 1603890151, 789710361, 438276282, 1205614745, 929435387, 273101854, 1330188497, 1927005651, 70974738, 512638222, 655376420, 1799812840)

================
CPU     371%
user    2.130
system  1.530
total   0.986
```

### level 2

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 19 3 --category 1 --level 2 --verbose 1 --block-size 20
SEED: 741431373315218258

solution: (720759561, 571519362, 1471239584, 2146374199, 1943406434, 1021387208, 911847040, 972284090, 1269846527, 2125700056, 582856260, 1326807684, 820744516, 83875554, 1033926054, 974403091)

================
CPU     337%
user    1.903
system  1.347
total   0.963
```

### level 3

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 30 14 --category 1 --level 3 --verbose 1 --block-size 20
SEED: 6296888057838770351

solution: (2056319062, 1175325906, 34344908, 300877541, 871395921, 953051611, 1276817066, 429446330, 716236050, 1498148665, 419137803, 983513800, 877501144, 162784581, 615817441, 1532450981)

================
CPU     362%
user    2.133
system  1.565
total   1.020
```

### level 4

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 51 21 --category 1 --level 4 --verbose 1 --block-size 40
SEED: 14464219551167228876

solution: (146123803, 1660954690, 553686861, 1592631770, 2039784960, 874444650, 1462760700, 1629573947, 927239148, 2020986341, 2134682761, 1440980008, 214415113, 823589071, 1178840115, 237668181)

================
CPU     311%
user    2.567
system  1.425
total   1.284
```

### level 5

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 75 24 --category 1 --level 5 --verbose 1 --block-size 40
SEED: 2380040916040384452

solution: (1683538635, 247852287, 1110454234, 2134379965, 524624671, 1135659173, 611030817, 917303344, 67431860, 844532212, 2121384016, 1630029172, 1311537205, 1289944654, 358213366, 2024311471)

================
CPU     194%
user    8.672
system  1.823
total   5.389
```

### level 6

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 120 26 --category 1 --level 6 --verbose 2 --block-size 30 --threads 10 --sieve
SEED: 17559565087081154795

E|v|         = 23619293351.5
E|b[0]|      = 20806089090.1
E|v|/E|b[0]| = 1.135

gh = 432893343225474777088.000000, goal_r0/gh = 1.289033, r0/gh = 3.069376
  31: ↑ 31 ↓  1  T:   0.61309s, TT:   0.61312s, r0:1.03855e+21 r0/gh:   2.39910
  34: ↑ 34 ↓  1  T:   0.67278s, TT:   1.29709s, r0:1.03855e+21 r0/gh:   2.39910
  37: ↑ 37 ↓  1  T:   1.03919s, TT:   2.34736s, r0:1.03855e+21 r0/gh:   2.39910
  40: ↑ 40 ↓  6  T:   1.20862s, TT:   3.56684s, r0:1.03855e+21 r0/gh:   2.39910
  43: ↑ 43 ↓  5  T:   1.70155s, TT:   5.27941s, r0:1.03855e+21 r0/gh:   2.39910
  46: ↑ 46 ↓  8  T:   2.22245s, TT:   7.51265s, r0:1.03855e+21 r0/gh:   2.39910
  49: ↑ 49 ↓ 18  T:   2.69158s, TT:  10.21570s, r0:1.03855e+21 r0/gh:   2.39910
  52: ↑ 52 ↓ 12  T:   4.12785s, TT:  14.35474s, r0:1.03855e+21 r0/gh:   2.39910
  55: ↑ 55 ↓ 14  T:   4.34259s, TT:  18.70839s, r0:9.68850e+20 r0/gh:   2.23808
  58: ↑ 58 ↓ 22  T:   4.46029s, TT:  23.18006s, r0:9.37393e+20 r0/gh:   2.16541
  61: ↑ 61 ↓ 29  T:   4.61750s, TT:  27.80868s, r0:9.37393e+20 r0/gh:   2.16541
  64: ↑ 64 ↓ 35  T:   5.35344s, TT:  33.17343s, r0:8.38208e+20 r0/gh:   1.93629
  67: ↑ 67 ↓ 28  T:   8.15652s, TT:  41.34195s, r0:8.38208e+20 r0/gh:   1.93629
  70: ↑ 70 ↓ 43  T:  10.30037s, TT:  51.65352s, r0:8.03424e+20 r0/gh:   1.85594
  73: ↑ 73 ↓ 45  T:  16.68387s, TT:  68.34896s, r0:6.19472e+20 r0/gh:   1.43100
  76: ↑ 76 ↓ 50  T:  28.59283s, TT:  96.95415s, r0:6.19472e+20 r0/gh:   1.43100
  79: ↑ 79 ↓ 78  T:  29.21822s, TT: 126.18474s, r0:1.95242e+20 r0/gh:   0.45102
svp: norm 118207.1 ,hf 0.00001
('asvp', 121) :: threads: 10, cputime: 859.8936s, walltime: 126.2022s, flast: 42.00, |db|: 2^15.52
solution: (512741665, 1096178369, 808807049, 608491186, 420891056, 1682771835, 358966452, 55989687, 890238631, 2137448551, 2058494244, 38743896, 96170410, 49379854, 1551639435, 1878181314)

================
CPU     656%
user    14:28.15
system  3.877
total   2:12.83
```

### level 7

for the multi-threaded version, the solution cannot be found after 10 threads run for about 12 hours

```sh
time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 150 27 --category 1 --level 7 --verbose 2 --block-size 30 --threads 10 --sieve
```

while sieving over GPU, a Tesla T4, we can get the results in 2 hours

POC at colab: <https://colab.research.google.com/drive/1MwqKIxTzkqJMMD8vx3MhHY0l-8g976rm?usp=sharing>

log: <https://drive.google.com/file/d/1Y1B7usFfgONYfj78Hnt5POUDI5xMHEXT/view?usp=sharing>

## category 2

### level 1

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 2 30 8 17 --category 2 --level 1 --verbose 1 --block-size 20 --check
SEED: 12928939533349146289

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: 7 0
checking z_i: maybe
checking c_i: True

================
CPU     102%
user    2.144
system  0.261
total   2.354
```

```json
{
    "modulus": 2147483647,
    "zbits": 17,
    "coefficients": [1596998372, 913674193],
    "initial_state": [25583676, 1935662022]
}
```

### level 2

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 2 60 15 23 --category 2 --level 2 --verbose 1 --block-size 20 --check
SEED: 608460742630951138

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: -26204 -811
checking z_i: maybe
checking c_i: True

================
CPU     101%
user    2.714
system  0.237
total   2.902
```

```json
{
    "modulus": 2147483647,
    "zbits": 23,
    "coefficients": [423368878, 1375517413],
    "initial_state": [1968687461, 159779378]
}
```

### level 3

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 3 68 17 21 --category 2 --level 3 --verbose 1 --block-size 20 --check
SEED: 15922487348339387766

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: -1610 -50
checking z_i: maybe
checking c_i: True

================
CPU     101%
user    2.914
system  0.296
total   3.160
```

```json
{
    "modulus": 2147483647,
    "zbits": 21,
    "coefficients": [233454232, 712694596, 1250324919],
    "initial_state": [1968687461, 159779378, 933973255]
}
```

### level 4

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 4 95 25 21 --category 2 --level 4 --verbose 1 --block-size 30 --check 
SEED: 18135446848677888257

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: 1540 9
checking z_i: maybe
checking c_i: True

================
CPU     100%
user    7.603
system  0.266
total   7.828
```

```json
{
    "modulus": 2147483647,
    "zbits": 21,
    "coefficients": [1724886998, 287764120, 669496309, 29431304],
    "initial_state": [23214526, 63888791, 632526187, 515165230]
}
```

### level 5

```sh
 % time sage recover_coefficients__kernel.sage 2113941029 5 85 23 18 --category 2 --level 5 --verbose 1 --block-size 30 --check
SEED: 13512413672535046818

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: -24 0
checking z_i: maybe
checking c_i: True

================
CPU     100%
user    6.491
system  0.256
total   6.697
```

```json
{
    "modulus": 2113941029,
    "zbits": 18,
    "coefficients": [241592126, 1225700761, 270381722, 1809937814, 545364186],
    "initial_state": [1265879737, 1938499214, 100295411, 164486483, 782938]
}
```

### level 6

```sh
 % time sage recover_coefficients__kernel.sage 2140900439 8 90 20 11 --category 2 --level 6 --verbose 1 --block-size 20 --check
SEED: 2774189050861033585

expect_vectors: 6

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    4.257
system  0.153
total   4.403
```

```json
{
    "modulus": 2140900439,
    "zbits": 11,
    "coefficients": [1468898684, 429201201, 1438911747, 1343646518, 197478154, 1760674261, 1954064960, 1521596057],
    "initial_state": [2100894922, 1278709211, 1882769045, 822236456, 503871792, 639865048, 740172666, 1225091197]
}
```

### level 7

```sh
 % time sage recover_coefficients__kernel.sage 2086596509 10 110 26 11 --category 2 --level 7 --verbose 1 --block-size 20 --check
SEED: 4768417682194502211

expect_vectors: 6

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    9.136
system  0.159
total   9.289
```

```json
{
    "modulus": 2086596509,
    "zbits": 11,
    "coefficients": [1111873952, 1210156476, 822665602, 1221543376, 50425289, 1211476215, 1888560914, 1679919063, 1715756131, 141246785],
    "initial_state": [1800774673, 1719445571, 1277627869, 633482595, 1079260842, 2060264416, 573852671, 1245793554, 141816054, 891093089]
}
```

### level 8

```sh
 % time sage recover_coefficients__kernel.sage 2123058169 12 110 28 8 --category 2 --level 8 --verbose 1 --block-size 20 --check
SEED: 12270745752118493764

expect_vectors: 5

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    4.516
system  0.145
total   4.652
```

```json
{
    "modulus": 2123058169,
    "zbits": 8,
    "coefficients": [1380518532, 572802739, 397998604, 1517367287, 1838517468, 1894991963, 186761507, 1691926163, 917271042, 1840254211, 1520485994, 1544456793],
    "initial_state": [1407121795, 1964561578, 1522896740, 706799210, 1022864344, 1609032902, 945477963, 1469966024, 232347115, 2025154173, 653675533, 603475034]
}
```

### level 9

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 14 128 32 8 --category 2 --level 9 --verbose 1 --block-size 32 --check
SEED: 8421578719646212792

expect_vectors: 5

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    9.820
system  0.190
total   10.004
```

```json
{
    "modulus": 2147483647,
    "zbits": 8,
    "coefficients": [755735009, 435105367, 1987422269, 141113323, 1831273687, 150474978, 1521781010, 88703098, 2128502238, 1314935750, 1897202874, 765777736, 1257457888, 851182418],
    "initial_state": [2043273873, 1844618476, 93566804, 974010131, 1813982016, 210537594, 1382311654, 397516522, 1569721206, 1648699957, 971748769, 1513210834, 1522349751, 1075882738]
}
```

## category 3

### level 1

```sh
 % time sage recover_modulus__kernel.sage 31 16 140 30 5 --category 3 --level 1 --verbose 1 --block-size 20 --check                                           
SEED: 10019631545210073147
FLAGS = 01

find the kernel (rank 110)

expect_vectors: 6

find the kernel (rank 2)

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     278%
user    1:25.33
system  0.483
total   30.772
```

```json
{
    "modulus": 2123847813,
    "zbits": 5,
    "coefficients": [2018484748, 1845189071, 1463275556, 1602150465, 194422107, 1025586312, 1724843525, 410393693, 106087782, 229852172, 1293380914, 235543842, 1642599451, 250756393, 239416449, 1593118903],
    "initial_state": [1422595032, 698794309, 2264898, 2123217555, 1515919515, 2048279701, 1227019859, 1625549320, 475257352, 1682624639, 1669847210, 510649846, 272012336, 1608958057, 1318317428, 1804116296]
}
```

### level 2

```sh
 % time sage recover_modulus__kernel.sage 31 16 190 40 10 --category 3 --level 2 --verbose 1 --block-size 30 --check            
SEED: 2732523400034528389
FLAGS = 01

find the kernel (rank 150)

expect_vectors: 6

find the kernel (rank 2)

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     318%
user    3:57.23
system  0.699
total   1:14.81
```

```json
{
    "modulus": 2146390813,
    "zbits": 10,
    "coefficients": [1709517653, 1473447434, 146866621, 1301261246, 1098198483, 1586144939, 880631859, 1190804449, 419206704, 377180855, 997067781, 668707083, 1991423106, 1200018419, 1124879071, 2081342702],
    "initial_state": [963105734, 789585151, 1238195227, 2028522939, 1124205863, 1618865668, 452174891, 77673612, 46901163, 1112184517, 2114462053, 259959215, 1976235589, 1517149832, 147549104, 73665604]
}
```

### level 3

```sh
 % time sage recover_modulus__kernel.sage 31 16 265 70 14 --category 3 --level 3 --verbose 2 --block-size 30 --threads 10 --sieve --timeout 300 --max-dim 80 --check
SEED: 11316827758390952469
FLAGS = 01

find the kernel (rank 195)

Loaded file 'svpchallenge-195.txt'
gh = 165086.056394, goal_r0/gh = 0.000000, r0/gh = 0.909447
'threads': 10,       :: n: 195, cputime  0.0090s, walltime:  0.0090s, flast: -1.00, |db|: 2^0.00
expect_vectors: 5

find the kernel (rank 2)

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     815%
user    43:59.48
system  3.875
total   5:23.99
```

```json
{
    "modulus": 2147385873,
    "zbits": 14,
    "coefficients": [1433408125, 1187896517, 1493677871, 674581842, 873135996, 1326711093, 1726922573, 833159114, 897246300, 1821464147, 1671306427, 570111420, 268554940, 1034987293, 1263557877, 894951295],
    "initial_state": [1976346479, 1550017240, 380413913, 160326487, 379904755, 622685282, 377075770, 1527069890, 946605673, 2096385369, 1160042120, 619378670, 50681874, 1753705306, 1683069046, 1512092243]
}
```
